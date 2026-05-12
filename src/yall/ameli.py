##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides the interface to the AMELI repository of tensor
# matrices on Zenodo.
#
##########################################################################
import json
from datetime import datetime, timedelta
from functools import total_ordering
import io
import logging
import math
from pathlib import Path
import requests
import tempfile
import zipfile

import numpy as np
import h5py
from requests import Timeout

logger = logging.getLogger("yall.ameli")

AMELI_PATH = Path(__file__).resolve().parent / "ameli"
MATRIX_PATH = Path(__file__).resolve().parent / "matrix"

ZENODO_API = "https://zenodo.org/api"
ZENODO_RECORD = f"{ZENODO_API}/records"
ZENODO_REFRESH_HOURS = 12
RECORDS = {
    "f1": 19130697,
    "f2": 19139480,
    "f3": 19144764,
    "f4": 19154321,
    "f5": 19154326,
    "f6": 19158643,
    "f7": 19158647,
    "f8": 19158658,
    "f9": 19158660,
    "f10": 19158667,
    "f11": 19158671,
    "f12": 19158675,
    "f13": 19158677,
}


##########################################################################
# Update of AMELI matrices
##########################################################################

@total_ordering
class Version:
    def __init__(self, version_str):
        self.original_str = version_str
        cleaned = version_str.lower().lstrip('v')
        parts = cleaned.split('.')

        try:
            self.major = int(parts[0]) if len(parts) > 0 else 0
            self.minor = int(parts[1]) if len(parts) > 1 else 0
            self.patch = int(parts[2]) if len(parts) > 2 else 0
        except (ValueError, IndexError) as e:
            raise ValueError(f"Invalid version format: {version_str}") from e

    def _as_tuple(self):
        return (self.major, self.minor, self.patch)

    def __eq__(self, other):
        if not isinstance(other, Version):
            return NotImplemented
        return self._as_tuple() == other._as_tuple()

    def __lt__(self, other):
        if not isinstance(other, Version):
            return NotImplemented
        return self._as_tuple() < other._as_tuple()

    def __repr__(self):
        return f"Version({self.major}.{self.minor}.{self.patch})"

    def __str__(self):
        return f"{self.major}.{self.minor}.{self.patch}"


def get_local_version(config_path):
    """ Read the first non-empty line from file "VERSION" in the configuration folder and return it as Version
    and datetime objects. Return None as version if the file is missing yet. """

    version_file = Path(config_path) / "VERSION"

    # Return None if the file doesn't exist
    if not version_file.exists():
        return None, None

    # Read file and split into lines
    content = version_file.read_text(encoding="utf-8").strip()
    assert content, "Empty VERSION file!"
    assert "\n" not in content, "VERSION file contains multiple lines!"

    # This will raise an Exception internally if the content is invalid
    version, timestamp = content.split(":", 1)
    version = Version(version.strip())
    timestamp = datetime.fromisoformat(timestamp.strip())
    return version, timestamp


def get_zenodo_version(concept_id):
    """ Return version string of given Zenodo concept ID or None in case of an error. """

    url = f"{ZENODO_RECORD}/{concept_id}"

    # Try to read Zenodo record, but ignore it if it not accessible yet
    try:
        response = requests.get(url, timeout=3)
        response.raise_for_status()
    except Timeout:
        logger.info(f"Zenodo timeout on record {concept_id}")
        return None
    except Exception as e:
        logger.warning(f"Warning: Access to record {concept_id} failed ({e})")
        return None

    data = response.json()
    version = data["metadata"]["version"]
    return Version(version)


def download_files(concept_id, filenames, config_path):
    """ Download multiple ZIP files from the given Zenodo record to a temporary folder. Extract them into the
     config_path if all downloads succeed, then clean up. """

    # Sanity check
    config_path = Path(config_path)
    assert config_path.exists()

    # Get record metadata
    url = f"{ZENODO_RECORD}/{concept_id}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    # Record version
    version = Version(data["metadata"]["version"])

    # Metadata of all files in the record
    files_in_record = data.get('files', [])

    # Use a temporary directory context manager. Folder and contents are deleted when the "with" block is left.
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Initialise list of temporary local file paths
        downloaded_paths = []

        # Download all files to the temporary folder
        for filename in filenames:
            assert filename.endswith(".zip")

            # Metadata of this file
            file_data = next((f for f in files_in_record if f['key'] == filename), None)
            if not file_data:
                message = f"Error: File '{filename}' not found in Zenodo record {concept_id}"
                logger.error(message)
                raise Exception(message)

            # Download the file as stream to stay memory-efficient for large files
            local_tmp_file = tmp_path / filename
            with requests.get(file_data["links"]["self"], stream=True) as r:
                r.raise_for_status()
                with open(local_tmp_file, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)

            # Append temporary file path
            downloaded_paths.append(local_tmp_file)

        # Extract files
        for zip_path in downloaded_paths:
            with zipfile.ZipFile(zip_path) as zip_ref:
                zip_ref.extractall(config_path)

            # Remove hashes file
            hashes = config_path / "hashes.json"
            if hashes.exists():
                hashes.unlink()

    # Update the VERSION file
    update_version(version, config_path)
    logger.debug(f"Zenodo record {concept_id} updated to {version} in {config_path}")


def update_version(version, path):
    """ Update VERSION file with current timestamp. """

    version_file = path / "VERSION"
    timestamp = datetime.now().isoformat(timespec='seconds')
    content = f"{version}: {timestamp}"
    version_file.write_text(content, encoding="utf-8")
    logger.debug(f"Zenodo version {version} timestamp updated in {path}")


def remove_vault(config):
    path = MATRIX_PATH / config
    if path.exists():
        for file in path.iterdir():
            file.unlink()
        path.rmdir()
    logger.debug(f"Local matrix vault for configuration {config} removed")

def update(config, force=False):
    logger.debug(f"Updating configuration {config} with force={force}")

    path = AMELI_PATH / config
    if not path.exists():
        path.mkdir(parents=True)

    local_version, timestamp = get_local_version(path)
    if local_version is None:
        logger.debug(f"Found no local version of configuration {config}")
    else:
        logger.debug(f"Found local version {local_version} of configuration {config} ({timestamp})")
    if local_version is not None and datetime.now() - timestamp < timedelta(hours=ZENODO_REFRESH_HOURS) and not force:
        logger.debug(f"Local version {local_version} is fresh")
        return

    concept_id = RECORDS[config]
    zenodo_version = get_zenodo_version(concept_id)
    logger.debug(f"Zenodo record {concept_id} has version {zenodo_version}")

    if zenodo_version is None:
        assert local_version is not None, f"Cannot download data from Zenodo record {concept_id}!"
        return

    if local_version is None or local_version < zenodo_version:
        filenames = ("product.zip", "sljm.zip", "slj.zip")
        download_files(concept_id, filenames, path)
        remove_vault(config)
    else:
        update_version(local_version, path)


##########################################################################
# Decode square roots of rationals
##########################################################################

def decode_scalar(s: int, n: int, d: int) -> float:
    """ Return (-1)^s * sqrt(n/d) for the given sign, numerator and denominator in optimised numerical precision. """

    assert n > 0
    assert d > 0
    assert s in (0, 1)

    # Size of n/d in bits
    current_diff = n.bit_length() - d.bit_length()

    # Mantissa of Python float is 53 bits and 60 bits provides a safety margin
    shift = 60 - current_diff // 2

    # Scale numerator or denominator
    if shift > 0:
        numerator = n << 2 * shift
        denominator = d
    else:
        numerator = n
        denominator = d << -2 * shift

    # Calculate scaled square root
    root_scaled = math.isqrt(numerator // denominator)

    # Scale back and avoid potential exponent overflow
    value = math.ldexp(root_scaled, -shift)
    if s:
        value = -value
    return value


# Vectorised decoder function
v_decode = np.frompyfunc(decode_scalar, 3, 1)


def decode_vector(s, n, d):
    """ Return (-1)^s * sqrt(n/d) for each element of the given lists. """

    assert len(s) == len(n) == len(d)
    return v_decode(s, n, d).astype(float)


def decode_uint_array(meta, name):
    """ Decode an array encoded by 'encode_large()' with given base-name from the dictionary meta and return it as
    list. Each processed item is removed from meta. """

    # Single storage array
    if name in meta:
        array = meta[name][:].astype(object)

    # Combine large-element array stored in multiple parts
    else:
        keys = [key for key in meta if key.startswith(f"{name}_part")]
        keys.sort()
        array = meta[keys[-1]][:].astype(object)
        for key in reversed(keys[:-1]):
            bits = meta[key].dtype.itemsize * 8
            array = (array << bits) | meta[key][:].astype(object)

    # Return restored array as list
    return array


##########################################################################
# Decode AMELI matrices
##########################################################################

def matrix_path(config, state_space, name):
    name = name.replace("/", "_").replace(",", "_")
    path = AMELI_PATH / config / state_space.lower() / f"{name}.zdc"
    if not path.exists():
        message = f"Matrix file '{path}' does not exist!"
        try:
            concept_id = RECORDS[config]
            url = f"{ZENODO_RECORD}/{concept_id}"
            message += f" Check {url}."
        except:
            pass
        logger.error(message)
        raise ValueError(message)
    return path


# def matrix_mtime(config_name, state_space, name):
#     path = matrix_path(config_name, state_space, name)
#     return path.stat().st_mtime


def read_matrix(path, item):
    """ Return a float representation of the given AMELI matrix. """

    with zipfile.ZipFile(path, "r") as z:
        with z.open(item) as f:
            data = io.BytesIO(f.read())

    root = h5py.File(data, "r")

    sign = decode_uint_array(root, "sign")
    numerator = decode_uint_array(root, "numerator")
    denominator = decode_uint_array(root, "denominator")
    values = decode_vector(sign, numerator, denominator)

    is_symmetric = root.attrs["isSymmetric"]
    num_states = root.attrs["numStates"]
    matrix = np.zeros((num_states, num_states), dtype=float)
    for row, col, index in zip(root["rows"], root["columns"], root["elements"]):
        value = values[index]
        matrix[row, col] = value
        if row != col and is_symmetric:
            matrix[col, row] = value
    return matrix


def get_ameli_matrix(name, config, coupling):
    state_space = coupling.name.lower()
    assert state_space in ("slj", "sljm", "product")
    path = matrix_path(config, state_space, name)
    return read_matrix(path, "data/matrix.hdf5")


def read_indices(path, item):
    """ Return a float representation of the given AMELI matrix. """

    with zipfile.ZipFile(path, "r") as z:
        with z.open(item) as f:
            data = io.BytesIO(f.read())
            return np.array(h5py.File(data, "r")["indices"])


def read_json(path, item):
    """ Return a JSON file from a data container. """

    with zipfile.ZipFile(path, "r") as z:
        with z.open(item) as f:
            return json.loads(f.read())


def read_transform(config):
    path = AMELI_PATH / config / "transform.zdc"
    meta = read_json(path, "data/transform.json")
    transform = {
        "electronPool": meta["row_states"]["electronPool"],
        "rowStates": read_indices(path, "data/row_states.hdf5"),
        "tensorChain": meta["col_states"]["tensorChain"],
        "irreducibleRepresentations": meta["col_states"]["irreducibleRepresentations"],
        "colStates": read_indices(path, "data/col_states.hdf5"),
        "transform": read_matrix(path, "data/matrix.hdf5"),
    }
    return transform
