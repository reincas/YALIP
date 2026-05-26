##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module acts as interface to the HDF5 cache of converted matrices.
# It also generates and diagonalises perturbation Hamiltonians.
#
##########################################################################

import h5py
import logging
import numpy as np

from .ameli import MATRIX_PATH, get_ameli_matrix

logger = logging.getLogger("yalip.matrix")

# Translation of alternative names to canonical names
ALT_NAMES = {
    "Lz": "L/0",
    "Sz": "S/0",
    "Jz": "J/0",
    "L2": "LL",
    "S2": "SS",
    "J2": "JJ",
    "C7": "CR",
    "C5": "CR",
    "H3/0": "LL",
    "H3/1": "C2",
    "H3/2": "CR",
    "H5fix": {"H5/0": 1.0, "H5/2": 0.56, "H5/4": 0.38},
    "H6fix": {"H6/2": 1.0, "H6/4": 0.75, "H6/6": 0.50},
}


def normalize_radial(radial):
    """ Convert an alternative parameter set into the standard set of radial integrals "H1"-"H6". """

    # Get keys and initialize the normalised set.
    keys = list(radial.keys())
    new_radial = {}

    # No conversion for standard parameters
    for key in list(keys):
        if key not in ["F_0", "F_2", "F_4", "F_6", "P_2", "P_4", "P_6", "E^0", "E^1", "E^2", "E^3"]:
            new_radial[key] = radial[key]
            keys.remove(key)

    # Convert "F_k" and "P_k" to "H1/k" and "H6/k", respectively
    for key in list(keys):
        if key == "F_0":
            new_radial["H1/0"] = radial[key]
        elif key == "F_2":
            new_radial["H1/2"] = radial[key] * 225
        elif key == "F_4":
            new_radial["H1/4"] = radial[key] * 1089
        elif key == "F_6":
            new_radial["H1/6"] = radial[key] * 184041 / 25
        elif key == "P_2":
            new_radial["H6/2"] = radial[key] * 225
        elif key == "P_4":
            new_radial["H6/4"] = radial[key] * 1089
        elif key == "P_6":
            new_radial["H6/6"] = radial[key] * 184041 / 25
        else:
            continue
        keys.remove(key)

    # The "E^i" parameters need a linear transformation to "H1/k". "E^0" or "H1/0" can be used as an alternative
    # to "base" to shift the whole energy level spectrum
    if "E^1" in keys:

        # Build transformation matrix
        A = np.array([[1, 9 / 7, 0, 0],
                      [0, 1 / 42, 143 / 42, 11 / 42],
                      [0, 1 / 77, -130 / 77, 4 / 77],
                      [0, 1 / 462, 5 / 66, -1 / 66]])
        A[1, :] *= 225
        A[2, :] *= 1089
        A[3, :] *= 184041 / 25

        # With offset parameter
        if "E^0" in keys:
            F0, F2, F4, F6 = A @ np.array([radial[f"E^{i}"] for i in range(4)])
            for i in range(4):
                keys.remove(f"E^{i}")

        # Without offset parameter
        else:
            F0 = None
            A = A[1:, 1:]
            F2, F4, F6 = A @ np.array([radial[f"E^{i}"] for i in range(1, 4)])
            for i in range(1, 4):
                keys.remove(f"E^{i}")

        # Store the converted parameters
        if F0 is not None:
            new_radial[f"H1/0"] = F0
        new_radial[f"H1/2"] = F2
        new_radial[f"H1/4"] = F4
        new_radial[f"H1/6"] = F6

    # There should be no remaining parameters
    if len(keys) != 0:
        raise ValueError(f"Unknown radial integrals: {', '.join(keys)}!")

    # Return normalised set of radial integrals
    return new_radial


def get_energies(config, radial, state_space, base_transform):
    """ Build and diagonalize the matrix of a perturbation hamiltonian operator as linear combination of the
    interaction hamiltonians and factors specified in the dictionary radial in the given coupling scheme."""

    assert state_space in ("sljm", "slj")
    assert isinstance(radial, dict)

    if any(k.startswith("Hcf/") for k in radial.keys()):
        assert state_space == "sljm", "Crystal field parameters require SLJM coupling!"

    # Build linear combination specified in the radial dictionary
    radial = normalize_radial(radial)
    names = {k: v for k, v in radial.items() if k != "base" and not k.startswith("JO/")}

    # Build the perturbation Hamiltonian
    if state_space == "sljm":
        H = get_matrix(names, config, "product")
    else:
        H = get_matrix(names, config, "slj")

    # Diagonalise the Hamiltonian and get energies and intermediate coupling vectors
    energies, transform = np.linalg.eigh(H)

    # Adjust the energy level of the ground state
    if "base" in radial:
        energies += radial["base"] - energies[0]

    # Back transformation to SLJM states
    if state_space == "sljm":
        transform = base_transform.T @ transform

    # Return the total hamiltonian as Matrix object in the given coupling scheme
    return energies, transform


##########################################################################
# AMELI interface
##########################################################################

def get_matrix(name, config, state_space):
    """ Return the matrix of the tensor operator or the weighted sum of tensor operators in the coupling scheme. """

    assert isinstance(name, (str, dict))
    assert isinstance(config, str)
    assert isinstance(state_space, str)

    # Create cache file folder
    path = MATRIX_PATH / config
    if not path.exists():
        path.mkdir(parents=True)

    # Linear combination of matrices
    if isinstance(name, dict):
        matrix = 0.0
        for sub_name, sub_weight in name.items():
            matrix += sub_weight * get_matrix(sub_name, config, state_space)
        return matrix

    # Use alternative name
    if name in ALT_NAMES:
        return get_matrix(ALT_NAMES[name], config, state_space)

    # Read matrix
    vault = path / f"{state_space}.hdf5"
    with h5py.File(vault, "a") as fp:

        # Convert and store AMELI matrix
        if name not in fp.keys():
            matrix = get_ameli_matrix(name, config, state_space)
            fp.create_dataset(name, data=matrix)

        # Read converted matrix from HDF5 cache
        else:
            matrix = np.array(fp[name])

    # Return float representation of the matrix
    return matrix
