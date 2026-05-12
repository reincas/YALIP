##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides the calculation of matrices of higher-order tensor
# operators in the determinantal product space of a given electron
# configuration based on unit tensor operators from the module "unit".
#
# The function get_matrix() returns a Matrix object containing a given
# tensor operator matrix in a given coupling scheme. The matrices of
# certain tensor operators, namely the perturbation hamiltonians are
# eventually stored in the HDF5 file cache in SLJM and SLJ coupling.
#
# Matrix objects support diagonalisation and elementary arithmetic
# operations to build for linear combinations of these matrices.
#
##########################################################################
from pathlib import Path

import h5py
import numpy as np

from .ameli import MATRIX_PATH, get_ameli_matrix
from .state import Coupling


def normalise_radial(radial):
    """ Convert an alternative parameter set into the standard set of radial integrals "H1"-"H6". """

    # Get keys and initialize the normalised set.
    keys = list(radial.keys())
    new_radial = {}

    # No conversion for parameter "base"
    if "base" in keys:
        new_radial["base"] = radial["base"]
        keys.remove("base")

    # No conversion for standard parameters
    for key in list(keys):
        if key[:2] in ("H1", "H2", "H3", "H4", "H5", "H6") or key[:3] == "Hcf":
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


def build_hamilton(ion, radial: dict, coupling=None):
    """ Build and return the matrix of a perturbation hamiltonian operator as linear combination of the interaction
    hamiltonians and factors specified in the dictionary radial in the given coupling scheme."""

    assert coupling is None or coupling in (Coupling.SLJM, Coupling.SLJ)
    assert isinstance(radial, dict)

    # Default is coupling of yall ion
    coupling = coupling or ion.coupling

    # Initialize empty hamiltonian matrix
    num_states = len(ion.states(coupling))
    matrix = np.zeros((num_states, num_states), dtype=float)

    # Build linear combination specified in the radial dictionary
    for name in radial:
        if name == "base":
            continue
        if name[:1] != "H":
            continue
        if name[:3] == "Hcf" and coupling != Coupling.SLJM:
            raise ValueError("Crystal field interaction requires SLJM coupling!")
        matrix = matrix + radial[name] * get_matrix(ion, name, coupling).array

    # Return the total hamiltonian as Matrix object in the given coupling scheme
    return matrix


def get_matrix(name, config, coupling):
    """ Return a Matrix object of the matrix of the tensor operator with given name and in the given coupling
    scheme (default: Coupling.Product). The matrices of the tensor operators in the list STORE are stored in the
    HDF5 file cache in SLJM and SLJ coupling. """

    assert isinstance(coupling, Coupling)
    assert coupling in (Coupling.Product, Coupling.SLJM, Coupling.SLJ)

    path = MATRIX_PATH / config
    if not path.exists():
        path.mkdir(parents=True)

    vault = path / f"{coupling.name.lower()}.hdf5"
    with h5py.File(vault, "a") as fp:
        if name not in fp.keys():
            matrix = get_ameli_matrix(name, config, coupling)
            fp.create_dataset(name, data=matrix)
        else:
            matrix = np.array(fp[name])
    return matrix
