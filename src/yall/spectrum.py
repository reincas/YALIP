##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger("yall.spectrum")

# Physical constants
CONST_e = 1.6022e-19  # C
CONST_eps0 = 8.8542e-12  # C / V m
CONST_me = 9.1095e-31  # kg
CONST_h = 6.6262e-34  # J s
CONST_c = 2.99792458e8  # m / s
CONST_gs = 2.00231924


@dataclass
class Reduced:
    """ Dataclass containing the reduced matrix elements required for the calculation of the line strength
    of electric and magnetic dipole transitions. """

    U2: np.ndarray
    U4: np.ndarray
    U6: np.ndarray
    LS: np.ndarray


@dataclass
class LineStrength:
    """ Dataclass containing the matrices of the line strengths of electric and magnetic dipole transitions. """

    Sed: np.ndarray
    Smd: np.ndarray


def delta_k(k):
    """ Return the wavenumber difference matrix delta_k[i, j] = |k[i] - k[j]| ."""

    k = np.asarray(k)
    return np.abs(k[:, np.newaxis] - k[np.newaxis, :])


class Cauchy:
    """ Function class to determine the spectral refractive index using the Schott Glass Dispersion formula
    n = A*L^-4 + B*L^-2 + C + D*L^2 + E*L^4. """

    def __init__(self, A, B, C, D, E):
        """ Store Cauchy coefficients based on the wavelength unit µm. """

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def refractive_index(self, k):
        """ Return the refractive index for a single wavenumber or a numpy array of wavenumbers in cm^-1. """

        # Convert cm^-1 to microns
        lam = 10000.0 / np.asanyarray(k, dtype=float)

        # Calculate refractive index based on Schott Glass Dispersion formula
        n = (self.A * lam ** -4 +
             self.B * lam ** -2 +
             self.C +
             self.D * lam ** 2 +
             self.E * lam ** 4)

        # Return refractive index
        if np.isscalar(k):
            return float(n)
        return n


class Sellmeier:
    """ Function class to determine the spectral refractive index using the Sellmeier formula
    n^2 = 1 + sum_i(B_i * L^2 / (L^2 - C_i). """

    def __init__(self, B1, B2, B3, C1, C2, C3):
        self.B = np.array([B1, B2, B3])
        self.C = np.array([C1, C2, C3])

    def refractive_index(self, k):
        """ Return the refractive index for a single wavenumber or a numpy array of wavenumbers in cm^-1. """

        # Handle scalar vs array input
        is_scalar = np.isscalar(k)

        # Convert cm^-1 to microns (L)
        lam2 = (10000.0 / np.asanyarray(k, dtype=float)) ** 2

        # Calculate refractive index based on the Sellmeier formula
        n2 = 1.0
        for i in range(len(self.B)):
            if self.B[i] != 0:
                n2 += (self.B[i] * lam2) / (lam2 - self.C[i])
        n = np.sqrt(n2)

        # Return refractive index
        if np.isscalar(k):
            return float(n)
        return n

def line_strengths(judd_ofelt: dict, reduced: Reduced, J: list) -> LineStrength:
    """ Calculate and return matrices containing the line strengths of all electric and magnetic dipole
    transitions. Rows refer to final and columns to initial states. """

    # Multiply each matrix column with factor with J of the initial state
    invJi = np.array([1 / (3 * (2 * float(j) + 1)) for j in J])

    result_ed = (judd_ofelt["JO/2"] * reduced.U2 +
                 judd_ofelt["JO/4"] * reduced.U4 +
                 judd_ofelt["JO/6"] * reduced.U6) * invJi
    result_md = reduced.LS * invJi

    # Remove diagonal elements
    np.fill_diagonal(result_ed, 0.0)
    np.fill_diagonal(result_md, 0.0)

    # Apply scaling factors
    result_ed *= CONST_e ** 2 / (4 * np.pi * CONST_eps0) * 1e-24
    result_md *= 1 / (4 * np.pi * CONST_eps0) * (CONST_e * CONST_h / (2 * np.pi * 2 * CONST_me * CONST_c)) ** 2
    return LineStrength(Sed=result_ed, Smd=result_md)
