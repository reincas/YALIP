##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np

from .ameli import update
from .matrix import get_matrix, get_energies, get_reduced
from .spectrum import Reduced, CONST_gs, line_strengths, Cauchy, Sellmeier, oscillator_strengths, radiative_rates
from .state import Coupling, init_states, StateListJ, StateListM

logger = logging.getLogger("yall.lanthanide")

# Symbols of the 15 Lanthanides. The configurations of the triply ionized atoms are 4f0 - 4f14
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]

##########################################################################
#
# Radial integrals in cm^(-1) for LaF3 from:
#
#   [5]  W. T. Carnall, H. Crosswhite, H. M. Crosswhite
#        "Energy level structure and transition probabilities of the
#        trivalent lanthanides in LaF3"
#        ANL-78-XX-95, Argonne National Laboratory Report, 1978
#        (available as microfiche)
#
#   [39] W. T. Carnall, G. L. Goodman, K. Rajnak, and R. S. Rana
#        "A systematic analysis of the spectra of the lanthanides doped
#        into single crystal LaF3"
#        J. Chem. Phys. 90 (1989), no. 7, pp. 3443-3457
#
#   [x]  Dimitar N. Petrov, B.M. Angelov
#        "Spin – Orbit interaction in Yb3+ – Ground level and nephelauxetic
#        effect in crystals"
#        Chem. Phys. 525 (2019), 110416
#        https://doi.org/10.1016/j.chemphys.2019.110416
#
#   [RC] Reinhard Caspary
#        "Applied Rare-Earth Spectroscopy for Fiber Laser Optimization"
#        Dissertation, Shaker, Aachen, 2002
#
##########################################################################

RADIAL = {
    "Ce3+":  # from [39]
        {"base": 217, "H2": 647.3},
    "Pr3+":  # from [39]
        {"base": 205, "H1/2": 68878, "H1/4": 50347, "H1/6": 32901, "H2": 751.7,
         "H3/0": 16.23, "H3/1": -566.6, "H3/2": 1371, "H5fix": 2.08, "H6fix": 88.6},
    "Pr3+/alt":  # from [5]
        {"base": 191, "H1/2": 69305, "H1/4": 50675, "H1/6": 32813, "H2": 750.8,
         "H3/0": 21, "H3/1": -842, "H3/2": 1625, "H5fix": 1.99, "H6fix": 200},
    "Pr3+/ZBLAN":  # from [RC]
        {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
         "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67},
    "Nd3+":  # from [5]
        {"base": 235, "H1/2": 73036, "H1/4": 52624, "H1/6": 35793, "H2": 884.9,
         "H3/0": 21.28, "H3/1": -583, "H3/2": 1443,
         "H4/2": 306, "H4/3": 41, "H4/4": 59, "H4/6": -283, "H4/7": 326, "H4/8": 298, "H5fix": 2.237, "H6fix": 213},
    "Pm3+":  # from [5]
        {"base": 120, "H1/2": 77000, "H1/4": 55000, "H1/6": 37500, "H2": 1022,
         "H3/0": 21.0, "H3/1": -560, "H3/2": 1400,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 2.49, "H6fix": 440},
    "Sm3+":  # from [5]
        {"base": 101, "H1/2": 79915, "H1/4": 57256, "H1/6": 40424, "H2": 1177.2,
         "H3/0": 20.07, "H3/1": -563, "H3/2": 1436,
         "H4/2": 288, "H4/3": 36, "H4/4": 56, "H4/6": -283, "H4/7": 333, "H4/8": 342, "H5fix": 2.76, "H6fix": 344},
    "Eu3+":  # from [5]
        {"base": 0, "H1/2": 84000, "H1/4": 60000, "H1/6": 42500, "H2": 1327,
         "H3/0": 20, "H3/1": -570, "H3/2": 1450,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.03, "H6fix": 300},
    "Gd3+":  # from [5]
        {"base": 2, "H1/2": 85587, "H1/4": 61361, "H1/6": 45055, "H2": 1503.5,
         "H3/0": 20, "H3/1": -590, "H3/2": 1450,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.32, "H6fix": 611},
    "Tb3+":  # from [5]
        {"base": 124, "H1/2": 91220, "H1/4": 65798, "H1/6": 43661, "H2": 1702.2,
         "H3/0": 19.81, "H3/1": -600, "H3/2": 1400,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.61, "H6fix": 583},
    "Dy3+":  # from [5]
        {"base": 175, "H1/2": 94877, "H1/4": 67470, "H1/6": 45745, "H2": 1912,
         "H3/0": 17.64, "H3/1": -608, "H3/2": 1498,
         "H4/2": 423, "H4/3": 50, "H4/4": 117, "H4/6": -334, "H4/7": 432, "H4/8": 353, "H5fix": 3.92, "H6fix": 771},
    "Ho3+":  # from [5]
        {"base": 9, "H1/2": 97025, "H1/4": 68885, "H1/6": 47744, "H2": 2144.2,
         "H3/0": 18.98, "H3/1": -579, "H3/2": 1570,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 4.25, "H6fix": 843},
    "Er3+":  # from [39]
        {"base": 219, "H1/2": 97483, "H1/4": 67904, "H1/6": 54010, "H2": 2376,
         "H3/0": 17.79, "H3/1": -582.1, "H3/2": 1800,
         "H4/2": 400, "H4/3": 43, "H4/4": 73, "H4/6": -271, "H4/7": 308, "H4/8": 299, "H5fix": 3.86, "H6fix": 594},
    "Er3+/alt":  # from [5]
        {"base": 217, "H1/2": 100274, "H1/4": 70555, "H1/6": 49900, "H2": 2381,
         "H3/0": 17.88, "H3/1": -599, "H3/2": 1719,
         "H4/2": 441, "H4/3": 42, "H4/4": 64, "H4/6": -314, "H4/7": 387, "H4/8": 363, "H5fix": 4.58, "H6fix": 852},
    "Er3+/ZBLAN":  # from [RC]
        {"base": 50.99, "H1/2": 97088.92, "H1/4": 68587.69, "H1/6": 55006.43, "H2": 2369.69,
         "H3/0": 17.81, "H3/1": -559.04, "H3/2": 1603.10,
         "H4/2": 471.61, "H4/3": 20.39, "H4/4": 18.81, "H4/6": -398.11, "H4/7": 199.03, "H4/8": 449.82,
         "H5fix": 4.66, "H6fix": 475.64},
    "Tm3+":  # from [39]
        {"base": 250, "H1/2": 100134, "H1/4": 69613, "H1/6": 55975, "H2": 2636,
         "H3/0": 17.26, "H3/1": -624.5, "H3/2": 1820, "H5fix": 3.81, "H6fix": 695},
    "Tm3+/alt":  # from [5]
        {"base": 175, "H1/2": 102459, "H1/4": 72424, "H1/6": 51380, "H2": 2640,
         "H3/0": 17, "H3/1": -737, "H3/2": 1700, "H5fix": 4.93, "H6fix": 729.6},
    "Tm3+/ZBLAN":  # from [RC]
        {"base": 149.80, "H1/2": 102403.01, "H1/4": 73241.80, "H1/6": 50320.22, "H2": 2583.47,
         "H3/0": 17.43, "H3/1": -841.95, "H3/2": 1820, "H5fix": -2.35, "H6fix": 12.78},
    "Yb3+":  # guess based on [x]
        {"base": 0, "H2": 2900},
}

# Judd-Ofelt parameters from [RC]
JUDD_OFELT = {
    "Pr3+/ZBLAN":  # from [RC]
        {"JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972},
    "Er3+/ZBLAN":  # from [RC]
        {"JO/2": 2.915, "JO/4": 1.464, "JO/6": 1.184},
    "Tm3+/ZBLAN":  # from [RC]
        {"JO/2": 2.920, "JO/4": 1.856, "JO/6": 0.670},
}

# Spectral refractive index coefficients
MATERIAL = {
    "Pb:ZBLAN": Cauchy(1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6),
    "ZBLAN": Cauchy(1.35123e-5, 2.94780e-3, 1.48965, -1.30933e-3, -3.23335e-6),
    "SiO2": Sellmeier(0.6961663, 0.4079426, 0.8974794, 0.0684043, 0.1162414, 9.896161)
}

class Lanthanide:
    """ Lanthanide ion with given number of 4f electrons. It provides all energy levels and states in intermediate
    coupling as well as the line strengths of all radiative electric and magnetic dipole transitions. """

    def __init__(self, num: int, coupling=None, radial=None):
        """ Initialize attributes, datastructures, and the file cache for the yall ion with given number of
        4f electrons. Calculate energy levels and states in intermediate coupling based on the given dictionary of
        radial integrals (default: Carnalls integrals for LaF3). The given coupling scheme may be SLJM or
        SLJ (default). """

        # Coupling scheme for the states. Default is SLJ.
        self.coupling = coupling or Coupling.SLJ

        assert isinstance(num, int)
        assert 0 < num < 14
        assert self.coupling in (Coupling.SLJ, Coupling.SLJM), f"Coupling {self.coupling} is not supported."

        # Number of electrons, name and electron configuration of the lanthanide ion
        self.num = num
        self.name = f"{LANTHANIDES[num]}3+"
        self.config = f"f{num}"

        # Update matrix data from the Zenodo repository
        update(self.config)

        # Initialize the electron states for SLJM and SLJ coupling
        self.state_dict = init_states(self.config)

        # Calculate energy levels and states in intermediate coupling based on the given radial integrals
        self.radial = None
        self.energies = None
        self.set_radial(radial or RADIAL[self.name])

        self._reduced_ = None

    def __enter__(self):
        """ Entry hook of the content management protocol. """

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ Exit hook of the content management protocol. """

        self.close()
        return False

    def close(self):

        pass

    def matrix(self, name, coupling=None):
        """ Return Matrix object of the given operator in the given coupling scheme (default: Product). """

        coupling = coupling or self.coupling
        return get_matrix(name, self.config, coupling)

    def reduced(self, name, coupling):
        """ Return array of reduced matrix elements of the given operator in the coupling scheme SLJ or
        intermediate. """

        assert coupling in (Coupling.SLJ, Coupling.Intermediate)
        matrix = get_reduced(name, self.config)
        if coupling == Coupling.Intermediate:
            V = self.states(coupling).transform
            matrix = V.T @ matrix @ V
        return matrix

    def states(self, coupling=None):
        """ Return StateList object of the given coupling scheme or the one selected when the Lanthanide object
        was initialized. """

        coupling = coupling or self.coupling
        return self.state_dict[coupling.name]

    def num_states(self, coupling=None):
        """ Return number of states in the given coupling scheme or the one selected when the Lanthanide object
        was initialized. """

        return len(self.states(coupling))

    def set_radial(self, radial):
        """ Normalise and store new dictionary of radial integrals and use them to calculate all energy levels and
        states in intermediate coupling. """

        assert isinstance(radial, dict)

        # Do nothing, if the given integrals are identical to the current ones
        if radial == self.radial:
            return

        # Diagonalise the Hamiltonian and get energies and intermediate coupling vectors
        self.radial = radial
        self.energies, transform = get_energies(self.config, self.radial, self.states(self.coupling))

        # Adjust the energy level of the ground state
        if "base" in self.radial:
            self.energies += self.radial["base"] - self.energies[0]

        # Build and store the StateListJ object of the electron states in intermediate coupling
        if self.coupling == Coupling.SLJ:
            self.state_dict["Intermediate"] = StateListJ(self.state_dict["SLJ"], self.energies, transform)
        else:
            self.state_dict["Intermediate"] = StateListM(self.state_dict["SLJM"], self.energies, transform)

        # Invalidate previously calculated reduced matrix elements
        self._reduced_ = None

    def line_reduced(self):
        """ Calculate and return all reduced matrix elements required for the calculation of line strengths of
        radiative electric or magnetic dipole transitions if they were not yet calculated. Rows refer to final and
        columns to initial states. """

        if not self._reduced_:
            self._reduced_ = Reduced(
                U2=np.power(self.reduced("U/2", Coupling.Intermediate), 2),
                U4=np.power(self.reduced("U/4", Coupling.Intermediate), 2),
                U6=np.power(self.reduced("U/6", Coupling.Intermediate), 2),
                LS=np.power(self.reduced([(1.0, "L"), (CONST_gs, "S")], Coupling.Intermediate), 2))
        return self._reduced_

    def line_strengths(self, judd_ofelt):
        """ Calculate and return matrices containing the line strengths in Jm^3 of all electric and magnetic dipole
        transitions. Rows refer to final and columns to initial states. """

        # Multiply each matrix column with factor with J of the initial state
        reduced = self.line_reduced()
        states = self.states(Coupling.Intermediate)
        J = list(map(float, states.J))
        return line_strengths(judd_ofelt, reduced, J)

    def oscillator_strengths(self, judd_ofelt, material):
        """ Calculate and return matrices containing the dimensionless oscillator strengths of all electric and
        magnetic dipole transitions. Rows refer to final and columns to initial states. """

        # Multiply each matrix column with factor with J of the initial state
        reduced = self.line_reduced()
        states = self.states(Coupling.Intermediate)
        J = list(map(float, states.J))
        return oscillator_strengths(judd_ofelt, reduced, J, self.energies, material)

    def radiative_rates(self, judd_ofelt, material):
        """ Calculate and return matrices containing the radiative emission rates in 1/s of all electric and magnetic
        dipole transitions. Rows refer to final and columns to initial states. """

        # Multiply each matrix column with factor with J of the initial state
        reduced = self.line_reduced()
        states = self.states(Coupling.Intermediate)
        J = list(map(float, states.J))
        return radiative_rates(judd_ofelt, reduced, J, self.energies, material)

    def str_levels(self, min_weight=0.0):
        """ Return a list containing an extensive description string for each state with energy and composition with
        weights in intermediate coupling. Only SLJ states with the given minimum weight are included. """

        for state in self.states(Coupling.Intermediate):
            yield f"  {state.energy:7.0f} | {state.long(min_weight)} >"

    def __str__(self):
        """ Return a short string representation of this Lanthanide object. """

        return f"{self.name} ({self.config})"
