##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides access to states in intermediate coupling, their
# energy levels and radiative transitions based on given radial integrals
# and Judd-Ofelt parameters through the class `Levels`.
#
##########################################################################

import logging
import numpy as np

from .matrix import get_energies
from .spectrum import line_strengths, Dipole, CONST_gs, oscillator_strengths, radiative_rates
from .states import Coupling, State, StateList, get_states

logger = logging.getLogger("yalip.levels")


class IntermediateState(State):
    """ Class for an electron state in an intermediate coupling of SLJM or SLJ states. """

    def __init__(self, energy, values, states):
        """ Store the energy level and the linear combination vector of the state. """

        assert isinstance(energy, float)
        assert isinstance(values, np.ndarray)
        assert isinstance(states, list)
        assert len(states) == len(values)

        # Energy level of the state
        self.energy = energy

        # Sort states by decreasing weight
        weights = np.abs(values * values.conjugate())
        indices = list(reversed(np.argsort(weights)))

        # Linear combination vector with respect of the base states
        self.values = values[indices]

        # Weight factor vector
        self.weights = weights[indices]
        assert abs(sum(self.weights) - 1.0) < 1e-7

        # List of base states
        self.states = [states[i] for i in indices]

        # Initialise attributes of the parent object
        ref = self.states[0]
        super().__init__(ref.config, ref.coupling, ref.J, ref.mult)

    def short(self):
        """ Return a short string representation of the state. """

        return self.states[0].short()

    def long(self, min_weight=0.0):
        """ Return a long string representation of the state. """

        indices = [i for i in range(len(self.states)) if self.weights[i] >= min_weight]
        return " + ".join([f"{self.weights[i]:.2f} {self.states[i].short()}" for i in indices])


class IntermediateList(StateList):
    """ Class containing a list of IntermediateState objects. """

    def __init__(self, base_states, energies, transform):
        """ Store the given SLJM states, the energies of all states and the matrix containing the linear combination
        (columns) vectors of all states. """

        assert isinstance(base_states, StateList)
        assert len(transform.shape) == 2 and transform.shape[0] == transform.shape[1]
        assert len(energies) == len(base_states) == transform.shape[0]

        # List of base states
        self.base_states = base_states
        coupling = base_states.coupling
        assert coupling in (Coupling.SLJM, Coupling.SLJ)

        # List of energy levels
        self.energies = energies

        # J quantum number of each state in intermediate coupling is taken from its main SLJ component
        if coupling == Coupling.SLJM:
            self.J = None
        else:
            weight = np.abs(transform * transform.conj())
            self.J = [base_states.J[i] for i in np.argmax(weight, axis=0)]

        # Build list of IntermediateState objects
        states = []
        for i in range(len(base_states)):
            if coupling == Coupling.SLJM:
                states.append(IntermediateState(energies[i], transform[:, i], list(base_states.states)))
            else:
                indices = np.array(np.argwhere(np.array(base_states.J) == self.J[i]).flat)
                states.append(IntermediateState(energies[i], transform[indices, i], [base_states[i] for i in indices]))

        # Initialise attributes of the parent object
        super().__init__(base_states.config, coupling, states, transform)

    def matrix(self, name):
        """ Return the matrix of the tensor operator or the weighted sum of tensor operators. """

        matrix = super().matrix(name)
        return self.transform.T @ matrix @ self.transform

    def reduced(self, name):
        """ Return the reduced matrix of the given tensor operator or the weighted sum of tensor operators. """

        assert self.coupling == Coupling.SLJ
        matrix = super().reduced(name)
        return self.transform.T @ matrix @ self.transform


class Levels:
    """ Intermediate coupling class providing energy levels and radiative transitions. """

    def __init__(self, config, coupling, radial, jo=None, material=None):

        assert isinstance(config, str)
        assert isinstance(coupling, Coupling)
        assert isinstance(radial, dict)
        assert coupling in (Coupling.SLJ, Coupling.SLJM)
        assert jo is None or isinstance(jo, dict)

        # Electron configuration
        self.config = config

        # Coupling scheme
        self.coupling = coupling

        # Radial integrals
        assert "base" in radial
        self.radial_integrals = radial.copy()

        # Judd-Ofelt parameters
        self.judd_ofelt = jo

        # Material object providing spectral refractive indices
        self.material = material

        # Basis states
        if any(key.startswith("Hcf/") for key in radial):
            assert self.coupling == Coupling.SLJM
        self.base_states = get_states(self.config, self.coupling)

        # Intermediate states
        state_space = self.base_states.state_space
        base_transform = self.base_states.transform
        energies, transform = get_energies(self.config, self.radial_integrals, state_space, base_transform)
        self.states = IntermediateList(self.base_states, energies, transform)

        # Matrix elements required to calculate dipole transition strengths
        if self.coupling == Coupling.SLJM:
            self.dipole = None
        else:
            self.dipole = Dipole(
                U2=np.power(self.reduced("U/2"), 2),
                U4=np.power(self.reduced("U/4"), 2),
                U6=np.power(self.reduced("U/6"), 2),
                LS=np.power(self.reduced({"L": 1.0, "S": CONST_gs}), 2))

    @property
    def energies(self):
        """ List of state energies. """

        return self.states.energies

    @property
    def mult(self):
        """ List of state multiplicities. """

        return self.states.mult

    def __len__(self):
        """ Returns the number of states. """

        return len(self.states.states)

    def __iter__(self):
        """ Generate each state. """

        for state in self.states.states:
            yield state

    def __getitem__(self, item):
        """ Return the state object with the given index. """

        return self.states.states[item]

    def set_base(self, base):
        """ Set the base state energy in cm^-1. """
        
        self.states.energies += base - self.states.energies[0]
        for i, state in enumerate(self.states.states):
            state.energy = self.states.energies[i]
        self.radial_integrals["base"] = base

    def matrix(self, name):
        """ Return the matrix of the tensor operator or the weighted sum of tensor operators. """

        assert isinstance(name, (str, dict))
        return self.states.matrix(name)

    def reduced(self, name):
        """ Return the reduced matrix of the given tensor operator or the weighted sum of tensor operators. """

        assert isinstance(name, (str, dict))
        assert self.coupling == Coupling.SLJ
        return self.states.reduced(name)

    def line_strengths(self):
        """ Calculate and return matrices containing the line strengths in Jm^3 of all electric and magnetic dipole
        transitions. Rows refer to final and columns to initial states. """

        if self.coupling == Coupling.SLJM:
            raise NotImplementedError("Not yet implemented for SLJM coupling!")
        return line_strengths(self.judd_ofelt, self.dipole, self.mult)

    def oscillator_strengths(self):
        """ Calculate and return matrices containing the dimensionless oscillator strengths of all electric and
        magnetic dipole transitions. Rows refer to final and columns to initial states. """

        if self.coupling == Coupling.SLJM:
            raise NotImplementedError("Not yet implemented for SLJM coupling!")
        return oscillator_strengths(self.judd_ofelt, self.dipole, self.mult, self.energies, self.material)

    def radiative_rates(self):
        """ Calculate and return matrices containing the radiative emission rates in 1/s of all electric and magnetic
        dipole transitions. Rows refer to final and columns to initial states. """

        if self.coupling == Coupling.SLJM:
            raise NotImplementedError("Not yet implemented for SLJM coupling!")
        return radiative_rates(self.judd_ofelt, self.dipole, self.mult, self.energies, self.material)

    def life_times(self):
        """ Calculate and return the radiative lifetime of each state in s. """

        if self.coupling == Coupling.SLJM:
            raise NotImplementedError("Not yet implemented for SLJM coupling!")
        A = self.radiative_rates()
        A = A.ed + A.md
        with np.errstate(divide='ignore'):
            tau = 1.0 / np.sum(A, axis=0)
        return tau

    def branching_ratios(self):
        """ Calculate and return matrices containing the radiative branching ratios of all dipole transitions.
        Rows refer to final and columns to initial states. """

        if self.coupling == Coupling.SLJM:
            raise NotImplementedError("Not yet implemented for SLJM coupling!")
        A = self.radiative_rates()
        A = A.ed + A.md
        with np.errstate(divide='ignore', invalid='ignore'):
            beta = A / np.sum(A, axis=0, keepdims=True)
        beta = np.triu(beta, k=1)
        return beta

    def str_levels(self, min_weight=0.0):
        """ Return a list containing an extensive description string for each state with energy and composition with
        weights in intermediate coupling. Only SLJ states with the given minimum weight are included. """

        for state in self.states:
            yield f"  {state.energy:7.0f} | {state.long(min_weight)} >"
