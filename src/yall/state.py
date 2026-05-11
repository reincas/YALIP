##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides classes for the representation of single states
# and lists of all states of a configuration in the determinantal product
# state space and in SLJM, SLJ and intermediate coupling. The main tasks
# of these objects is the generation of short and long string
# representations of the states and the provision of transformation
# matrices between the different coupling schemes.
#
##########################################################################

from enum import Enum
import numpy as np

from .ameli import read_transform

SPECTRAL = "spdfghiklmnoqrtuvwxyz"


class Coupling(Enum):
    """ This enumeration class is used to mark the four coupling schemes used in the Lanthanide package: determinantal
    product state coupling, SLJM coupling, SLJ coupling, and intermediate coupling. """

    Product = 0
    SLJM = 1
    SLJ = 2
    J = 3
    M = 4


def norm_magnetic(value):
    if value.endswith("/2"):
        v = int(value[:-2])
    else:
        v = int(value)
    if v > 0 and not value.startswith("+"):
        value = "+" + value
    else:
        value = str(value)
    return value


class StateList:
    """ Abstract class for a list of electron states in a certain coupling scheme. The abstract defines some common
    methods acting on the common attribute states, which contains the list of State objects. """

    states = []

    def __len__(self):
        """ Returns the number of states. """

        return len(self.states)

    def __iter__(self):
        """ Generate each state. """

        for state in self.states:
            yield state

    def __getitem__(self, item: int):
        """ Return the state object with the given index. """

        return self.states[item]

    def short(self):
        """ Return a list containing the short names of all states. """

        return [state.short() for state in self.states]

    def long(self):
        """ Return a list containing the long names of all states. """

        return [state.long() for state in self.states]


##########################################################################
# Product states
##########################################################################

class StateProduct:
    """ Class for a determinantal product state. """

    def __init__(self, quantum):
        """ Store the quantum numbers of the electrons of the state. """

        self.quantum = quantum

    def short(self):
        """ Return a short string representation of the state. """

        quantum = []
        for single in self.quantum:
            ml = norm_magnetic(single["ml"])
            ms = "d" if single["ms"] == "-1/2" else "u"
            quantum.append(ml + ms)
        return " ".join(quantum)

    def long(self):
        """ Return a long string representation of the state. """

        quantum = []
        for single in self.quantum:
            l = SPECTRAL[single["l"]]
            ml = norm_magnetic(single["ml"])
            s = single["s"]
            ms = norm_magnetic(single["ms"])
            quantum.append(f"{l},{ml},{s},{ms}")
        return " ".join(quantum)

    def __getitem__(self, item: int):
        """ Return the quantum number dictionary of the state with the given index. """

        return self.quantum[item]

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListProduct(StateList):
    """ Class containing a list of StateProduct objects representing an electron configuration. """

    def __init__(self, values, pool):
        """ Store the given array of electron indices and build the list of StateProduct objects. """

        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2

        # Pool of single electron quantum numbers
        self.pool = pool

        # Array of indices into the electron pool
        self.values = values

        # List of state objects
        self.states = [StateProduct([self.pool[i] for i in state]) for state in values]

        # No transformation matrix
        self.transform = None


##########################################################################
# SLJM states
##########################################################################

class StateSLJM:
    """ Class for an electron state in SLJM coupling following the chain of symmetry operators. """

    def __init__(self, quantum):
        self.quantum = quantum
        self.J = quantum["J2"]

    def short(self):
        """ Return a short string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']} {self['J2']} {self['Jz']}"

    def long(self):
        """ Return a long string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']} {self['sen']} {self['C7']}{self['C2']}{self['tau']} {self['J2']} {self['Jz']}"

    def __getitem__(self, sym: str):
        """ Return the irreducible representation of the given symmetry operator. """

        if sym not in self.quantum:
            raise KeyError(f"Unknown symmetry {sym}!")
        return self.quantum[sym]

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListSLJM(StateList):
    """ Class containing a list of StateSLJM objects representing an electron configuration. """

    def __init__(self, values, chain, repr, transform):
        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2
        assert values.shape[1] == len(chain)
        assert set(chain) == set(repr.keys())

        self.values = values
        self.chain = chain
        self.repr = repr
        self.transform = transform

        self.states = [StateSLJM({sym: repr[sym][state[i]] for i, sym in enumerate(chain)}) for state in values]
        self.J = [state.J for state in self.states]

    def stretched(self):
        return [i for i, state in enumerate(self.states) if state["Jz"].replace("+", "") == state["J2"]]


##########################################################################
# SLJ states
##########################################################################

class StateSLJ:
    """ Class for an electron state in SLJ coupling following the chain of symmetry operators. """

    def __init__(self, quantum):
        self.quantum = quantum
        self.J = quantum["J2"]

    def short(self):
        """ Return a short string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']} {self['J2']}"

    def long(self):
        """ Return a long string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']} {self['sen']} {self['C7']}{self['C2']}{self['tau']} {self['J2']}"

    def __getitem__(self, sym: str):
        """ Return the irreducible representation of the given symmetry operator. """

        if sym not in self.quantum:
            raise KeyError(f"Unknown symmetry {sym}!")
        return self.quantum[sym]

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListSLJ(StateList):
    """ Class containing a list of StateSLJ objects representing an electron configuration. """

    def __init__(self, values, chain, repr, transform):
        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2
        assert values.shape[1] == len(chain)
        assert set(chain) == set(repr.keys())

        self.values = values
        self.chain = chain
        self.repr = repr
        self.transform = transform

        self.states = [StateSLJM({sym: repr[sym][state[i]] for i, sym in enumerate(chain)}) for state in values]
        self.J = [state.J for state in self.states]


##########################################################################
# Intermediate SLJ states
##########################################################################

class StateJ:
    """ Class for an electron state in an intermediate coupling of SLJ states. """

    def __init__(self, energy, values, states):
        """ Store the energy level of the state and the vector (values) for the linear combination of the given
        SLJ states. """

        assert isinstance(energy, float)
        assert isinstance(values, np.ndarray)
        assert isinstance(states, list)
        assert len(states) == len(values)
        assert len(set(state["J2"].J for state in states)) == 1

        # Energy level of the state
        self.energy = energy

        # Sort states by weight
        weights = np.abs(values * values.conjugate())
        indices = list(reversed(np.argsort(weights)))

        # Linear combination vector and weight factor vector
        self.values = values[indices]
        self.weights = weights[indices]
        assert abs(sum(self.weights) - 1.0) < 1e-7

        # The state in intermediate coupling is a linear combination of this list of related SLJ states
        self.states = [states[i] for i in indices]

        # Common quantum number J of the total angular momentum of all related SLJ states
        self.J = self.states[0]["J2"].J

    def short(self):
        """ Return a short string representation of the state. """

        return self.states[0].short()

    def long(self, min_weight=0.0):
        """ Return a long string representation of the state. """

        indices = [i for i in range(len(self.states)) if self.weights[i] >= min_weight]
        return " + ".join([f"{self.weights[i]:.2f} {self.states[i].short()}" for i in indices])

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListJ(StateList):
    """ Class containing a list of StateJ objects representing an electron state in an intermediate SLJ coupling. """

    def __init__(self, slj_states, energies, transform):
        """ Store the given StateListSLJ object with the list of related SLJ states, the energies of all states and
        the transformation matrix representing the linear combination of the SLJ states. """

        assert isinstance(slj_states, StateListSLJ)
        assert len(transform.shape) == 2 and transform.shape[0] == transform.shape[1]
        assert len(energies) == len(slj_states) == transform.shape[0]

        # Store the SLJ states, energy levels and the transformation matrix from SLJ to intermediate coupling
        self.slj_states = slj_states
        self.energies = energies
        self.transform = transform

        # Matrix of weight factors for the SLJ components or each state in intermediate coupling
        weight = np.abs(self.transform * self.transform.conj())

        # J quantum number of each state in intermediate coupling is taken from its main SLJ component
        self.J = [self.slj_states.J[i] for i in np.argmax(weight, axis=0)]

        # Build the list of StateJ objects
        self.states = []
        for i in range(len(self.slj_states)):
            # Indices of all SLJ states with the same quantum number J as the current state
            indices = np.array(np.argwhere(self.slj_states.J == self.J[i]).flat)

            # Combination vector and SLJ states for the linear combination of the current state
            values = self.transform[indices, i]
            slj_states = [self.slj_states[i] for i in indices]

            # Add StateJ object representing the current state in intermediate coupling
            self.states.append(StateJ(energies[i], values, slj_states))


##########################################################################
# Intermediate SLJM states
##########################################################################

class StateM:
    """ Class for an electron state in an intermediate coupling of SLJM states. """

    def __init__(self, energy, values, states):
        """ Store the energy level of the state and the vector (values) for the linear combination of the given
        SLJM states. """

        assert isinstance(energy, float)
        assert isinstance(values, np.ndarray)
        assert isinstance(states, list)
        assert len(states) == len(values)

        # Energy level of the state
        self.energy = energy

        # Sort states by weight
        weights = np.abs(values * values.conjugate())
        indices = list(reversed(np.argsort(weights)))

        # Linear combination vector and weight factor vector
        self.values = values[indices]
        self.weights = weights[indices]
        assert abs(sum(self.weights) - 1.0) < 1e-7

        # The state in intermediate coupling is a linear combination of this list of related SLJM states
        self.states = [states[i] for i in indices]

    def short(self):
        """ Return a short string representation of the state. """

        return self.states[0].short()

    def long(self, min_weight=0.0):
        """ Return a long string representation of the state. """

        indices = [i for i in range(len(self.states)) if self.weights[i] >= min_weight]
        return " + ".join([f"{self.weights[i]:.2f} {self.states[i].short()}" for i in indices])

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListM(StateList):
    """ Class containing a list of StateJM objects representing an electron state in an intermediate coupling
    of SLJM states. """

    def __init__(self, sljm_states, energies, transform):
        """ Store the given StateListSLJM object with the list of related SLJM states, the energies of all states and
        the transformation matrix representing the linear combination of the SLJM states. """

        assert isinstance(sljm_states, StateListSLJM)
        assert len(transform.shape) == 2 and transform.shape[0] == transform.shape[1]
        assert len(energies) == len(sljm_states) == transform.shape[0]

        # Store the SLJM states, energy levels and the transformation matrix from SLJM to intermediate coupling
        self.sljm_states = sljm_states
        self.energies = energies
        self.transform = transform

        # Build the list of StateJM objects
        self.states = []
        for i in range(len(self.sljm_states)):
            # Add StateM object representing the current state in intermediate coupling
            values = self.transform[:, i]
            sljm_states = list(self.sljm_states.states)
            self.states.append(StateM(energies[i], values, sljm_states))


##########################################################################
# HDF5 cache interface
##########################################################################

def init_states(config_name):
    data = read_transform(config_name)

    Product_States = StateListProduct(data["rowStates"], data["electronPool"])

    values = data["colStates"]
    chain = data["tensorChain"]
    transform = data["transform"]
    reprs = data["irreducibleRepresentations"]
    SLJM_states = StateListSLJM(values, chain, reprs, transform)

    indices = SLJM_states.stretched()
    assert chain[-1] == "Jz"
    values = values[indices,:-1]
    chain = chain[:-1]
    del reprs["Jz"]
    transform = transform[:, indices]
    SLJ_states = StateListSLJ(values, chain, reprs, transform)

    # Return StateList dictionary
    return {
        Coupling.Product.name: Product_States,
        Coupling.SLJM.name: SLJM_states,
        Coupling.SLJ.name: SLJ_states,
    }
