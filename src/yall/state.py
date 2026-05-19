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
import logging
import numpy as np

from .ameli import read_transform

logger = logging.getLogger("yall.state")

SPECTRAL = "spdfghiklmnoqrtuvwxyz"


class Coupling(Enum):
    """ This enumeration class is used to mark the four coupling schemes used in the Lanthanide package: determinantal
    product state coupling, SLJM coupling, SLJ coupling, and intermediate coupling. """

    Product = 0
    SLJM = 1
    SLJ = 2
    Intermediate = 999 # Todo ############################################################


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


##########################################################################
# Abstract states
##########################################################################

class State:
    """ Abstract class for an electron state in a certain coupling scheme. """

    coupling: Coupling
    is_intermediate: bool

    def __init__(self, config, mult):
        """ Store electron configuration and state multiplicity. """

        assert isinstance(config, str)
        assert isinstance(mult, int)

        # Electron configuration
        self.config = config

        # State multiplicity
        self.mult = mult


class StateList:
    """ Abstract class for a list of electron states in a certain coupling scheme. """

    coupling: Coupling
    is_intermediate: bool

    def __init__(self, config, states, transform):
        """ Store electron configuration, list of state objects, transformations matrix and list of state
        multiplicities. """

        # Electron configuration
        self.config = config

        # List of state objects
        self.states = states

        # Transformation matrix
        self.transform = transform

        # List of state multiplicities
        self.mult = [state.mult for state in self.states]

    def __len__(self):
        """ Returns the number of states. """

        return len(self.states)

    def __iter__(self):
        """ Generate each state. """

        for state in self.states:
            yield state

    def __getitem__(self, item):
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

class StateProduct(State):
    """ Class for a determinantal product state. """

    coupling = Coupling.Product
    is_intermediate = False

    def __init__(self, config, quantum):
        """ Store electron configuration and single electron quantum numbers. """

        # Single electron quantum numbers
        self.quantum = quantum

        # Initialise attributes of the parent object
        super().__init__(config, 1)

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

    coupling = Coupling.Product
    is_intermediate = False

    def __init__(self, config, values, pool):
        """ Store electron configuration, single electron pool indices, and pool of single electron quantum numbers.
        Build the list of product states. """

        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2

        # Array of indices into the electron pool
        self.values = values

        # Pool of single electron quantum numbers
        self.pool = pool

        # List of state objects
        states = [StateProduct(config, [self.pool[i] for i in state]) for state in values]

        # Initialise attributes of the parent object
        super().__init__(config, states, None)


##########################################################################
# SLJM states
##########################################################################

class StateSLJM(State):
    """ Class for an electron state in SLJM coupling following the chain of symmetry operators. """

    coupling = Coupling.SLJM
    is_intermediate = False

    def __init__(self, config, quantum):
        """ Store electron configuration and dictionary of quantum numbers. """

        # Quantum number dictionary
        self.quantum = quantum

        # Quantum number of total angular momentum as string
        self.J = quantum["J2"]

        # Initialise attributes of the parent object
        super().__init__(config, 1)

    def short(self):
        """ Return a short string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']}_{self['J2']}_{self['Jz']}"

    def long(self):
        """ Return a long string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']}_{self['sen']}_{self['C7']}{self['C2']}{self['tau']}_{self['J2']}_{self['Jz']}"

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

    coupling = Coupling.SLJM
    is_intermediate = False

    def __init__(self, config, values, chain, repr, transform):
        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2
        assert values.shape[1] == len(chain)
        assert set(chain) == set(repr.keys())

        # Electron configuration
        self.config = config

        # Dictionary of indices into the lists of irreducible representations for all states
        self.values = values

        # Chain of symmetry operator names
        self.chain = chain

        # Dictionary of lists of all irreducible representations
        self.repr = repr

        # List of SLJM states
        states = [StateSLJM(config, {sym: repr[sym][s[i]] for i, sym in enumerate(chain)}) for s in values]

        # Values of total angular momentum as list of strings
        self.J = [state.J for state in states]

        # Initialise attributes of the parent object
        super().__init__(config, states, transform)

    def stretched(self):
        """ Return a list indices containing all stretched states. """

        return [i for i, state in enumerate(self.states) if state["Jz"].replace("+", "") == state["J2"]]


##########################################################################
# SLJ states
##########################################################################

class StateSLJ(State):
    """ Class for an electron state in SLJ coupling following the chain of symmetry operators. """

    coupling = Coupling.SLJ
    is_intermediate = False

    def __init__(self, config, quantum):
        """ Store electron configuration and dictionary of quantum numbers. """

        # Quantum number dictionary
        self.quantum = quantum

        # Quantum number of total angular momentum as string
        self.J = quantum["J2"]

        # State multiplicity
        if self.J.endswith("/2"):
            mult = int(self.J[:-2]) + 1
        else:
            mult = 2 * int(self.J) + 1

        # Initialise attributes of the parent object
        super().__init__(config, mult)

    def short(self):
        """ Return a short string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']}_{self['J2']}"

    def long(self):
        """ Return a long string representation of the state. """

        return f"{self['S2']}{self['L2']}{self['num']}_{self['sen']}_{self['C7']}{self['C2']}{self['tau']}_{self['J2']}"

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

    coupling = Coupling.SLJ
    is_intermediate = False

    def __init__(self, config, values, chain, repr, transform):
        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2
        assert values.shape[1] == len(chain)
        assert set(chain) == set(repr.keys())

        # Electron configuration
        self.config = config

        # Dictionary of indices into the lists of irreducible representations for all states
        self.values = values

        # Chain of symmetry operator names
        self.chain = chain

        # Dictionary of lists of all irreducible representations
        self.repr = repr

        # List of SLJ states
        states = [StateSLJ(config, {sym: repr[sym][s[i]] for i, sym in enumerate(chain)}) for s in values]

        # Values of total angular momentum as list of strings
        self.J = [state.J for state in states]

        # Initialise attributes of the parent object
        super().__init__(config, states, transform)


##########################################################################
# Intermediate SLJ states
##########################################################################

class StateJ(State):
    """ Class for an electron state in an intermediate coupling of SLJ states. """

    coupling = Coupling.SLJ
    is_intermediate = True

    def __init__(self, energy, values, states):
        """ Store the energy level of the state and the linear combination vector of the given SLJ states. """

        assert isinstance(energy, float)
        assert isinstance(values, np.ndarray)
        assert isinstance(states, list)
        assert len(states) == len(values)
        assert len(set(state.J for state in states)) == 1

        # Energy level of the state
        self.energy = energy

        # Sort states by decreasing weight
        weights = np.abs(values * values.conjugate())
        indices = list(reversed(np.argsort(weights)))

        # Linear combination vector
        self.values = values[indices]

        # Weight factor vector
        self.weights = weights[indices]
        assert abs(sum(self.weights) - 1.0) < 1e-7

        # The state in intermediate coupling is a linear combination of this list of related SLJ states
        self.states = [states[i] for i in indices]

        # Quantum number of total angular momentum as string
        self.J = self.states[0].J

        # Initialise attributes of the parent object
        super().__init__(self.states[0].config, self.states[0].mult)

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

    coupling = Coupling.SLJ
    is_intermediate = True

    def __init__(self, slj_states, energies, transform):
        """ Store the given SLJ states, the energies of all states and the matrix containing the linear combination
        (columns) vectors of all states. """

        assert isinstance(slj_states, StateListSLJ)
        assert len(transform.shape) == 2 and transform.shape[0] == transform.shape[1]
        assert len(energies) == len(slj_states) == transform.shape[0]

        # List of SLJ states
        self.slj_states = slj_states

        # List of energy levels
        self.energies = energies

        # Matrix of weight factors for the SLJ components or each state in intermediate coupling
        weight = np.abs(transform * transform.conj())

        # J quantum number of each state in intermediate coupling is taken from its main SLJ component
        self.J = [slj_states.J[i] for i in np.argmax(weight, axis=0)]

        # Build list of StateJ objects
        states = []
        for i in range(len(slj_states)):
            # Indices of all SLJ states with the same quantum number J as the current state
            indices = np.array(np.argwhere(np.array(slj_states.J) == self.J[i]).flat)

            # Add StateJ object representing the current state in intermediate coupling
            states.append(StateJ(energies[i], transform[indices, i], [slj_states[i] for i in indices]))

        # Initialise attributes of the parent object
        super().__init__(slj_states.config, states, transform)


##########################################################################
# Intermediate SLJM states
##########################################################################

class StateM(State):
    """ Class for an electron state in an intermediate coupling of SLJM states. """

    coupling = Coupling.SLJM
    is_intermediate = True

    def __init__(self, energy, values, states):
        """ Store the energy level of the state and the linear combination vector of the given SLJM states. """

        assert isinstance(energy, float)
        assert isinstance(values, np.ndarray)
        assert isinstance(states, list)
        assert len(states) == len(values)

        # Energy level of the state
        self.energy = energy

        # Sort states by decreasing weight
        weights = np.abs(values * values.conjugate())
        indices = list(reversed(np.argsort(weights)))

        # Linear combination vector
        self.values = values[indices]

        # Weight factor vector
        self.weights = weights[indices]
        assert abs(sum(self.weights) - 1.0) < 1e-7

        # The state in intermediate coupling is a linear combination of this list of related SLJM states
        self.states = [states[i] for i in indices]

        # Initialise attributes of the parent object
        super().__init__(self.states[0].config, self.states[0].mult)
        assert self.mult == 1

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
    """ Class containing a list of StateM objects representing an electron state in an intermediate SLJM coupling. """

    coupling = Coupling.SLJM
    is_intermediate = True

    def __init__(self, sljm_states, energies, transform):
        """ Store the given SLJM states, the energies of all states and the matrix containing the linear combination
        (columns) vectors of all states. """

        assert isinstance(sljm_states, StateListSLJM)
        assert len(transform.shape) == 2 and transform.shape[0] == transform.shape[1]
        assert len(energies) == len(sljm_states) == transform.shape[0]

        # List of SLJM states
        self.sljm_states = sljm_states

        # List of energy levels
        self.energies = energies

        # Build list of StateJ objects
        states = []
        for i in range(len(sljm_states)):
            values = transform[:, i]
            sljm_states = list(sljm_states.states)
            states.append(StateM(energies[i], values, sljm_states))

        # Initialise attributes of the parent object
        super().__init__(sljm_states.config, states, transform)
        assert set(self.mult) == {1}


##########################################################################
# HDF5 cache interface
##########################################################################

def init_states(config):
    data = read_transform(config)

    Product_States = StateListProduct(config, data["rowStates"], data["electronPool"])

    values = data["colStates"]
    chain = data["tensorChain"]
    transform = data["transform"]
    reprs = data["irreducibleRepresentations"]
    SLJM_states = StateListSLJM(config, values, chain, reprs, transform)

    indices = SLJM_states.stretched()
    assert chain[-1] == "Jz"
    values = values[indices, :-1]
    chain = chain[:-1]
    del reprs["Jz"]
    transform = transform[:, indices]
    SLJ_states = StateListSLJ(config, values, chain, reprs, transform)

    # Return StateList dictionary
    return {
        Coupling.Product.name: Product_States,
        Coupling.SLJM.name: SLJM_states,
        Coupling.SLJ.name: SLJ_states,
    }
