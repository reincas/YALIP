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

import logging
import numpy as np

from . import Coupling
from .ameli import get_ameli_transform
from .matrix import get_matrix

logger = logging.getLogger("yall.states")

SPECTRAL = "spdfghiklmnoqrtuvwxyz"


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

    def __init__(self, config, coupling, J=None, mult=1):
        """ Store electron configuration and state multiplicity. """

        assert isinstance(config, str)
        assert isinstance(coupling, Coupling)
        assert J is None or isinstance(J, str)
        assert isinstance(mult, int)

        # Electron configuration
        self.config = config

        # Coupling scheme
        self.coupling = coupling

        # Quantum number of total angular momentum as string
        self.J = J

        # State multiplicity
        self.mult = mult

    @property
    def state_space(self):
        return self.coupling.name.lower()


class StateList:
    """ Abstract class for a list of electron states in a certain coupling scheme. """

    def __init__(self, config, coupling, states, transform=None):
        """ Store electron configuration, list of state objects, transformations matrix and list of state
        multiplicities. """

        assert isinstance(config, str)
        assert isinstance(coupling, Coupling)
        assert isinstance(states, list)

        # Electron configuration
        self.config = config

        # Coupling scheme
        self.coupling = coupling

        # List of state objects
        self.states = states

        # Transformation matrix
        self.transform = transform

        # Values of total angular momentum as list of strings
        self.J = [state.J for state in self.states]

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
        """ Generate the short names of all states. """

        for state in self.states:
            yield state.short()

    def long(self):
        """ Generate the long names of all states. """

        for state in self.states:
            yield state.long()

    @property
    def state_space(self):
        return self.coupling.name.lower()

    def matrix(self, name):
        """ Return the matrix of the tensor operator or the weighted sum of tensor operators. """

        return get_matrix(name, self.config, self.state_space)

    def reduced(self, name):
        """ Return the reduced matrix of the given tensor operator or the weighted sum of tensor operators. """

        assert self.coupling == Coupling.SLJ
        return get_matrix(name, self.config, "slj_reduced")


##########################################################################
# Product states
##########################################################################

class StateProduct(State):
    """ Class for a determinantal product state. """

    def __init__(self, config, quantum):
        """ Store electron configuration and single electron quantum numbers. """

        # Single electron quantum numbers
        self.quantum = quantum

        # Initialise attributes of the parent object
        super().__init__(config, Coupling.Product)

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
            l = SPECTRAL[int(single["l"])]
            ml = norm_magnetic(single["ml"])
            ms = norm_magnetic(single["ms"])
            quantum.append(f"{l}_{ml}_{ms}")
        return " ".join(quantum)

    def __getitem__(self, item: int):
        """ Return the quantum number dictionary of the state with the given index. """

        return self.quantum[item]

    def __str__(self):
        """ Return a long string representation of the state. """

        return self.long()


class StateListProduct(StateList):
    """ Class containing a list of StateProduct objects representing an electron configuration. """

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
        super().__init__(config, Coupling.Product, states)


##########################################################################
# SLJM states
##########################################################################

class StateSLJM(State):
    """ Class for an electron state in SLJM coupling following the chain of symmetry operators. """

    def __init__(self, config, quantum):
        """ Store electron configuration and dictionary of quantum numbers. """

        # Quantum number dictionary
        self.quantum = quantum

        # Initialise attributes of the parent object
        super().__init__(config, Coupling.SLJM, quantum["J2"])

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

    def __init__(self, config, values, chain, repr, transform):
        assert isinstance(values, np.ndarray)
        assert len(values.shape) == 2
        assert values.shape[1] == len(chain)
        assert set(chain) == set(repr.keys()), f"{set(chain)} != {set(repr.keys())}"

        # Dictionary of indices into the lists of irreducible representations for all states
        self.values = values

        # Chain of symmetry operator names
        self.chain = chain

        # Dictionary of lists of all irreducible representations
        self.repr = repr

        # List of SLJM states
        states = [StateSLJM(config, {sym: repr[sym][s[i]] for i, sym in enumerate(chain)}) for s in values]

        # Initialise attributes of the parent object
        super().__init__(config, Coupling.SLJM, states, transform)

    def stretched(self):
        """ Return a list indices containing all stretched states. """

        return [i for i, state in enumerate(self.states) if state["Jz"].replace("+", "") == state["J2"]]


##########################################################################
# SLJ states
##########################################################################

class StateSLJ(State):
    """ Class for an electron state in SLJ coupling following the chain of symmetry operators. """

    def __init__(self, config, quantum):
        """ Store electron configuration and dictionary of quantum numbers. """

        # Quantum number dictionary
        self.quantum = quantum

        # State multiplicity
        J = quantum["J2"]
        if J.endswith("/2"):
            mult = int(J[:-2]) + 1
        else:
            mult = 2 * int(J) + 1

        # Initialise attributes of the parent object
        super().__init__(config, Coupling.SLJ, J, mult)

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

        # Initialise attributes of the parent object
        super().__init__(config, Coupling.SLJ, states, transform)


##########################################################################
# AMELI interface
##########################################################################

def get_states(config, coupling):
    """ Return StateList object for the given coupling scheme. """

    assert isinstance(config, str)
    assert coupling in (Coupling.Product, Coupling.SLJM, Coupling.SLJ)

    # Read AMELI transformation data container
    data = get_ameli_transform(config)

    # Product states
    if coupling == Coupling.Product:
        values = data["rowStates"]
        pool = data["electronPool"]
        return StateListProduct(config, values, pool)

    # SLJM states
    values = data["colStates"]
    chain = data["tensorChain"]
    transform = data["transform"]
    reprs = data["irreducibleRepresentations"]
    sljm_states = StateListSLJM(config, values, chain, reprs, transform)
    if coupling == Coupling.SLJM:
        return sljm_states

    # SLJ states
    indices = sljm_states.stretched()
    assert chain[-1] == "Jz"
    values = values[indices, :-1]
    chain = chain[:-1]
    transform = transform[:, indices]
    reprs = reprs.copy()
    del reprs["Jz"]
    return StateListSLJ(config, values, chain, reprs, transform)


class States:
    """ Class of an electron configuration in a given coupling scheme providing access to the state objects and
    matrices. """

    def __init__(self, config, coupling):
        """ Store electron configuration, coupling scheme, and a StateList object. """

        assert isinstance(config, str)
        assert isinstance(coupling, Coupling)

        self.config = config
        self.coupling = coupling
        self.states = get_states(config, coupling)

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

    def short(self):
        """ Generate the short names of all states. """

        for state in self.states.states:
            yield state.short()

    def long(self):
        """ generate the long names of all states. """

        for state in self.states.states:
            yield state.long()

    def matrix(self, name):
        """ Return the matrix of the tensor operator or the weighted sum of tensor operators. """

        return self.states.matrix(name)

    def reduced(self, name):
        """ Return the reduced matrix of the given tensor operator or the weighted sum of tensor operators. """

        return self.states.reduced(name)

    def __str__(self):
        """ Return a string representation of the object. """

        return f"States({self.config}, {self.coupling.name})"