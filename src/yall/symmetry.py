##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides Symmetry classes for symmetry representations of
# electron states in SLJM or SLJ coupling. They offer a uniform interface
# which is mainly used to generate string representations of these states.
#
# The class SymmetryList holds a list of equivalent Symmetry classes. It
# is used in the states module to determine the transformation matrices
# from the determinantal product states to SLJM or SLJ states.
#
##########################################################################

import math

#from .halfint import HalfInt

# Letters used to represent orbital angular momenta
#CHR_ORBITAL = "SPDFGHIKLMNOQRTUVWXYZ"


# def casimir_Rk(w: tuple, d: int) -> int:
#     """ This function returns the eigenvalue of the Casimir operator of the rotational group Rd in d dimensions for
#     the given representation w. The representation is a tuple of integers. Its length is given by the rank of the
#     operator which is d/2 for even d and (d-1)/2 for odd d. """
#
#     assert len(w) == d // 2
#
#     sum = 0
#     for i in range(len(w)):
#         sum += w[i] * (w[i] + d - 2 - 2 * i)
#     return sum // 2
#
#
# def casimir_G2(u: tuple):
#     """ This function returns the integer eigenvalue of the Casimir operator of the special group G2 for the given
#     representation u. Due to the rank 2 of the operator, the representation is a 2-tuple of integers. """
#
#     return u[0] * u[0] + u[1] * u[1] + u[0] * u[1] + 5 * u[0] + 4 * u[1]
#
#
#
# def casimir_dict():
#     """ Return a dictionary mapping eigenvalues of the Casimir operators G(R5), G(R7) and G(G2) to their
#     representation tuples of integers. """
#
#     # Casimir operator of the rotational group in 5 dimensions R5
#     R5 = {}
#     for i in range(3):
#         for j in range(i + 1):
#             w = (i, j)
#             R5[casimir_Rk(w, 5)] = w
#
#     # Casimir operator of the rotational group in 7 dimensions R7
#     R7 = {}
#     for i in range(3):
#         for j in range(i + 1):
#             for k in range(j + 1):
#                 w = (i, j, k)
#                 R7[casimir_Rk(w, 7)] = w
#
#     # Casimir operator of the special group G2
#     G2 = {}
#     for i in range(5):
#         for j in range(i + 1):
#             u = (i, j)
#             G2[casimir_G2(u)] = u
#
#     # Return mapping tuple
#     return {"GR/5": R5, "GR/7": R7, "GG/2": G2}
#
#
# # Store Casimir mapping tuple as static module constant
# CASIMIR = casimir_dict()


class Symmetry:
    """ Abstract class for the representation of the state of a certain symmetry defined by its subclasses. Common
    attributes are the operator name (self.name), the representation symbol (self.symbol), and it irreducible
    representation (self.repr). """

    name: str
    symbol: str
    # value: float
    # key: int
    repr: str

    def __str__(self):
        """ Return the irreducible representation of the symmetry state. """

        return self.repr

    # def __lt__(self, other):
    #     """ Return True if the symmetry state key is strictly less than that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key < other.key
    #
    # def __le__(self, other):
    #     """ Return True if the symmetry state key is less than or equal to that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key <= other.key
    #
    # def __gt__(self, other):
    #     """ Return True if the symmetry state key is strictly greater than that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key > other.key
    #
    # def __ge__(self, other):
    #     """ Return True if the symmetry state key is greater than or equal to that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key >= other.key
    #
    # def __eq__(self, other):
    #     """ Return True if the symmetry state key is equal to that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key == other.key
    #
    # def __ne__(self, other):
    #     """ Return True if the symmetry state key is not equal to that of the given symmetry state. """
    #
    #     if type(self) != type(other):
    #         return NotImplemented
    #     return self.key != other.key


class SymmetryS2(Symmetry):
    """ Symmetry class of the symmetry group equivalent to the squared total electron spin. For electrons with
    the orbital quantum number l this is the group of unitary transformations in 2l+1 dimensions, thus U7 for
    f electrons. The multiplicity 2S+1 is used as representation. """

    name = "S2"
    symbol = "2S+1"

    def __init__(self, repr):
        self.repr = repr
        # self.value = value
        # self.key = self._to_int_(math.sqrt(4 * self.value + 1))
        # self.str_value = (str(self.key))
        #
        # if self.key % 2 == 0:
        #     self.S = HalfInt(self.key - 1)
        # else:
        #     self.S = (self.key - 1) // 2


class SymmetryGR7(Symmetry):
    """ Symmetry class of the rotational group in 7 dimensions R7 or SO(7). We use the Casimir operator G(R7) as
    equivalent tensor operator and its 3-tuple W as representation. """

    name = "GR/7"
    symbol = "W"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(5 * self.value)
        self.str_value = "(%d%d%d)" % CASIMIR["GR/7"][self.key]


class SymmetryGR5(Symmetry):
    """ Symmetry class of the rotational group in 5 dimensions R5 or SO(5). We use the Casimir operator G(R5) as
    equivalent tensor operator and its 2-tuple W as representation. """

    name = "GR/5"
    symbol = "W"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(3 * self.value)
        self.str_value = "(%d%d)" % CASIMIR["GR/5"][self.key]


class SymmetryGG2(Symmetry):
    """ Symmetry class of the special group G2. We use the Casimir operator G(G2) as equivalent tensor operator and
    its 2-tuple W as representation. """

    name = "GG/2"
    symbol = "U"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(12 * self.value)
        self.str_value = "[%d%d]" % CASIMIR["GG/2"][self.key]


class SymmetryL2(Symmetry):
    """ Symmetry class of the rotational group in 3 dimensions R3 or SO(3). We use the tensor operator of the squared
    total orbital angular momentum L2 as equivalent tensor operator and its uppercase letter as representation. """

    name = "L2"
    symbol = "L"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_((math.sqrt(4 * self.value + 1) - 1) / 2)
        self.str_value = CHR_ORBITAL[self.key]
        self.L = self.key

class SymmetryJ2(Symmetry):
    """ Symmetry class of the tensor operator of the squared total angular momentum J2 and its integer or half-integer
    value as representation. """

    name = "J2"
    symbol = "J"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(math.sqrt(4 * self.value + 1) - 1)
        self.str_value = f"{self.key // 2}" if self.key % 2 == 0 else f"{self.key}/2"

        if self.key % 2 == 0:
            self.J = self.key // 2
        else:
            self.J = HalfInt(self.key)


class SymmetryJz(Symmetry):
    """ Symmetry class of the component q=0 of the tensor operator of the total angular momentum Jz and its integer
    or half-integer value as representation. """

    name = "Jz"
    symbol = "M"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(2 * self.value)
        self.str_value = f"{(self.key // 2):+d}" if self.key % 2 == 0 else f"{self.key:+d}/2"

        if self.key % 2 == 0:
            self.M = self.key // 2
        else:
            self.M = HalfInt(self.key)

class SymmetryTau(Symmetry):
    """ This pseudo Symmetry class is used for the arbitrarily labeled pairs of states which match in all real
    symmetry representations in the chain of operators S2, G(R7), G(G2), L2, and J2. The representation is either
    "a" or "b" for these matching states or "*" for unique states. """

    name = "tau"
    symbol = "𝜏"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(self.value)
        self.str_value = "*ab"[self.key]


class SymmetryNum(Symmetry):
    """ This pseudo Symmetry class is used as a short-cut for the full representation of different states, which
    match in S2, L2, and J2. The representation is an integer number. """

    name = "num"
    symbol = "#"

    def __init__(self, value):
        self.value = value
        self.key = self._to_int_(self.value)
        self.str_value = str(self.key)


# Static dictionary mapping operator names to the respective classes
SYMMETRY = {
    "S2": SymmetryS2,
    "GR/7": SymmetryGR7,
    "GR/5": SymmetryGR5,
    "GG/2": SymmetryGG2,
    "L2": SymmetryL2,
    "J2": SymmetryJ2,
    "Jz": SymmetryJz,
    "tau": SymmetryTau,
    "num": SymmetryNum,
}


class SymmetryList:
    """ This class keeps a list of equivalent Symmetry objects. """

    def __init__(self, values, name):
        """ Convert the given vector of eigenvalues, representing states, to a list of Symmetry objects with the
        given operator name. """

        if not name in SYMMETRY:
            raise ValueError(f"Unknown operator {name}")

        self.name = name
        self.values = list(values)
        self.object_class = SYMMETRY[name]
        self.objects = [self.object_class(value) for value in values]
        self.keys = [sym.key for sym in self.objects]

    def count(self):
        """ Return a dictionary which contains the number of all states with the same symmetry representation. """

        keys = sorted(set([obj.key for obj in self.objects]))
        return dict([(self.find(key).str_value, self.count_key(key)) for key in keys])

    def find(self, key):
        """ Find and return one Symmetry object with the given key value. """

        for obj in self.objects:
            if obj.key == key:
                return obj
        raise ValueError(f"Unknown symmetry key {key}!")

    def count_key(self, key: int) -> int:
        """ Return the number of object with the given key value. """

        return len([obj for obj in self.objects if obj.key == key])

    def split_syms(self, last_slices=None):
        """ Take an optional list of consecutive slices covering all states. Default is a single slice. Return a
        new list of consecutive slices, which is generated by splitting each slice in last_slices in such a way, that
        each of the new slices contains only Symmetry objects with the same key value. """

        # Default is one slice covering all states
        if not last_slices:
            last_slices = [slice(0, len(self))]

        # Make sure that the slices cover all states and are consecutive
        assert last_slices[0].start == 0 and last_slices[-1].stop == len(self)
        for i in range(1, len(last_slices)):
            assert last_slices[i - 1].stop == last_slices[i].start

        # Build new list of slices by splitting last_slices
        new_slices = []
        i = 0
        for last_slice in last_slices:

            # Find splitting index
            for j in range(last_slice.start + 1, last_slice.stop + 1):

                # Each of the new slices covers only Symmetry objects with the same key value
                if j == last_slice.stop or self.objects[j].key != self.objects[i].key:
                    new_slices.append(slice(i, j))
                    i = j

            # Start index of the next slice
            i = last_slice.stop

        # Return new list of slices
        return new_slices

    def __getitem__(self, item: int):
        """ Return Symmetry object addressed by the given state index. """

        return self.objects[item]

    def __add__(self, other):
        """ Use the operator "+" to concatenate two SymmetryList objects. """

        # Other list is empty
        if not other:
            return self

        # Other object must also be a SymmetryList
        if not isinstance(other, SymmetryList):
            return NotImplemented

        # The operator name of both symmetry lists must be the same
        if self.name != other.name:
            raise ValueError(f"Cannot add lists of different symmetries!")

        # Return SymmetryList with concatenated lists of Symmetry objects
        return SymmetryList(self.values + other.values, self.name)

    def __radd__(self, other):
        """ Use the operator "+" to concatenate two SymmetryList objects. """

        return self.__add__(other)

    def __len__(self):
        """ Return the number of Symmetry objects, which is equivalent to the number of states. """

        return len(self.objects)

    def __str__(self):
        """ Return a string representation showing the number of states for every key. """

        count = self.count()
        return ", ".join(f"{self.find(key)}: {count[key]}" for key in sorted(count.keys()))
