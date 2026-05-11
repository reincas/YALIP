##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This Module provides the class HalfInt which represents half-integer
# values as needed for angular momentum quantum numbers. Basic mathematical
# operations and comparisons are supported.
#
##########################################################################

class HalfInt:
    def __init__(self, numerator: int):
        assert isinstance(numerator, int)
        self.numerator = numerator

    def __neg__(self):
        numerator = -self.numerator
        if numerator % 2:
            return HalfInt(numerator)
        return numerator // 2

    def __add__(self, other):
        if isinstance(other, HalfInt):
            numerator = self.numerator + other.numerator
        elif isinstance(other, int):
            numerator = self.numerator + 2 * other
        else:
            return NotImplemented
        if numerator % 2:
            return HalfInt(numerator)
        return numerator // 2

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, HalfInt):
            numerator = self.numerator - other.numerator
        elif isinstance(other, int):
            numerator = self.numerator - 2 * other
        else:
            return NotImplemented
        if numerator % 2:
            return HalfInt(numerator)
        return numerator // 2

    def __rsub__(self, other):
        if isinstance(other, int):
            numerator = 2 * other - self.numerator
        else:
            return NotImplemented
        if numerator % 2:
            return HalfInt(numerator)
        return numerator // 2

    def __mul__(self, other):
        if other == 2:
            return self.numerator
        if isinstance(other, int):
            numerator = self.numerator * other
        else:
            return NotImplemented
        if numerator % 2:
            return HalfInt(numerator)
        return numerator // 2

    def __rmul__(self, other):
        return self.__mul__(other)

    def __lt__(self, other):
        if isinstance(other, HalfInt):
            return self.numerator < other.numerator
        elif isinstance(other, int):
            return self.numerator < 2 * other
        return NotImplemented

    def __le__(self, other):
        if isinstance(other, HalfInt):
            return self.numerator <= other.numerator
        elif isinstance(other, int):
            return self.numerator <= 2 * other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, HalfInt):
            return self.numerator > other.numerator
        elif isinstance(other, int):
            return self.numerator > 2 * other
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, HalfInt):
            return self.numerator >= other.numerator
        elif isinstance(other, int):
            return self.numerator >= 2 * other
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, HalfInt):
            return self.numerator == other.numerator
        elif isinstance(other, int):
            return self.numerator == 2 * other
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return f"HalfInt({self.numerator})"

    def __str__(self):
        return f"{self.numerator}/2"
