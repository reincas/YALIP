##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides the functions wigner3j and wigner6j calculating
# Wigner 3-j and 6-j symbols for integer (int) or half-integer (HalfInt)
# arguments.
#
##########################################################################

import math
from typing import Union
from functools import lru_cache

from .halfint import HalfInt


##########################################################################
# Wigner 3-j symbol
##########################################################################

def wigner3j(j1: Union[int, HalfInt], j2: Union[int, HalfInt], j3: Union[int, HalfInt],
             m1: Union[int, HalfInt], m2: Union[int, HalfInt], m3: Union[int, HalfInt]) -> float:
    """ Calculate Wigner 3-j symbol for integer or half-integer (HalfInt) parameters. """

    return wigner3j_double(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3)


@lru_cache(maxsize=None)
def wigner3j_double(dj1: int, dj2: int, dj3: int, dm1: int, dm2: int, dm3: int) -> float:
    """ Calculate Wigner 3-j symbol for doubled (integer) values. """

    ##############################################
    # Use shortcut evaluation if all mi == 0
    if dm1 == dm2 == dm3 == 0:

        # Validate integer condition for all jis
        if not (dj1 % 2 == dj2 % 2 == dj3 % 2 == 0):
            return 0.0

        # Validate even sum of the jis
        j = (dj1 + dj2 + dj3) >> 1
        if j % 2 != 0:
            return 0.0

        # Factorial arguments as single values
        t1 = j - dj1
        t2 = j - dj2
        t3 = j - dj3
        t4 = j + 1

        # Validate triangle condition for the jis
        if t1 < 0 or t2 < 0 or t3 < 0:
            return 0.0

        # Calculate first factor
        result = math.sqrt(math.factorial(t1) * math.factorial(t2) * math.factorial(t3) / math.factorial(t4))

        # Half values as factorial arguments for the second factor
        ht1 = t1 >> 1
        ht2 = t2 >> 1
        ht3 = t3 >> 1
        ht4 = t4 >> 1

        # Multiply second factor
        result *= math.factorial(ht4) / (math.factorial(ht3) * math.factorial(ht2) * math.factorial(ht1))

        # Convert result to float and apply sign factor
        result = float(result)
        if ht4 % 2:
            result = -result
        return result

    ##############################################
    # Full evaluation of the Wigner 3-j symbol

    # Validate that the sum all mis is zero
    if dm1 + dm2 + dm3 != 0:
        return 0.0

    # Validate triangular condition for jis
    dt1 = dj2 + dj3 - dj1
    dt2 = dj3 + dj1 - dj2
    dt3 = dj1 + dj2 - dj3
    if dt1 < 0 or dt2 < 0 or dt3 < 0:
        return 0.0

    # Validate mi <= ji
    dt1m = dj1 - dm1
    dt2m = dj2 - dm2
    dt3m = dj3 - dm3
    if dt1m < 0 or dt2m < 0 or dt3m < 0:
        return 0.0

    # Validate mi >= -ji
    dt1p = dj1 + dm1
    dt2p = dj2 + dm2
    dt3p = dj3 + dm3
    if dt1p < 0 or dt2p < 0 or dt3p < 0:
        return 0.0

    # Validate that ji and mi are either both integers
    # or both half-integers
    if not (dt1p % 2 == dt2p % 2 == dt3p % 2 == 0):
        return 0.0

    # Validate that the sum of jis is integer
    dt = dj1 + dj2 + dj3 + 2
    if dt % 2 != 0:
        return 0.0

    # Factorial arguments as single values
    t1 = dt1 >> 1
    t2 = dt2 >> 1
    t3 = dt3 >> 1
    t1m = dt1m >> 1
    t2m = dt2m >> 1
    t3m = dt3m >> 1
    t1p = dt1p >> 1
    t2p = dt2p >> 1
    t3p = dt3p >> 1
    t = dt >> 1

    # Calculate common factor
    factor = math.sqrt(math.factorial(t1) * math.factorial(t2) * math.factorial(t3) *
                       math.factorial(t1m) * math.factorial(t2m) * math.factorial(t3m) *
                       math.factorial(t1p) * math.factorial(t2p) * math.factorial(t3p) /
                       math.factorial(t))

    # Boundaries of the sum factor
    s1 = -(dj3 - dj2 + dm1) // 2
    s2 = -(dj3 - dj1 - dm2) // 2
    zmin = max(0, s1, s2)
    zmax = min(t1m, t2p, t3)

    # Calculate the sum factor
    values = []
    for z in range(zmin, zmax + 1):

        # Factorial arguments as single values
        tm1 = t1m - z
        tm2 = t2p - z
        tm3 = t3 - z
        tp1 = z - s1
        tp2 = z - s2

        # Calculate summand
        value = 1.0 / (math.factorial(tm1) * math.factorial(tm2) * math.factorial(tm3) *
                       math.factorial(tp1) * math.factorial(tp2) * math.factorial(z))

        # Apply sign factor and append to the list
        if z % 2:
            value = -value
        values.append(value)

    # Multiply common factor and apply common sign factor
    result = factor * sum(sorted(values))
    if ((dj1 - dj2 - dm3) >> 1) % 2:
        result = -result
    return result


##########################################################################
# Wigner 6-j symbol
##########################################################################

def wigner6j(j1: Union[int, HalfInt], j2: Union[int, HalfInt], j3: Union[int, HalfInt],
             l1: Union[int, HalfInt], l2: Union[int, HalfInt], l3: Union[int, HalfInt]) -> float:
    """ Calculate Wigner 6-j symbol for integer or half-integer (HalfInt) parameters. """

    return wigner6j_double(2 * j1, 2 * j2, 2 * j3, 2 * l1, 2 * l2, 2 * l3)


@lru_cache(maxsize=None)
def wigner6j_double(dj1: int, dj2: int, dj3: int, dl1: int, dl2: int, dl3: int) -> float:
    """ Calculate Wigner 6-j symbol for doubled (integer) values. """

    # Common factor
    factor = 1.0
    factor *= triangle(dj1, dj2, dj3)
    if factor <= 0.0:
        return 0.0
    factor *= triangle(dl1, dl2, dj3)
    if factor <= 0.0:
        return 0.0
    factor *= triangle(dl1, dj2, dl3)
    if factor <= 0.0:
        return 0.0
    factor *= triangle(dj1, dl2, dl3)
    if factor <= 0.0:
        return (0.0)

    # Determine the minimum z value based on the factorial arguments
    arg1 = (dj1 + dj2 + dj3) >> 1
    arg2 = (dj1 + dl2 + dl3) >> 1
    arg3 = (dl1 + dj2 + dl3) >> 1
    arg4 = (dl1 + dl2 + dj3) >> 1
    zmin = max(arg1, arg2, arg3, arg4)

    # Determine the maximum z value based on the factorial arguments
    arg5 = (dj1 + dj2 + dl1 + dl2) >> 1
    arg6 = (dj2 + dj3 + dl2 + dl3) >> 1
    arg7 = (dj3 + dj1 + dl3 + dl1) >> 1
    zmax = min(arg5, arg6, arg7)

    # Calculate sum
    values = []
    for z in range(zmin, zmax + 1):
        value = math.factorial(z + 1) / (math.factorial(z - arg1) * math.factorial(z - arg2) *
                                         math.factorial(z - arg3) * math.factorial(z - arg4) *
                                         math.factorial(arg5 - z) * math.factorial(arg6 - z) *
                                         math.factorial(arg7 - z))
        if z % 2 != 0:
            value = -value
        values.append(value)

    # Return result
    return sum(sorted(values)) * math.sqrt(factor)


def triangle(a, b, c):
    """ Triangle function for the Wigner 6-j symbol. """

    # Factorial arguments as double values
    arg1 = a + b - c
    arg2 = a - b + c
    arg3 = -a + b + c
    arg4 = a + b + c + 2

    # No argument may be negative
    if arg1 < 0 or arg2 < 0 or arg3 < 0:
        return 0.0

    # Either none or two of the arguments may be half-integer
    if arg4 % 2 != 0:
        return 0.0

    # Factorial arguments as single values
    arg1 >>= 1
    arg2 >>= 1
    arg3 >>= 1
    arg4 >>= 1

    # Return value
    return math.factorial(arg1) * math.factorial(arg2) * math.factorial(arg3) / math.factorial(arg4)
