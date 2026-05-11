##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

from .halfint import HalfInt
from .wigner import wigner3j, wigner6j
from .single import SingleElements, single_elements, init_single
from .unit import ORBITAL, SPIN, MAGNETIC, LEN_SHELL, product_states, calc_unit, get_unit
from .symmetry import SYMMETRY, Symmetry, SymmetryS2, SymmetryGR7, SymmetryGR5, SymmetryGG2, SymmetryL2, \
    SymmetryJ2, SymmetryJz, SymmetryTau, SymmetryNum, SymmetryList
from .state import Coupling, StateListProduct, StateListSLJM, StateListSLJ, StateListJ, StateListM, \
    StateProduct, StateSLJM, StateSLJ, StateJ, StateM, init_states
from .matrix import normalise_radial, build_hamilton, reduced_matrix, Matrix, get_matrix
from .lanthanide import LANTHANIDES, RADIAL, JUDD_OFELT, Lanthanide, CONST_e, CONST_eps0, CONST_me, CONST_h, CONST_c


def create_all():
    for num in (1, 13, 2, 12, 3, 11, 4, 10, 5, 9, 6, 8, 7):
        with Lanthanide(num) as ion:
            print(ion)
            ion.line_reduced()
