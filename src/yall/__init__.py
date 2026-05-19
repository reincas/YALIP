##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

from enum import Enum


class Coupling(Enum):
    """ This enumeration class is used to mark the four coupling schemes used in the Lanthanide package: determinantal
    product state coupling, SLJM coupling, SLJ coupling, and intermediate coupling. """

    Product = 0
    SLJM = 1
    SLJ = 2


from .lanthanide import LANTHANIDES, RADIAL, JUDD_OFELT, MATERIAL
from .spectrum import CONST_e, CONST_eps0, CONST_me, CONST_h, CONST_c, CONST_gs
from .states import States
from .intermediate import Intermediate
