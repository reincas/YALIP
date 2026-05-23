##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

from enum import Enum


class Coupling(Enum):
    Product = 0
    SLJM = 1
    SLJ = 2


from .lanthanide import LANTHANIDES, RADIAL, JUDD_OFELT, MATERIAL
from .spectrum import CONST_e, CONST_eps0, CONST_me, CONST_h, CONST_c, CONST_gs, Cauchy, Sellmeier
from .states import States
from .levels import Levels
from .fits import Fit
