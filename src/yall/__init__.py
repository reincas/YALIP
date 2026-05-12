##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

from .state import Coupling
from .lanthanide import LANTHANIDES, RADIAL, JUDD_OFELT, Lanthanide, CONST_e, CONST_eps0, CONST_me, CONST_h, CONST_c

# def create_all():
#     for num in (1, 13, 2, 12, 3, 11, 4, 10, 5, 9, 6, 8, 7):
#         with Lanthanide(num) as ion:
#             print(ion)
#             ion.line_reduced()
#
# ameli.update()