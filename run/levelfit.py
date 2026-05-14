##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np

from yall import Coupling, Lanthanide, RADIAL
from yall.fit import str_compare, format_params

logger = logging.getLogger("levelfit")

KMEAS = [
    [1, '3H_5', 2365, 4],
    [2, '3H_6', 4485, 18],
    [3, '3F_2', 5105, 5],
    [4, '3F_3', 6467, 7],
    [5, '3F_4', 6958, 10],
    [6, '1G_4', 9883, 7],
    [7, '1D_2', 17026, 20],
    [8, '3P_0', 20859, 12],
    [(9, 10), ('3P_1', '1I_6'), 21505, 20],
    [11, '3P_2', 22645, 10]
]

STAGE = [
    ["base", "H1/2", "H1/4", "H1/6", "H2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2"],
    ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2", "H5fix", "H6fix"],
]


def init_logger(file_name=None, level=logging.INFO):
    root = logging.getLogger()
    root.setLevel(level)
    log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    if file_name is not None:
        file_h = logging.FileHandler(file_name, mode="a")
        file_h.setFormatter(log_format)
        file_h.setLevel(level)
        root.addHandler(file_h)

    console_h = logging.StreamHandler()
    console_h.setFormatter(log_format)
    console_h.setLevel(level)
    root.addHandler(console_h)


if __name__ == '__main__':
    init_logger(level=logging.DEBUG)

    num = 2
    coupling = Coupling.SLJ
    radial = RADIAL["Pr3+"]
    # radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
    #     "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
    lines = KMEAS

    with Lanthanide(num, coupling, radial) as ion:

        opt_params = ion.level_fit(lines, STAGE, radial)
        for line in ion.str_compare_lines(lines):
            logger.debug(line)
