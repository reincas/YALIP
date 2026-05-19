##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging

from yall import RADIAL, MATERIAL, LANTHANIDES, Fit

logger = logging.getLogger("run_fit")

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

STAGES = [
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

    name = "Pr3+"
    radial = RADIAL[name]
    lines = KMEAS

    jo = None
    material = MATERIAL["Pb:ZBLAN"]

    num = LANTHANIDES.index(name[:2])
    config = f"f{num}"

    opt = Fit(config, radial, jo, material)
    opt.level_fit(lines, STAGES)
    for line in opt.str_compare():
        logger.debug(line)
