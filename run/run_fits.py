##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging

from yalip import RADIAL, MATERIAL, Coupling, Fits, ion2config

logger = logging.getLogger("run_fit")

MEAS = {
    "Pr3+": [
        [1,      '3H_5',             2365,  4, 196.2e-8,  2.8e-8],
        [2,      '3H_6',             4485, 18,  77.4e-8,  5.1e-8],
        [3,      '3F_2',             5105,  5, 283.3e-8,  6.2e-8],
        [4,      '3F_3',             6467,  7, 637.0e-8, 16.2e-8],
        [5,      '3F_4',             6958, 10, 244.2e-8, 13.1e-8],
        [6,      '1G_4',             9883,  7,  26.4e-8,  0.5e-8],
        [7,      '1D_2',            17026, 20, 199.2e-8, 10.7e-8],
        [8,      '3P_0',            20859, 12, 180.0e-8, 17.4e-8],
        [(9, 10), ('3P_1', '1I_6'), 21505, 20, 542.6e-8, 30.5e-8],
        [11,      '3P_2',           22645, 10, 926.8e-8, 28.7e-8],
    ],
    "Er3+": [
        [1,  "4I_13/2",   6585,   3, 186.3e-8,  2.9e-8],
        [2,  "4I_11/2",  10246,   8,  55.6e-8,  3.2e-8],
        [3,  "4I_9/2",   12480,   6,  28.7e-8,  0.8e-8],
        [4,  "4F_9/2",   15322,  17, 185.6e-8, 13.5e-8],
        [5,  "4S_3/2",   18454,  21,  49.0e-8,  5.7e-8],
        [6,  "2H2_11/2", 19209,   7, 475.2e-8, 14.3e-8],
        [7,  "4F_7/2",   20538,   5, 155.3e-8,  3.9e-8],
        [8,  "4F_5/2",   22224,  11,  49.1e-8,  2.7e-8],
        [9,  "4F_3/2",   22619,  14,  25.4e-8,  2.0e-8],
        [10, "4F_9/2",   24618,  11,  56.8e-8,  3.2e-8],
        [11, "4G_11/2",  26476,  16, 895.3e-8, 59.0e-8],
        [12, "4G_9/2",   27480,  39, 148.8e-8, 35.5e-8],
        [13, "2K_15/2",  27894, 478,  40.5e-8, 31.3e-8],
        [14, "4G_7/2",   28149, 588,  30.5e-8, 30.5e-8],
        [15, "2P_3/2",   31608,  46,   6.9e-8,  1.3e-8],
        [17, "2K_13/2",  33096,  24,  10.5e-8,  0.9e-8],
        [19, "4G_7/2",   34111,  15,  15.8e-8,  1.0e-8],
        [20, "2D1_5/2",  34920,  20,   9.0e-8,  1.2e-8],
        [21, "2H2_9/2",  36496,  19,  36.5e-8,  3.2e-8],
        [22, "4D_5/2",   38655, 279,  14.3e-8, 19.9e-8],
        [23, "4D_7/2",   39324,  25, 753.7e-8, 53.9e-8],
    ],
    "Tm3+": [
        [1, "3F_4",   5862,   4, 186.9e-8,  2.3e-8],
        [2, "3H_5",   8303,  14, 156.0e-8,  7.3e-8],
        [3, "3H_4",  12702,   5, 201.3e-8,  3.8e-8],
        [4, "3F_3",  14574,  18, 253.9e-8, 16.4e-8],
        [5, "3F_2",  15180, 119,  20.4e-8,  7.7e-8],
        [6, "1G_4",  21366,  42,  74.7e-8,  7.1e-8],
        [7, "1D_2",  28018,  68, 195.1e-8, 29.8e-8],
        [8, "1I_6",  35039,  26,  68.3e-8,  2.7e-8],
        [10, "3P_1", 36576,  24,  26.2e-8,  2.0e-8],
        [11, "3P_2", 38344,  16, 245.8e-8,  8.1e-8],
    ],
}

STAGES = {
    "Pr3+": [
        ["base", "H1/2", "H1/4", "H1/6", "H2"],
        ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2"],
        ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2", "H5fix", "H6fix"],
    ],
    "Er3+": [
        ["base", "H1/2", "H1/4", "H1/6", "H2"],
        ["base", "H1/2", "H1/4", "H1/6", "H2",
         "H3/0", "H3/1", "H3/2", "H4/2", "H4/3", "H4/4", "H4/6", "H4/7", "H4/8"],
        ["base", "H1/2", "H1/4", "H1/6", "H2",
         "H3/0", "H3/1", "H3/2", "H4/2", "H4/3", "H4/4", "H4/6", "H4/7", "H4/8",
         "H5fix", "H6fix"],
    ],
    "Tm3+": [
        ["base", "H1/2", "H1/4", "H1/6", "H2"],
        ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2"],
        ["base", "H1/2", "H1/4", "H1/6", "H2", "H3/0", "H3/1", ":H3/2", "H5fix", "H6fix"],
    ],
}


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

    name = "Er3+"

    config = ion2config(name)
    coupling = Coupling.SLJ
    radial = RADIAL[name]
    material = MATERIAL["Pb:ZBLAN"]
    opt = Fits(config, coupling, radial, material)

    lines = MEAS[name]
    stages = STAGES[name]
    opt.run(lines, stages)

    for line in opt.table():
        print(line)
