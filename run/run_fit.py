##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging

import numpy as np

from yall import RADIAL, MATERIAL, LANTHANIDES, Levels, Coupling
from yall.states import get_states
from yall.fit import LevelFit, format_params, str_compare

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


def level_fit(base_states, radial, lines, opt_names):
    if base_states.coupling != Coupling.SLJ:
        raise NotImplementedError(f"Energy level fit not yet implemented for {base_states.coupling} coupling!")

    radial = radial.copy()

    mult = np.array(base_states.mult)
    matrices = {name: base_states.matrix(name) for name in radial.keys() if name != "base"}

    if not isinstance(opt_names[0], (list, tuple)):
        opt_names = [opt_names]

    opt = LevelFit({}, matrices, mult, lines)
    for i, names in enumerate(opt_names):
        raw_names = [n[1:] if n.startswith(":") else n for n in names]
        assert len(set(raw_names)) == len(raw_names)
        opt.params = {n: radial[n] for n in raw_names}

        p = format_params(opt.params, 6)
        dk = opt.get_sigma()
        logger.info(f"Stage {i}: Initial dk: {dk:.2f}, parameters: {p}")

        names = [n for n in names if n != "base" and not n.startswith(":")]
        opt.fit(names)
        radial |= opt.params

        p = format_params(opt.params, 6)
        dk = opt.get_sigma()
        logger.info(f"Stage {i}: Final dk: {dk:.2f}, parameters: {p}")

    return opt.params


if __name__ == '__main__':
    init_logger(level=logging.DEBUG)

    name = "Pr3+"
    radial = RADIAL[name]
    lines = KMEAS
    material = MATERIAL["Pb:ZBLAN"]

    num = LANTHANIDES.index(name[:2])
    config = f"f{num}"
    states = get_states(config, Coupling.SLJ)

    opt_params = level_fit(states, radial, lines, STAGE)
    ion = Levels(config, opt_params, material)
    for line in str_compare(lines, ion.states):
        logger.debug(line)
