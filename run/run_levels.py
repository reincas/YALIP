##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np

from yall import MATERIAL, Coupling, Levels

logger = logging.getLogger("run_levels")


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


if __name__ == "__main__":
    init_logger(level=logging.DEBUG)

    material = MATERIAL["Pb:ZBLAN"]
    radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
              "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}
    jo = {"JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972}

    num = 2
    config = f"f{num}"
    coupling = Coupling.SLJ

    ion = Levels(config, coupling, radial, jo, material)
    logger.info("List of states in intermediate coupling:")
    for state in ion.str_levels(min_weight=0.05):
        logger.info(f"    {state}")

    reduced = ion.dipole
    R = np.column_stack((reduced.U2[1:, 0], reduced.U4[1:, 0], reduced.U6[1:, 0], reduced.LS[1:, 0]))
    logger.info("Squared reduced matrix elements (U2, U4, U6, LS):")
    for line in R:
        line = "  ".join([f"{v:7.4f}" for v in line])
        logger.info(f"    {line}")

    f = ion.oscillator_strengths()
    f = np.column_stack((f.ed[1:, 0], f.md[1:, 0])) * 1e8
    logger.info("GSA oscillator strengths (ed, md) in 1e-8:")
    for line in f:
        line = "  ".join([f"{v:7.1f}" for v in line])
        logger.info(f"    {line}")

    A = ion.radiative_rates()
    i = -1
    A = np.column_stack((A.ed[i - 1::-1, i], A.md[i - 1::-1, i]))
    logger.info("Radiative emission rates (ed, md) to ground state in 1/s:")
    for line in A:
        line = "  ".join([f"{v:7.0f}" for v in line])
        logger.info(f"    {line}")

    t = ion.life_times()[1:]
    logger.info("Radiative life times in ms:")
    line = "  ".join([f"{v:.2f}" for v in t * 1000])
    logger.info(f"    {line}")

    beta = ion.branching_ratios()
    i = -1
    beta = beta[i - 1::-1, i]
    logger.info("Branching ratios:")
    line = "  ".join([f"{v:.3f}" for v in beta])
    logger.info(f"    {line}")
