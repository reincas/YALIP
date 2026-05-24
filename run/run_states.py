##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging

from yalip import Coupling, States, ion2config

logger = logging.getLogger("run_states")


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

    name = "Pr3+"
    config = ion2config(name)
    coupling = Coupling.SLJM
    states = States(config, coupling)
    logger.info(states)
    logger.info(f"Number of states: {len(states)}")

    matrix = states.matrix("H1/2")
    logger.info(f"Matrix shape: {matrix.shape}")

    logger.info(f"List of {coupling.name} coupling states:")
    for state in states:
        logger.info(f"    {state}")

