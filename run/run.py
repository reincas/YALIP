import logging

import numpy as np
from scidatacontainer import Container

from yall import Lanthanide
from yall.state import Coupling, init_states
from yall.ameli import update, AMELI_PATH

logger = logging.getLogger("run")


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

    radial = {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
              "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67}

    num = 2
    config = f"f{num}"
    update(config)
    coupling = Coupling.SLJM

    with Lanthanide(num, coupling, radial) as ion:
        print(ion)
        for state in ion.str_levels(min_weight=0.05):
            print(state)
