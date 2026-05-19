##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import numpy as np
from typing import cast, Any

from spectrum import jo_factors
from yall import MATERIAL, Levels

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

    ion = Levels(config, radial, jo, material)
    for state in ion.str_levels(min_weight=0.05):
        print(state)

    reduced = ion.dipole
    R = np.column_stack((reduced.U2[1:, 0], reduced.U4[1:, 0], reduced.U6[1:, 0], reduced.LS[1:, 0]))
    np.set_printoptions(formatter=cast(Any, {'float': '{:7.4f}'.format}), linewidth=120)
    print(R)

    f = ion.oscillator_strengths()
    f = np.column_stack((f.ed[1:, 0], f.md[1:, 0])) * 1e8
    np.set_printoptions(formatter=cast(Any, {'float': '{:7.1f}'.format}), linewidth=120)
    print(f)

    A = ion.radiative_rates()
    A = np.column_stack((A.ed[-2::-1, -1], A.md[-2::-1, -1]))
    np.set_printoptions(formatter=cast(Any, {'float': '{:7.0f}'.format}), linewidth=120)
    print(A)

    print("**********")
    x = np.array([ion.judd_ofelt[f"JO/{lam}"] for lam in (2, 4, 6)])
    factor_ed, factor_md = jo_factors(ion.mult[0], ion.energies, material)
    A = np.column_stack((reduced.U2[:, 0], reduced.U4[:, 0], reduced.U6[:, 0])) * factor_ed[:]
    b = reduced.LS[1:, 0] * factor_md
    x = np.array([ion.judd_ofelt[f"JO/{lam}"] for lam in (2, 4, 6)])
    f = np.column_stack((A @ x, b)) * 1e8
    np.set_printoptions(formatter=cast(Any, {'float': '{:7.1f}'.format}), linewidth=120)
    print(f[1:,0])

