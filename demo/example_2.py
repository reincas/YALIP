##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This script generates the comparison table XI for the second application
# example in the AMELI paper.
#
##########################################################################

import logging
import math
from bisect import bisect

from yalip import Coupling, Levels, ion2config

logger = logging.getLogger("example_2")

# Literature data
DATA = {
    "Nd3+:LaF3": {
        "ref": "[Carnall1978]",
        "name": "Nd3+",
        "radial": {"H1/2": 73036, "H1/4": 52624, "H1/6": 35793,
                   "H2": 884.9, "H3/0": 21.28, "H3/1": -583, "H3/2": 1443,
                   "H4/2": 306, "H4/3": 41, "H4/4": 59, "H4/6": -283, "H4/7": 326, "H4/8": 298,
                   "H5/0": 2.237, "H5/2": 1.248, "H5/4": 0.84, "H6/2": 213, "H6/4": 160, "H6/6": 106.5,
                   "Hcf/2,0": 216, "Hcf/4,0": 1225, "Hcf/6,0": 1506, "Hcf/6,6": 770},  # p. 57-58
        "terms": ["4I_9/2_-9/2", "4I_9/2_+5/2", "4I_9/2_+3/2", "4I_9/2_+1/2", "4I_9/2_-7/2", "4I_11/2_-11/2",
                  "4I_11/2_+5/2", "4I_11/2_+3/2", "4I_11/2_+1/2", "4I_11/2_-7/2", "4I_11/2_-9/2", "4I_13/2_+13/2",
                  "4I_13/2_+5/2", "4I_13/2_+3/2", "4I_13/2_+1/2", "4I_13/2_-7/2", "4I_13/2_-9/2", "4I_13/2_-11/2",
                  "4I_15/2_+15/2", "4I_15/2_-7/2", "4I_15/2_-9/2", "4I_15/2_+5/2", "4I_15/2_+1/2", "4I_15/2_+3/2",
                  "4I_15/2_-11/2", "4I_15/2_+13/2", "4F_3/2_+1/2", "4F_3/2_+3/2", "2H2_9/2_-7/2", "4F_5/2_+3/2",
                  "4F_5/2_+1/2", "4F_5/2_+5/2", "4F_5/2_+1/2", "2H2_9/2_+3/2", "2H2_9/2_-9/2", "2H2_9/2_+5/2",
                  "4F_7/2_+3/2", "4F_7/2_-7/2", "4F_7/2_+1/2", "4S_3/2_+1/2", "4S_3/2_+3/2", "4F_7/2_+5/2",
                  "4F_9/2_+1/2", "4F_9/2_+3/2", "4F_9/2_+3/2", "4F_9/2_+5/2", "4F_9/2_-7/2", "2H2_11/2_+5/2",
                  "2H2_11/2_-11/2", "2H2_11/2_+3/2", "2H2_11/2_-7/2", "2H2_11/2_+1/2", "2H2_11/2_-9/2", "4G_5/2_+3/2",
                  "4G_5/2_+1/2", "4G_5/2_+5/2", "4G_7/2_+5/2", "2G1_7/2_+3/2", "4G_5/2_+5/2", "4G_5/2_+1/2",
                  "4G_7/2_+5/2", "4G_7/2_+1/2", "4G_7/2_+3/2", "4G_7/2_-7/2", "2K_13/2_+13/2", "4G_9/2_+5/2",
                  "2K_13/2_+1/2", "4G_9/2_-7/2", "4G_9/2_-9/2", "2K_13/2_-11/2", "4G_9/2_+3/2", "2K_13/2_+3/2",
                  "4G_9/2_+1/2", "2K_13/2_+5/2", "2K_13/2_-9/2", "2K_13/2_-7/2", "2G1_9/2_-7/2", "2G1_9/2_+3/2",
                  "2G1_9/2_-7/2", "2G1_9/2_+1/2", "2G1_9/2_-9/2", "2D1_3/2_+3/2", "2D1_3/2_+1/2", "4G_11/2_-7/2",
                  "2K_15/2_+15/2", "4G_11/2_+5/2", "4G_11/2_-9/2", "4G_11/2_-11/2", "2K_15/2_+13/2", "4G_11/2_-9/2",
                  "2K_15/2_+1/2", "2K_15/2_+5/2", "2K_15/2_+3/2", "2K_15/2_-11/2", "2K_15/2_-9/2", "2K_15/2_-7/2",
                  "4G_11/2_+1/2", "2P_1/2_+1/2", "2D1_5/2_+3/2", "2D1_5/2_+1/2", "2D1_5/2_+5/2", "2P_3/2_+1/2",
                  "2P_3/2_+3/2", "4D_3/2_+1/2", "4D_3/2_+3/2", "4D_5/2_+5/2", "4D_5/2_+1/2", "2I_11/2_-11/2",
                  "4D_5/2_+3/2", "4D_1/2_+1/2", "2I_11/2_-9/2", "2I_11/2_-7/2", "2I_11/2_+5/2", "2I_11/2_+3/2",
                  "2I_11/2_+1/2", "2L_15/2_+13/2", "2L_15/2_+15/2", "2L_15/2_+1/2", "2L_15/2_+3/2", "2L_15/2_-11/2",
                  "2L_15/2_+5/2", "4D_7/2_+5/2", "2L_15/2_-9/2", "4D_7/2_-7/2", "2L_15/2_-7/2", "4D_7/2_+3/2",
                  "4D_7/2_+1/2", "2I_13/2_-11/2", "2I_13/2_-9/2", "2I_13/2_-7/2", "2I_13/2_+13/2", "2I_13/2_+5/2",
                  "2I_13/2_+3/2", "2I_13/2_+1/2", "2L_17/2_+15/2", "2L_17/2_+17/2", "2L_17/2_+1/2", "2L_17/2_+3/2",
                  "2L_17/2_+13/2", "2L_17/2_+5/2", "2L_17/2_-11/2", "2L_17/2_-9/2", "2L_17/2_-7/2", "2H1_9/2_-7/2",
                  "2H1_9/2_+1/2", "2H1_9/2_-9/2", "2H1_9/2_+5/2", "2H1_9/2_+3/2", "2D2_3/2_+3/2", "2D2_3/2_+1/2",
                  "2H1_11/2_-9/2", "2H1_11/2_+1/2", "2D2_5/2_+5/2", "2H1_11/2_-7/2", "2H1_11/2_+1/2", "2H1_11/2_+3/2",
                  "2H1_11/2_-11/2", "2D2_5/2_+3/2", "2H1_11/2_+5/2", "2F2_5/2_+5/2", "2F2_5/2_+1/2", "2F2_5/2_+3/2",
                  "2F2_7/2_-7/2", "2F2_7/2_+3/2", "2F2_7/2_+1/2", "2F2_7/2_+5/2", "2G2_9/2_+5/2", "2G2_9/2_-9/2",
                  "2G2_9/2_+3/2", "2G2_9/2_-7/2", "2G2_9/2_+1/2", "2G2_7/2_-7/2", "2G2_7/2_+3/2", "2G2_7/2_+5/2",
                  "2G2_7/2_+1/2", "2F1_7/2_+5/2", "2F1_7/2_-7/2", "2F1_7/2_+3/2", "2F1_7/2_+1/2", "2F1_5/2_+5/2",
                  "2F1_5/2_+3/2", "2F1_5/2_+1/2"],
        "kcalc": [3, 38, 142, 294, 500,
                  1963, 2042, 2075, 2098, 2203, 2227,
                  3901, 3983, 4043, 4102, 4126, 4205, 4275,
                  5820, 5838, 5997, 6171, 6187, 6293, 6420, 6545,
                  11596, 11626,
                  12585, 12589, 12630, 12678, 12704, 12763, 12854, 12873,
                  13514, 13583, 13673, 13693, 13695, 13711,
                  14847, 14861, 14886, 14924, 14957,
                  16028, 16046, 16059, 16060, 16095, 16140,
                  17308, 17311, 17360, 17491, 17505, 17564, 17611,
                  19139, 19245, 19271, 19324,
                  19567, 19632, 19645, 19687, 19693, 19737, 19738, 19790, 19845, 19917, 19927, 19971,
                  21150, 21183, 21199, 21235, 21267,
                  21339, 21352,
                  21535, 21620, 21621, 21731, 21774, 21779, 21790, 21824, 21827, 21857, 21901, 21931, 21946,
                  21995,
                  23455,
                  23996, 23999, 24057,
                  26394, 26416,
                  28361, 28369,
                  28495, 28528, 28634, 28686, 28938, 29459, 29475, 29565, 29659, 29767,
                  30271, 30346, 30411, 30451, 30533, 30534, 30602, 30615, 30646, 30701, 30712, 30792,
                  30850, 30895, 30955, 31002, 31041, 31070, 31079, 31767, 31836, 31926, 31968, 31990, 32013,
                  32048,
                  32093, 32126,
                  33035, 33137, 33168, 33226, 33258,
                  33612, 33631,
                  34274, 34374, 34445, 34519, 34551, 34573, 34686, 34709, 34818,
                  38723, 38778, 38815,
                  40113, 40126, 40187, 40254,
                  47871, 47888, 47964, 48006, 48055,
                  48861, 48869, 48979, 49065,
                  66548, 66705, 66793, 66859,
                  67857, 67858, 68075],  # p. 63
    },
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


def latex_level(intermediate):
    """ Return a LaTeX string representation of the main component of the state in intermediate coupling. """

    state = intermediate.states[0]
    num = f"({state['num']})" if state['num'] else ""
    M = "$%s$" % state['Jz'].replace("-", "").replace("+", "")
    level = r"$\prescript{%s}{}{\mathrm{%s}}%s_{%s}$" % (state["S2"], state["L2"], num, state["J2"])
    return level, M


class Table:
    def __init__(self):
        self.rows = []
        self.rows.append([r"\multicolumn{2}{|c}{Level}",
                          r"\multicolumn{1}{c}{$k^\mathrm{ref}$}",
                          r"\multicolumn{1}{c}{$k^\mathrm{calc}$}",
                          r"\multicolumn{1}{c|}{$\Delta_\mathrm{abs}$}"])
        self.rows.append([r"\multicolumn{1}{|l}{$\prescript{2S+1}{}{L}_{J}$}",
                          r"\multicolumn{1}{r}{$M$}",
                          r"\multicolumn{1}{c}{\unit{cm^{-1}}}",
                          r"\multicolumn{1}{c}{\unit{cm^{-1}}}",
                          r"\multicolumn{1}{c|}{\unit{cm^{-1}}}"])
        self.kdiff = []
        self.indices = None

    def add_row(self, state, k_ref, mark):
        self.kdiff.append(k_ref - state.energy)

        kref = f"${k_ref:.0f}$"
        kcalc = f"${state.energy:.1f}$"
        kdiff = f"${k_ref - state.energy:.1f}$"
        level, M = latex_level(state)

        row = [level, M, kref, kcalc, kdiff]
        if mark:
            i, letter = mark
            assert row[i][-1] == "$"
            row[i] = row[i][:-1] + f"^{letter}$"

        self.rows.append(row)

    def pop(self, index):
        if index < 0:
            index += len(self.rows)
        assert index not in self.indices
        self.indices.append(index)
        return self.rows[index]

    def get_code(self, col_blocks):
        krmse = math.sqrt(sum(diff ** 2 for diff in self.kdiff) / len(self.kdiff))

        num_cols = 5
        num_rows = math.ceil((len(self.rows) + 1) / col_blocks)
        for i in range(num_rows * col_blocks - len(self.rows) - 1):
            self.rows.append(num_cols * [""])
        self.rows.append([r"\multicolumn{1}{l}{RMS:}"] + (num_cols - 2) * [""] + [f"{krmse:.1f}"])

        self.indices = []
        rows = []
        for i in range(num_rows):
            line = []
            for j in range(col_blocks):
                line += self.pop(i + j * num_rows)
            rows.append(line)
        assert set(self.indices) == set(range(len(self.rows))), set(range(len(self.rows))) - set(self.indices)
        self.indices = None

        lines = []
        lines.append(r"\begin{tabular}{%s}" % ("|" + col_blocks * "lrrrr|"))
        lines.append(r"  \hline")
        for i, row in enumerate(rows):
            if i == 2:
                lines.append(r"  \cline{%d-%d}" % (1, num_cols))
            if i == len(rows) - 1:
                lines.append(r"  \cline{%d-%d}" % (num_cols * (col_blocks - 1) + 1, num_cols * col_blocks))
            lines.append("  " + " & ".join(row) + r" \\[-0.09ex]")
        lines.append(r"  \hline")
        lines.append(r"\end{tabular}")
        return "\n".join(lines)


if __name__ == "__main__":
    init_logger(level=logging.DEBUG)

    # Load literature data
    key = "Nd3+:LaF3"
    data = DATA[key]
    name = data["name"]
    config = ion2config(name)
    radial = data["radial"]
    coupling = Coupling.SLJM
    terms = data["terms"]
    k_calc = data["kcalc"]
    assert all(a < b for a, b in zip(k_calc[:-1], k_calc[1:]))

    ### BEGIN HOT FIX ##################################
    err_term = "2I_11/2_-11/2"
    k_old = 28634
    k_new = 29634
    i = k_calc.index(k_old)
    assert terms[i] == err_term
    j = bisect(k_calc, k_new) - 1
    terms.insert(j, terms.pop(i))
    k_calc.insert(j, k_calc.pop(i))
    assert k_calc[j] == k_old
    k_calc[j] = k_new
    assert all(a < b for a, b in zip(k_calc[:-1], k_calc[1:]))
    ### END HOT FIX ####################################

    # Calculate energy levels and intermediate states
    ion = Levels(config, coupling, radial)

    # Kramer's degeneracy test
    assert len(ion) % 2 == 0
    diff = [abs(ion.energies[2 * i] - ion.energies[2 * i + 1]) for i in range(len(ion) // 2)]
    assert max(diff) < 1e-10
    logger.info(f"Number of Kramer-doublets: {len(k_calc)}")

    # Adjust base level
    energies = [ion.energies[2 * i] for i in range(len(ion) // 2)]
    mean = sum(a - b for a, b in zip(k_calc, energies)) / len(k_calc)
    logger.info(f"Mean energy deviation: {mean:.1f}")
    ion.set_base(ion.energies[0] + mean)
    logger.info(f"Adjusted parameter 'base': {radial['base']} -> {ion.radial_integrals['base']:.2f}")

    # Pick one state from each Kramer's-doublet
    energies = [ion.energies[2 * i] for i in range(len(ion) // 2)]
    states = [ion.states[2 * i] for i in range(len(ion) // 2)]

    logger.info("Generating table...")
    diff_max = -1
    table = Table()
    for i in range(len(states)):
        state = states[i]
        name = state.short()
        kcalc = energies[i]
        term = terms[i]
        kref = k_calc[i]
        kdiff = kref - kcalc
        if abs(kdiff) > diff_max:
            diff_max = abs(kdiff)

        mark = None
        match = name.replace("-", "").replace("+", "") == term.replace("-", "").replace("+", "")
        if not match:
            mark = (1, "a")
        if term == err_term:
            mark = (2, "b")
        table.add_row(state, kref, mark)
    logger.info(f"Maximum energy difference: {diff_max:.1f}")

    # Result summary
    mean = sum(a - b for a, b in zip(k_calc, energies)) / len(k_calc)
    logger.info(f"Mean energy deviation: {mean:.1f}")
    sigma = math.sqrt(sum((a - b) ** 2 for a, b in zip(k_calc, energies)) / len(k_calc))
    logger.info(f"RMS energy: {sigma:.1f}")

    # Store table
    tex_file = "example_2.tex"
    with open(tex_file, "w") as fp:
        fp.write(table.get_code(3))
    logger.info(f"Generated LaTeX table file '{tex_file}'")
