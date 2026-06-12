##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This script generates the comparison table IX for the first application
# example in the AMELI paper.
#
##########################################################################

import logging
import math

from yalip import Coupling, Levels, ion2config, Sellmeier

logger = logging.getLogger("example_1")

# Literature data
DATA = {
    "Er3+:LaCl3": {
        "ref": "[Hehlen2013]",
        "name": "Er3+",
        "material": Sellmeier(B1=0.62646, C1=math.sqrt(6.1295e-2),
                              B2=1.5212, C2=math.sqrt(6.1087e-3),
                              B3=0.2465, C3=math.sqrt(1.4086e-2)),
        "term": ["4I_15/2", "4I_13/2", "4I_11/2", "4I_9/2", "4F_9/2", "4S_3/2", "2H2_11/2", "4F_7/2", "4F_5/2",
                 "4F_3/2", "2G1_9/2", "4G_11/2"],
        "radial": {"F_2": 433.22, "F_4": 66.887, "F_6": 7.2952, "H2": 2385.9},
        "kmeas": [0, 6520, 10165, 12375, 15172, 18383, 19068, 20414, 22099, 22510, 24467, 26290],
        "kcalc": [0, 6527, 10145, 12296, 15139, 18311, 19185, 20304, 21949, 22303, 24430, 26415],
        "omega": {"JO/2": 5.449, "JO/4": 2.077, "JO/6": 0.6873},
        "fmeas": [0, 1.738e-6, 0.5010e-6, 0.4858e-6, 3.043e-6, 0.3794e-6, 10.881e-6, 2.467e-6, 0.3939e-6, 0.2130e-6,
                  0.5976e-6, 22.682e-6],
        "fed": [0, 1.039e-6, 0.5270e-6, 0.5733e-6, 2.641e-6, 0.3502e-6, 11.56e-6, 1.897e-6, 0.4348e-6, 0.2505e-6,
                0.6142e-6, 21.59e-6],
        "fmd": [0, 0.5757e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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


def latex_state(intermediate, min_weight=0.0, max_num=None):
    """ Return a LaTeX string representation of the state. """

    result = []
    if max_num is None:
        weights = intermediate.weights[:3]
        states = intermediate.states[:3]
    else:
        assert isinstance(max_num, int)
        assert max_num > 1
        weights = intermediate.weights[:max_num]
        states = intermediate.states[:max_num]
    sum_weight = 0.0
    for weight, state in zip(weights, states):
        if weight >= min_weight:
            num = f"({state['num']})" if state['num'] else ""
            result.append(
                f"{weight:.2f} \\prescript{{{state['S2']}}}{{}}{{\\mathrm{{{state['L2']}}}}}{num}_{{{state['J2']}}}")
            sum_weight += weight
    result = " + ".join(result)
    if 1 - sum_weight >= min_weight:
        result += " \\ldots"
    return "$" + result + "$"


def latex_level(intermediate):
    """ Return a LaTeX string representation of the main component of the state in intermediate coupling. """

    state = intermediate.states[0]
    num = f"({state['num']})" if state['num'] else ""
    return f"$\\prescript{{{state['S2']}}}{{}}{{\\mathrm{{{state['L2']}}}}}{num}_{{{state['J2']}}}$"


class Table:
    def __init__(self, min_weight=0.0):
        self.min_weight = abs(float(min_weight))

        head = ["\\multicolumn{1}{|c}{$k^\\mathrm{ref}$}",
                "\\multicolumn{1}{c}{$k^\\mathrm{calc}$}",
                "\\multicolumn{1}{c|}{$\\Delta_\\mathrm{abs}$}",
                "\\multicolumn{1}{|c}{$f_{ed}^\\mathrm{ref}$}",
                "\\multicolumn{1}{c}{$f_{md}^\\mathrm{ref}$}",
                "\\multicolumn{1}{c}{$f_{ed}^\\mathrm{calc}$}",
                "\\multicolumn{1}{c}{$f_{md}^\\mathrm{calc}$}",
                "\\multicolumn{1}{c|}{$\\Delta_\\mathrm{rel}$}",
                "\\multicolumn{1}{|c|}{Intermediate State}"]
        unit = ["\\multicolumn{1}{|c}{\\unit{cm^{-1}}}",
                "\\multicolumn{1}{c}{\\unit{cm^{-1}}}",
                "\\multicolumn{1}{c|}{\\unit{cm^{-1}}}",
                "\\multicolumn{1}{|c}{$10^{-8}$}",
                "\\multicolumn{1}{c}{$10^{-8}$}",
                "\\multicolumn{1}{c}{$10^{-8}$}",
                "\\multicolumn{1}{c}{$10^{-8}$}",
                "\\multicolumn{1}{c|}{\\%}",
                ""]
        self.head = [
            "\\begin{tabular}{|rrr|rrrrr|l|}",
            "  \\hline",
            f"  {" &".join(head)} \\\\",
            f"  {" &".join(unit)} \\\\",
            "  \\hline"]

        self.rows = []
        self.kdiff = []
        self.fdiff = []

    def add_row(self, state, k_ref, fed_calc, fmd_calc, fed_ref, fmd_ref):
        fmin = 1e-15

        if k_ref is None:
            kref = kdiff = ""
            fedref = fmdref = fdiff = ""
        else:
            kref = f"{k_ref:.0f}"
            kdiff = f"{k_ref - state.energy:.1f}"
            self.kdiff.append(k_ref - state.energy)

            fedref = "0" if abs(fed_ref) < fmin else f"{fed_ref * 1e8:.2f}"
            fmdref = "0" if abs(fmd_ref) < fmin else f"{fmd_ref * 1e8:.2f}"
            if fed_ref + fmd_ref:
                fdiff = (fed_ref + fmd_ref) / (fed_calc + fmd_calc) - 1
                self.fdiff.append(fdiff)
                fdiff = f"{100 * fdiff:.1f}"
            else:
                fdiff = ""

        kcalc = f"{state.energy:.1f}"
        fedcalc = "0" if abs(fed_calc) < fmin else f"{fed_calc * 1e8:.2f}"
        fmdcalc = "0" if abs(fmd_calc) < fmin else f"{fmd_calc * 1e8:.2f}"

        # level = latex_level(state)
        level = latex_state(state, min_weight=0.01, max_num=3)
        self.rows.append(
            f"  {kref} & {kcalc} & {kdiff} & {fedref} & {fmdref} & {fedcalc} & {fmdcalc} & {fdiff} &\n  {level} \\\\")

    def get_code(self):
        krmse = math.sqrt(sum(diff ** 2 for diff in self.kdiff) / len(self.kdiff))
        frmse = math.sqrt(sum(diff ** 2 for diff in self.fdiff) / len(self.fdiff))
        lines = self.head + self.rows + [
            "  \\hline",
            f"  \\multicolumn{{1}}{{|l}}{{RMS:}} & & {krmse:.1f} & & & & & {100 * frmse:.1f} & \\\\",
            "  \\hline",
            "\\end{tabular}"]
        return "\n".join(lines)


if __name__ == "__main__":
    init_logger(level=logging.DEBUG)

    # Load literature data
    key = "Er3+:LaCl3"
    data = DATA[key]
    name = data["name"]
    config = ion2config(name)
    radial = data["radial"]
    jo = data["omega"]
    material = data["material"]
    coupling = Coupling.SLJ
    k_calc = data["kcalc"]
    f_ed = data["fed"]
    f_md = data["fmd"]

    # Calculate energy levels and intermediate states
    ion = Levels(config, coupling, radial, jo, material)

    # Adjust base level
    mean = sum(a - b for a, b in zip(k_calc, ion.energies[:len(k_calc)])) / len(k_calc)
    logger.info(f"Mean energy deviation: {mean:.1f}")
    ion.set_base(ion.energies[0] + mean)
    logger.info(f"Adjusted parameter 'base': {radial['base']} -> {ion.radial_integrals['base']:.2f}")

    logger.info("Generating table...")
    f = ion.oscillator_strengths()
    diff_rel = []
    table = Table()
    for i in range(len(ion)):
        if i < len(k_calc):
            kref = k_calc[i]
            fed = f_ed[i]
            fmd = f_md[i]
            if i > 0:
                diff_rel.append((fed + fmd) / (f.ed[i, 0] + f.md[i, 0]) - 1)
        else:
            kref = fed = fmd = None
        table.add_row(ion[i], kref, f.ed[i, 0], f.md[i, 0], fed, fmd)

    # Result summary
    mean = sum(a - b for a, b in zip(k_calc, ion.energies[:len(k_calc)])) / len(k_calc)
    logger.info(f"Mean energy deviation: {mean:.1f}")
    sigma_k = math.sqrt(sum((a - b) ** 2 for a, b in zip(k_calc, ion.energies[:len(k_calc)])) / len(k_calc))
    logger.info(f"RMS energy: {sigma_k:.1f}")
    sigma_f = math.sqrt(sum(d ** 2 for d in diff_rel) / len(diff_rel))
    logger.info(f"RMS oscillator strength: {sigma_f * 100:.1f} %")

    # Store table
    tex_file = "example_1.tex"
    with open(tex_file, "w") as fp:
        fp.write(table.get_code())
    logger.info(f"Generated LaTeX table file '{tex_file}'")
