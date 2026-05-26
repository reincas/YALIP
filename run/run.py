##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import logging
import math

from yalip.matrix import normalize_radial
from yalip.fits import format_params
from yalip import Coupling, Levels, ion2config, Sellmeier, Fits

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


def double(input):
    return [element for element in input for _ in range(2)]


if __name__ == "__main__":
    init_logger(level=logging.DEBUG)

    DATA = {
        "Pr3+:LaF3/1": {
            "ref": "[Carnall1969/1970]",
            "name": "Pr3+",
            "radial": normalize_radial({"base": 225, "E^1": 4567.2, "E^2": 21.954, "E^3": 467.75,
                                        "H2": 721.90, "H3/0": 15.294, "H3/1": -669.02, "H3/2": 1413.7,
                                        "H5/0": 1.231, "H5/2": 0.690, "H5/4": 0.468, "P_2": 2.697, "P_4": 0.343,
                                        "P_6": 0.324}),
            "calc": [232, 2317, 4502, 5190, 6586, 7045, 9992, 17047, 20905, 21514, 21534, 22748, 46986],
            "stage": [":H1/2", ":H1/4", ":H1/6", ":H2", ":H3/0", ":H3/1", ":H3/2", "H5/0", "H5/2", "H5/4", ":H6/2",
                      ":H6/4", ":H6/6"],
            "fit": False,
        },
        "Tm3+:LaF3/0": {  # Error 1S0 74450 -> 74407
            "ref": "[Carnall1970]",
            "name": "Tm3+",
            "radial": {"base": 122, "E^1": 6737.5, "E^2": 33.643, "E^3": 681.22,
                       "H2": 2633.0, "H3/0": 13.124, "H3/1": -743.02, "H3/2": 1992.2},
            "calc": [119, 5835, 8320, 12701, 14586, 15150, 21366, 28098, 34866, 35604, 36562, 38315, 74407],
            "fit": False,
        },
        "Tm3+:LaF3/1": {  # Level 1S0 missing
            "ref": "[Carnall1970]",
            "name": "Tm3+",
            "radial": normalize_radial({"base": 106, "E^1": 6911.8, "E^2": 33.728, "E^3": 675.28,
                                        "H2": 2593.9, "H3/0": 9.475, "H3/1": -601.09, "H3/2": 1395.1,
                                        "H5/0": 5.002, "H5/2": 2.801, "H5/4": 1.901, "P_2": 3.915, "P_4": 0.090,
                                        "P_6": 0.070}),
            "calc": [100 + 6.7, 5858 - 9.2, 8336 - 1.8, 12711 - 4.4, 14559 + 3.4, 15173 + 2.1, 21352 + 4.2, 28061 + 0.8,
                     34866 - 0.25, 35604 - 0.16, 36559 + 0.08, 38344 + 1.4, 76505],
            # "calc": [106, 5858, 8336, 12701, 14559, 15173, 21352, 28061,
            #         34866, 35604, 36559, 38344, 76503],
            # "calc": [119, 5835, 8320, 12701, 14586, 15150, 21366, 28098, 34866, 35604, 36562, 38315, 76517],
            "stage": [":H1/2", ":H1/4", ":H1/6", ":H2", "H3/0", "H3/1", "H3/2", ":H5/0", ":H5/2", ":H5/4", ":H6fix"],
            "fit": False,
        },
        "Pr3+:LaF3/2": {
            "ref": "[Carnall1978]",
            "name": "Pr3+",
            "radial": {"base": 190, "H1/2": 69305, "H1/4": 50675, "H1/6": 32813,
                       "H2": 750.8, "H3/0": 21, "H3/1": -842, "H3/2": 1625,
                       "H5/0": 1.99, "H5/2": 1.11, "H5/4": 0.75, "H6fix": 200}, # p. 57
            "calc": [191, 2303, 4495, 5196, 6595, 7009, 10012, 17052, 20935, 21555, 21743, 22690, 46986],  # p. 60
            "fit": False,
        },
        "Nd3+:LaF3/2": {
            "ref": "[Carnall1978]",
            "name": "Nd3+",
            "radial": {"base": 235, "H1/2": 73036, "H1/4": 52624, "H1/6": 35793,
                       "H2": 884.9, "H3/0": 21.28, "H3/1": -583, "H3/2": 1443,
                       "H4/2": 306, "H4/3": 41, "H4/4": 59, "H4/6": -283, "H4/7": 326, "H4/8": 298,
                       "H5/0": 2.237, "H5/2": 1.248, "H5/4": 0.84, "H6fix": 213}, # p. 57
            "calc": [235, 2114, 4098, 6148, 11621, 12660, 12768, 13619, 13691, 14899, 16105, 17428, 17469, 19293,
                     19709, 19785, 21276, 21425, 21714, 21780, 23458, 24004, 26424],  # p. 68  Error: 21276 missing
            "stage": ["H1/2", "H1/4", "H1/6", "H2",
                      ":H3/0", ":H3/1", ":H3/2", ":H4/2", ":H4/3", ":H4/4", ":H4/6", ":H4/7", ":H4/8",
                      ":H5/0", ":H5/2", "H5/4", ":H6fix"],
            "fit": False,
        },
        "Nd3+:LaF3/2c": {
            "ref": "[Carnall1978]",
            "name": "Nd3+",
            "radial": {"base": 3, "H1/2": 73036, "H1/4": 52624, "H1/6": 35793,
                       "H2": 884.9, "H3/0": 21.28, "H3/1": -583, "H3/2": 1443,
                       "H4/2": 306, "H4/3": 41, "H4/4": 59, "H4/6": -283, "H4/7": 326, "H4/8": 298,
                       "H5/0": 2.237, "H5/2": 1.248, "H5/4": 0.84, "H6fix": 213,
                       "Hcf/2,0": 216 / 2, "Hcf/4,0": 1225 / 2, "Hcf/6,0": 1506 / 2, "Hcf/6,6": 770},  # p. 57-58
            "calc": double([3, 38, 142, 294, 500,
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
                     21150, 21187, 21199, 21235, 21267,
                     21339, 21352,
                     21535, 21620, 21621, 21731, 21774, 21779, 21790, 21824, 21827, 21857, 21901, 21931, 21946, 21995,
                     23455,
                     23996, 23999, 24057,
                     26394, 26416,
                     28361, 28369,
                     28495, 28528, 28686, 28938, 29459, 29475, 29565, 29634, 29659, 29767, # 28634 -> 29634
                     30271, 30346, 30411, 30451, 30533, 30534, 30602, 30615, 30646, 30701, 30712, 30792,
                     30850, 30895, 30955, 31002, 31041, 31070, 31079, 31767, 31836, 31926, 31968, 31990, 32013, 32048,
                     32093, 32126,
                     33035, 33137, 33168, 33226, 33258,
                     33612, 33631,
                     34274, 34374, 34445, 34519, 34551, 34573, 34686, 34709, 34818,
                     38723, 38778, 38815,
                     40113, 40126, 40187, 40254,
                     47871, 47888, 47964, 48006, 48055,
                     48861, 48869, 48979, 49065,
                     66548, 66705, 66793, 66859,
                     67857, 67858, 68075]),  # p. 63
            "stage": ["base", ":H1/2", ":H1/4", ":H1/6", ":H2",
                      ":H3/0", ":H3/1", ":H3/2", ":H4/2", ":H4/3", ":H4/4", ":H4/6", ":H4/7", ":H4/8",
                      ":H5/0", ":H5/2", ":H5/4", ":H6fix",
                      "Hcf/2,0", "Hcf/4,0", "Hcf/6,0", "Hcf/6,6"],
            "fit": False,
        },
        "Pm3+:LaF3/2": {
            "ref": "[Carnall1978]",
            "name": "Pm3+",
            "radial": {"base": 120, "H1/2": 77000, "H1/4": 55000, "H1/6": 37500,
                       "H2": 1010, "H3/0": 21, "H3/1": -560, "H3/2": 1400, # Error H2 1022 -> 1010
                       "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310,
                       "H5/0": 2.49, "H5/2": 1.39, "H5/4": 0.61, "H6fix": 440}, # p. 57  Error H5/4: 0.94 - > 0.61
            "calc": [120, 1612, 3239, 4951, 6714, 12638, 13080, 13933, 14486, 14887, 16223, 16939, 18053, 18075,
                     18255, 18565, 19862, 20307, 20554, 21935, 22475, 22807, 23140, 23772, 24216, 24702, 24840,
                     24907, 25811, 25895, 25907],  # p. 95
            "stage": [":H1/2", ":H1/4", ":H1/6", ":H2",
                      ":H3/0", ":H3/1", ":H3/2", ":H4/2", ":H4/3", ":H4/4", ":H4/6", ":H4/7", ":H4/8",
                      ":H5/0", ":H5/2", "H5/4", ":H6fix"],
            "fit": False,
        },
        "Tm3+:LaF3/2": {  # Error: "H2": 2640 -> 2587
            "ref": "[Carnall1978]",
            "name": "Tm3+",
            "radial": {"base": 170, "H1/2": 102459, "H1/4": 72424, "H1/6": 51380,
                       "H2": 2587, "H3/0": 17, "H3/1": -737, "H3/2": 1700, # Error H2: 2640 -> 2587
                       "H5/0": 4.93, "H5/2": 2.72, "H5/4": 1.37, "H6fix": 729.6}, # p. 57
            "calc": [175, 5818, 8391, 12721, 14597, 15181, 21314, 28001, 34975, 35579, 36615, 38268, 75300],  # p. 190
            "stage": [":H1/2", ":H1/4", ":H1/6", "H2", ":H3/0", ":H3/1", ":H3/2", ":H5/0", ":H5/2", ":H5/4", ":H6fix"],
            "fit": False,
        },
        "Nd3+:LaF3/3": {
            "ref": "[LizarazoFerro2026]",
            "name": "Nd3+",
            "radial": {"base": 3, "H1/2": 73040, "H1/4": 52790, "H1/6": 35770,
                       "H2": 885, "H3/0": 21.4, "H3/1": -590, "H3/2": 1430,
                       "H4/2": 292, "H4/3": 36, "H4/4": 60, "H4/6": -288, "H4/7": 339, "H4/8": 305,
                       "H5fix": 2.2, "H6fix": 210,
                       "Hcf/2,0": -250, "Hcf/4,0": 500, "Hcf/6,0": 640, "Hcf/2,2": -50, "Hcf/4,2": 500, "Hcf/6,2": -830,
                       "Hcf/4,4": 560, "Hcf/6,4": -400, "Hcf/6,6": -830},
            "calc": [235, 2114, 4098, 6148, 11621, 12660, 12768, 13619, 13691, 14899, 16105, 17428, 17469, 19293,
                     19709, 19785, 21425, 21714, 21780, 23458, 24004, 26424],  # p. 68
            "stage": [":H1/2", ":H1/4", ":H1/6", ":H2",
                      ":H3/0", ":H3/1", ":H3/2", "H4/2", "H4/3", "H4/4", "H4/6", "H4/7", "H4/8",
                      ":H5/0", ":H5/2", "H5/4", ":H6fix"],
            "fit": False,
        },
    }

    # jo = {"JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972}

    # LaF3 ordinary index [Laiho1983]
    material_o = Sellmeier(1 / 0.6505, 0.0, 0.0, math.sqrt(5.044e-3 / 0.6505), 0.0, 0.0)
    # LaF3 extraordinary index [Laiho1983]
    material_e = Sellmeier(1 / 0.6600, 0.0, 0.0, math.sqrt(5.093e-3 / 0.6600), 0.0, 0.0)
    ################## Refractive index in [Carnall1978]

    jo = None
    material = None

    data = DATA["Tm3+:LaF3/2"]
    name = data["name"]
    config = ion2config(name)
    radial = data["radial"]
    coupling = Coupling.SLJM if any(k.startswith("Hcf/") for k in radial) else Coupling.SLJ
    calc = data["calc"]
    assert all(a <= b for a, b in zip(calc[:-1], calc[1:]))

    if data["fit"]:
        lines = [[i, "", calc[i], 1.0, 0.0, 1.0] for i in range(1, len(calc))]
        opt = Fits(config, coupling, radial)
        ion = opt.run(lines, data["stage"])
    else:
        ion = Levels(config, coupling, radial, jo, material)

    mean = sum(a - b for a, b in zip(calc, ion.energies[:len(calc)])) / len(calc)
    print(f"Mean deviation: {mean:.1f}")
    #mean = radial["base"] - ion.energies[0]
    ion.set_base(ion.energies[0] + mean)

    print("List of states in intermediate coupling:")
    for i, state in enumerate(ion.str_levels(min_weight=0.05)):
        if i < len(calc):
            print(f"    {calc[i] - ion.energies[i]:6.1f} {calc[i]:6.0f} {state}")
        else:
            # pass
            print(f"    {"":6s} {"":6s} {state}")

    mean = sum(a - b for a, b in zip(calc, ion.energies[:len(calc)])) / len(calc)
    print(f"Mean deviation: {mean:.1f}")
    sigma = math.sqrt(sum((a - b) ** 2 for a, b in zip(calc, ion.energies[:len(calc)]))) / len(calc)
    print(f"Sigma: {sigma:.1f}")

    print(format_params(radial, 6))
    print(format_params(ion.radial_integrals, 6))
    # print(", ".join([f'"{k}"' for k in radial.keys()]))
    print(len(calc), len(ion.energies))

    # reduced = ion.dipole
    # R = np.column_stack((reduced.U2[1:, 0], reduced.U4[1:, 0], reduced.U6[1:, 0], reduced.LS[1:, 0]))
    # print("Squared reduced matrix elements (U2, U4, U6, LS):")
    # for line in R:
    #     line = "  ".join([f"{v:7.4f}" for v in line])
    #     print(f"    {line}")

    # f = ion.oscillator_strengths()
    # f = np.column_stack((f.ed[1:, 0], f.md[1:, 0])) * 1e8
    # print("GSA oscillator strengths (ed, md) in 1e-8:")
    # for line in f:
    #     line = "  ".join([f"{v:7.1f}" for v in line])
    #     print(f"    {line}")
    #
    # A = ion.radiative_rates()
    # i = len(ion) - 1
    # A = np.column_stack((A.ed[i - 1::-1, i], A.md[i - 1::-1, i]))
    # print("Radiative emission rates (ed, md) to ground state in 1/s:")
    # for line in A:
    #     line = "  ".join([f"{v:7.0f}" for v in line])
    #     print(f"    {line}")
    #
    # t = ion.life_times()[1:]
    # print("Radiative life times in ms:")
    # line = "  ".join([f"{v:.2f}" for v in t * 1000])
    # print(f"    {line}")
    #
    # beta = ion.branching_ratios()
    # i = len(ion) - 1
    # beta = beta[i - 1::-1, i]
    # print("Branching ratios:")
    # line = "  ".join([f"{v:.3f}" for v in beta])
    # print(f"    {line}")
