##########################################################################
# Copyright (c) 2026 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides a collection of convenient parameter sets from the literature.
#
#   [1]  W. T. Carnall, G. L. Goodman, K. Rajnak, and R. S. Rana
#        "A systematic analysis of the spectra of the lanthanides doped
#        into single crystal LaF3"
#        J. Chem. Phys. 90 (1989), no. 7, pp. 3443-3457
#
#   [2]  W. T. Carnall, H. Crosswhite, H. M. Crosswhite
#        "Energy level structure and transition probabilities of the
#        trivalent lanthanides in LaF3"
#        ANL-78-XX-95, Argonne National Laboratory Report, 1978
#        (available as microfiche)
#
#   [3]  Reinhard Caspary
#        "Applied Rare-Earth Spectroscopy for Fiber Laser Optimization"
#        Dissertation, Technical university Braunschweig, published with Shaker, Aachen, 2002
#
#   [4]  Dimitar N. Petrov, B.M. Angelov
#        "Spin – Orbit interaction in Yb3+ – Ground level and nephelauxetic
#        effect in crystals"
#        Chem. Phys. 525 (2019), 110416
#        https://doi.org/10.1016/j.chemphys.2019.110416
#
#   [5]  L. Wetenkamp
#        Charakterisierung von laseraktiv dotierten Schwermetallfluorid-Gläsern und Faserlasern
#        Dissertation, Technical university Braunschweig, 1991
#
#   [6]  I. H. Malitson
#        "Interspecimen comparison of the refractive index of fused silica"
#        J. Opt. Soc. Am. 55, 1205-1208 (1965)
#
##########################################################################


import logging

from .spectrum import Cauchy, Sellmeier

logger = logging.getLogger("yalip.lanthanide")

# Symbols of the 15 Lanthanides. The configurations of the triply ionized atoms are 4f0 - 4f14
LANTHANIDES = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]

# Radial integrals
RADIAL = {
    "Ce3+":  # from [1]
        {"base": 217, "H2": 647.3},
    "Pr3+":  # from [1]
        {"base": 205, "H1/2": 68878, "H1/4": 50347, "H1/6": 32901, "H2": 751.7,
         "H3/0": 16.23, "H3/1": -566.6, "H3/2": 1371, "H5fix": 2.08, "H6fix": 88.6},
    "Pr3+/alt":  # from [2]
        {"base": 191, "H1/2": 69305, "H1/4": 50675, "H1/6": 32813, "H2": 750.8,
         "H3/0": 21, "H3/1": -842, "H3/2": 1625, "H5fix": 1.99, "H6fix": 200},
    "Pr3+/ZBLAN":  # from [3]
        {"base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
         "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67},
    "Nd3+":  # from [2]
        {"base": 235, "H1/2": 73036, "H1/4": 52624, "H1/6": 35793, "H2": 884.9,
         "H3/0": 21.28, "H3/1": -583, "H3/2": 1443,
         "H4/2": 306, "H4/3": 41, "H4/4": 59, "H4/6": -283, "H4/7": 326, "H4/8": 298, "H5fix": 2.237, "H6fix": 213},
    "Pm3+":  # from [2]
        {"base": 120, "H1/2": 77000, "H1/4": 55000, "H1/6": 37500, "H2": 1022,
         "H3/0": 21.0, "H3/1": -560, "H3/2": 1400,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 2.49, "H6fix": 440},
    "Sm3+":  # from [2]
        {"base": 101, "H1/2": 79915, "H1/4": 57256, "H1/6": 40424, "H2": 1177.2,
         "H3/0": 20.07, "H3/1": -563, "H3/2": 1436,
         "H4/2": 288, "H4/3": 36, "H4/4": 56, "H4/6": -283, "H4/7": 333, "H4/8": 342, "H5fix": 2.76, "H6fix": 344},
    "Eu3+":  # from [2]
        {"base": 0, "H1/2": 84000, "H1/4": 60000, "H1/6": 42500, "H2": 1327,
         "H3/0": 20, "H3/1": -570, "H3/2": 1450,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.03, "H6fix": 300},
    "Gd3+":  # from [2]
        {"base": 2, "H1/2": 85587, "H1/4": 61361, "H1/6": 45055, "H2": 1503.5,
         "H3/0": 20, "H3/1": -590, "H3/2": 1450,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.32, "H6fix": 611},
    "Tb3+":  # from [2]
        {"base": 124, "H1/2": 91220, "H1/4": 65798, "H1/6": 43661, "H2": 1702.2,
         "H3/0": 19.81, "H3/1": -600, "H3/2": 1400,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 3.61, "H6fix": 583},
    "Dy3+":  # from [2]
        {"base": 175, "H1/2": 94877, "H1/4": 67470, "H1/6": 45745, "H2": 1912,
         "H3/0": 17.64, "H3/1": -608, "H3/2": 1498,
         "H4/2": 423, "H4/3": 50, "H4/4": 117, "H4/6": -334, "H4/7": 432, "H4/8": 353, "H5fix": 3.92, "H6fix": 771},
    "Ho3+":  # from [2]
        {"base": 9, "H1/2": 97025, "H1/4": 68885, "H1/6": 47744, "H2": 2144.2,
         "H3/0": 18.98, "H3/1": -579, "H3/2": 1570,
         "H4/2": 330, "H4/3": 41.5, "H4/4": 62, "H4/6": -295, "H4/7": 360, "H4/8": 310, "H5fix": 4.25, "H6fix": 843},
    "Er3+":  # from [1]
        {"base": 219, "H1/2": 97483, "H1/4": 67904, "H1/6": 54010, "H2": 2376,
         "H3/0": 17.79, "H3/1": -582.1, "H3/2": 1800,
         "H4/2": 400, "H4/3": 43, "H4/4": 73, "H4/6": -271, "H4/7": 308, "H4/8": 299, "H5fix": 3.86, "H6fix": 594},
    "Er3+/alt":  # from [2]
        {"base": 217, "H1/2": 100274, "H1/4": 70555, "H1/6": 49900, "H2": 2381,
         "H3/0": 17.88, "H3/1": -599, "H3/2": 1719,
         "H4/2": 441, "H4/3": 42, "H4/4": 64, "H4/6": -314, "H4/7": 387, "H4/8": 363, "H5fix": 4.58, "H6fix": 852},
    "Er3+/ZBLAN":  # from [3]
        {"base": 50.99, "H1/2": 97088.92, "H1/4": 68587.69, "H1/6": 55006.43, "H2": 2369.69,
         "H3/0": 17.81, "H3/1": -559.04, "H3/2": 1603.10,
         "H4/2": 471.61, "H4/3": 20.39, "H4/4": 18.81, "H4/6": -398.11, "H4/7": 199.03, "H4/8": 449.82,
         "H5fix": 4.66, "H6fix": 475.64},
    "Tm3+":  # from [1]
        {"base": 250, "H1/2": 100134, "H1/4": 69613, "H1/6": 55975, "H2": 2636,
         "H3/0": 17.26, "H3/1": -624.5, "H3/2": 1820, "H5fix": 3.81, "H6fix": 695},
    "Tm3+/alt":  # from [2]
        {"base": 175, "H1/2": 102459, "H1/4": 72424, "H1/6": 51380, "H2": 2640,
         "H3/0": 17, "H3/1": -737, "H3/2": 1700, "H5fix": 4.93, "H6fix": 729.6},
    "Tm3+/ZBLAN":  # from [3]
        {"base": 149.80, "H1/2": 102403.01, "H1/4": 73241.80, "H1/6": 50320.22, "H2": 2583.47,
         "H3/0": 17.43, "H3/1": -841.95, "H3/2": 1820, "H5fix": -2.35, "H6fix": 12.78},
    "Yb3+":  # guess based on [4]
        {"base": 0, "H2": 2900},
}

# Judd-Ofelt parameters from [3]
JUDD_OFELT = {
    "Pr3+/ZBLAN":  # from [3]
        {"JO/2": 1.981, "JO/4": 4.645, "JO/6": 6.972},
    "Er3+/ZBLAN":  # from [3]
        {"JO/2": 2.915, "JO/4": 1.464, "JO/6": 1.184},
    "Tm3+/ZBLAN":  # from [3]
        {"JO/2": 2.920, "JO/4": 1.856, "JO/6": 0.670},
}

# Spectral refractive index coefficients
MATERIAL = {
    # ZBLAN doped with 2.5 % PbF_2 from [3]
    "Pb:ZBLAN": Cauchy(1.35123e-5, 2.94780e-3, 1.49985, -1.30933e-3, -3.23335e-6),
    # ZBLAN from [5]
    "ZBLAN": Cauchy(1.35123e-5, 2.94780e-3, 1.48965, -1.30933e-3, -3.23335e-6),
    # Silica glass from [6]
    "SiO2": Sellmeier(0.6961663, 0.4079426, 0.8974794, 0.0684043, 0.1162414, 9.896161)
}


def config2ion(config):
    """ Convert configuration string into ion string. """

    assert isinstance(config, str)
    assert config[:1].lower() == "f"
    num = int(config[1:])
    assert num > 0 and num < len(LANTHANIDES)
    return f"{LANTHANIDES[num]}3+"


def ion2config(ion):
    """ Convert ion string into configuration string. """

    assert isinstance(ion, str)
    assert ion[2:] == "3+"
    num = LANTHANIDES.index(ion[:2])
    return f"f{num}"
