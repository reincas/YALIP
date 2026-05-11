##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################

import math

from lanthanide import product_states, init_single, get_matrix, ORBITAL, SPIN, \
    SymmetryList, StateListProduct, init_states

SYMS = {
    1: {
        "S2": {'2': 14},
        "GR/7": {'(100)': 14},
        "GG/2": {'[10]': 14},
        "L2": {'F': 14},
        "tau": {'*': 14},
        "J2": {'5/2': 6, '7/2': 8},
        "num": {'0': 14},
        "Jz": {'-7/2': 1, '-5/2': 2, '-3/2': 2, '-1/2': 2, '+1/2': 2, '+3/2': 2, '+5/2': 2, '+7/2': 1},
    },
    2: {
        "S2": {'1': 28, '3': 63},
        "GR/7": {'(000)': 1, '(110)': 63, '(200)': 27},
        "GG/2": {'[00]': 1, '[10]': 21, '[11]': 42, '[20]': 27},
        "L2": {'S': 1, 'P': 9, 'D': 5, 'F': 21, 'G': 9, 'H': 33, 'I': 13},
        "tau": {'*': 91},
        "J2": {'0': 2, '1': 3, '2': 15, '3': 7, '4': 27, '5': 11, '6': 26},
        "num": {'0': 91},
        "Jz": {'-6': 2, '-5': 3, '-4': 6, '-3': 7, '-2': 10, '-1': 11, '+0': 13, '+1': 11, '+2': 10, '+3': 7, '+4': 6,
               '+5': 3, '+6': 2},
    },
    12: {
        "S2": {'1': 28, '3': 63},
        "GR/7": {'(000)': 1, '(110)': 63, '(200)': 27},
        "GG/2": {'[00]': 1, '[10]': 21, '[11]': 42, '[20]': 27},
        "L2": {'S': 1, 'P': 9, 'D': 5, 'F': 21, 'G': 9, 'H': 33, 'I': 13},
        "tau": {'*': 91},
        "J2": {'0': 2, '1': 3, '2': 15, '3': 7, '4': 27, '5': 11, '6': 26},
        "num": {'0': 91},
        "Jz": {'-6': 2, '-5': 3, '-4': 6, '-3': 7, '-2': 10, '-1': 11, '+0': 13, '+1': 11, '+2': 10, '+3': 7, '+4': 6,
               '+5': 3, '+6': 2},
    },
    13: {
        "S2": {'2': 14},
        "GR/7": {'(100)': 14},
        "GG/2": {'[10]': 14},
        "L2": {'F': 14},
        "tau": {'*': 14},
        "J2": {'5/2': 6, '7/2': 8},
        "num": {'0': 14},
        "Jz": {'-7/2': 1, '-5/2': 2, '-3/2': 2, '-1/2': 2, '+1/2': 2, '+3/2': 2, '+5/2': 2, '+7/2': 1},
    },
}


class DummyLanthanide:
    def __init__(self, num: int):
        self.num = num
        self.l = ORBITAL
        self.s = SPIN
        self.product = product_states(num)
        self.single = init_single(None, None, self.product)

    def matrix(self, name):
        return get_matrix(self, name)


def run_state(num, num_slj):
    ion = DummyLanthanide(num)

    states = StateListProduct(ion.product)
    assert len(states) == math.comb(14, num)

    states = states.to_SLJM(ion)
    assert len(states) == math.comb(14, num)
    assert type(states[1]).__name__ == "StateSLJM"
    for i, name in enumerate(states.chain):
        result = SymmetryList(states.values[:, i], name).count()
        #print(f'"{name}": {result},')
        assert result == SYMS[num][name]

    states = states.to_SLJ()
    assert len(states) == num_slj
    assert type(states[1]).__name__ == "StateSLJ"

    states = init_states(None, None, ion)
    assert set(states.keys()) == {'Product', 'SLJM', 'SLJ'}


def test_state_Ce():
    run_state(1, 2)


def test_state_Pr():
    run_state(2, 13)


def test_state_Tm():
    run_state(12, 13)


def test_state_Yb():
    run_state(13, 2)
