##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# Compare calculated LS states with the states listed in the book from
# Nielson and Koster.
#
##########################################################################

import pytest

from lanthanide import Lanthanide, Coupling
from data_states import SOURCES, STATES

def seniority(S2, W):
    l = 3
    num2 = W[1:-1].count("2")
    num1 = W[1:-1].count("1")
    num0 = W[1:-1].count("0")
    assert W == f"({num2 * "2"}{num1 * "1"}{num0 * "0"})"

    c = num2
    d = num1
    S2 = int(S2) - 1
    if d == S2:
        v = 2 * c + d
    elif 2 * l + 1 - 2 * c - d == S2:
        v = 2 * c + S2
    else:
        raise RuntimeError("Unknown seniority number!")
    return v


@pytest.mark.parametrize("num", range(1, 14))
def test_states(num):
    with Lanthanide(num) as ion:
        states = set()
        for state in ion.states(Coupling.SLJ):
            S2 = state["S2"].repr
            L2 = state["L2"].repr
            W = state["GR/7"].repr
            sen = seniority(S2, W)
            U = f"({state["GG/2"].repr[1:-1]})"
            LSnum = state["num"].repr
            LSnum = "" if LSnum == "0" else LSnum
            LStau = state["tau"].repr
            LStau = "" if LStau == "*" else LStau.upper()
            state_str = f"{S2}{L2}{LSnum} {sen} {W} {U}{LStau}"
            states.add(state_str)

    num_ref = num if num <= 7 else 14 - num
    states_ref = STATES[f"f{num_ref:02d}"]
    assert states_ref["source"] in SOURCES
    assert states_ref["l"] == 3
    assert states_ref["num"] == num_ref
    states_ref = set(states_ref["states"])

    if states != states_ref:
        print(sorted(states - states_ref))
        print(sorted(states_ref - states))
    assert states == states_ref
