##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# Compare calculated energy levels with data from Carnalls 1989 paper
# based on the full set of radial integrals including crystal field
# parameters for the C2v symmetry group for the configurations f2.a, f3,
# f11 and f12 in LaF3.
#
# Stark spritting is reproduced correctly, but it seems that Carnalls
# calculations suffered from low numerical precision.
#
##########################################################################

import re
import numpy as np
import pytest

from lanthanide import Lanthanide, Coupling, HalfInt
from test_basic import correct
from data_crystal import SOURCES, RADIAL


def split_term(term):
    m = re.match(r"(^[0-9]+)([A-Z])([0-9]+(/2)?)(\(([1-9]+)\))?$", term)
    assert bool(m)
    S = int(m[1])
    L = m[2]
    J = HalfInt(int(m[3][:-2])) if m[3][-2:] == "/2" else int(m[3])
    serial = m[6]
    if serial is None:
        serial = 0
    else:
        serial = int(serial)
    return S, L, J, serial


def run_crystal(data, maxdev, show=False):

    # Data correction list
    corrections = data.get("correct", [])

    # SLJ coupling
    coupling = Coupling.SLJM

    # Number of f electrons
    assert "num" in data
    num = data["num"]

    # Energy levels
    assert isinstance(data["energies"], list)
    energies = list(data["energies"])
    correct("energies", energies, corrections)
    indices = np.argsort(energies)
    energies = sorted(map(float, energies))
    assert len(set(indices)) == len(indices)
    assert set(indices) == set(range(len(energies)))

    # Radial integrals
    assert "radial" in data
    assert isinstance(data["radial"], dict)
    radial = dict(data["radial"])
    correct("radial", radial, corrections)

    # Adjust energy offset
    radial["base"] = energies[0]

    # Collect SLJ term names
    assert "terms" in data
    assert isinstance(data["terms"], list)
    assert len(data["terms"]) == len(energies)
    terms = list(data["terms"])
    correct("terms", terms, corrections)
    terms = [terms[i] for i in indices]

    # Initialize Lanthanide object
    success = True
    with Lanthanide(num, radial=radial, coupling=coupling) as ion:

        # Equalise energy offset
        dev = []
        for i in range(len(energies)):
            if num % 2 == 0:
                state = ion.intermediate[i]
            else:
                state = ion.intermediate[2 * i]
            dev.append(state.energy - energies[i])
        radial["base"] -= sum(dev) / len(dev)
        ion.set_radial(radial)

        # Compare energy levels
        dev = []
        for i in range(len(energies)):
            if num % 2 == 0:
                state = ion.intermediate[i]
            else:
                state = ion.intermediate[2 * i]
            energy_ref = energies[i]
            state_ref = terms[i]
            energy_calc = state.energy
            state_calc = state.long(min_weight=0.01)
            diff = energy_calc - energy_ref
            dev.append(diff * diff)
            if show:
                print(f"{i:3d}  /  ref {energy_ref:6.0f} | {state_ref:16s} >   /{diff:6.0f}/   calc {energy_calc:6.0f} | {state_calc} >")

            # Compare term assignments
            state_ref = set(term.strip() for term in state_ref.split(","))
            state_calc = set()
            j = 0
            while len(state_calc) < len(state_ref):
                state_calc.add(state.states[j].short().split("[")[0].split("(")[0])
                j += 1
            if state_ref != state_calc:
                print(f'("terms", {list(indices).index(i)}, "{state_ref.pop()}", "{state_calc.pop()}"),')
                success = False
        dev = np.sqrt(sum(dev)) / len(dev)
        if show:
            print(f"rms deviation: {dev:.1f} cm-1")
        assert dev <= maxdev
    assert success

@pytest.mark.parametrize("name", [key for key in RADIAL if "blocking" not in RADIAL[key]])
def test_crystal(name):
    maxdev = {
        "Pr3+:LaF3": 1.2,
        "Nd3+:LaF3": 0.5,
        "Er3+:LaF3": 3.8,
        "Tm3+:LaF3": 11.0,
    }[name]

    # Select data set
    assert name in RADIAL
    data = RADIAL[name]

    # Test source link
    assert "source" in data
    assert data["source"] in SOURCES

    # Run test
    run_crystal(data, maxdev)
