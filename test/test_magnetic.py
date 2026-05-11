##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# Compare 1st order spin-spin and spin-other-orbit (H5) and 2nd order
# spin-orbit (H6) matrix elements from the literature for the f2.a and f12
# configurations with results from the Lanthanide package.
#
# Note: We use equation (3) in [24] for the comparison. This equation
# contains a typo. The correct sign factor is (-1)^(S'+L'+J) with L'.
#
##########################################################################

import itertools
import pytest

from lanthanide import Lanthanide, Coupling, wigner6j
from data_magnetic import SOURCES, MAGNETIC


@pytest.mark.parametrize("data_key", MAGNETIC.keys())
def test_magnetic(data_key):
    # Select data set
    assert data_key in MAGNETIC
    data = MAGNETIC[data_key]

    # Test source link
    assert "source" in data
    assert data["source"] in SOURCES

    # Number of f electrons
    assert "num" in data
    num = data["num"]

    # Rank of double tensor operators
    assert "t" in data
    t = data["t"]

    # Name of tensor operator
    assert "name" in data
    name = data["name"]

    # Reduced SL matrix elements
    assert "reduced" in data
    elements = data["reduced"]

    # Compare tensor operator matrix to the calculation
    success = True
    with Lanthanide(num) as ion:
        states = ion.states(Coupling.SLJ)
        array = ion.matrix(name, Coupling.SLJ).array

        # Initialize set of SL states and list of phase signs
        sl_states = set()
        phases = []

        # Compare every matrix element
        for i in range(array.shape[0]):
            for j in range(i + 1):

                # Quantum numbers S, L, J of final and initial states
                Sa = states[i]["S2"].S
                La = states[i]["L2"].L
                Ja = states[i]["J2"].J
                Sb = states[j]["S2"].S
                Lb = states[j]["L2"].L
                Jb = states[j]["J2"].J

                # Test for zero if Ja != Jb
                if Ja != Jb:
                    assert abs(array[i, j]) < 1e-9
                    continue

                # Factor of Wigner 6-j symbol
                factor = wigner6j(Sb, Lb, Ja, La, Sa, t)
                if abs(factor) < 1e-9:
                    assert abs(array[i, j]) < 1e-9
                    continue

                # Get respective reduced SL matrix element
                LSa = str(states[i]["S2"]) + str(states[i]["L2"])
                LSb = str(states[j]["S2"]) + str(states[j]["L2"])
                if f"{LSa} {LSb}" in elements:
                    reduced = elements[f"{LSa} {LSb}"]
                    element = f"< {LSa} | {LSb} >"
                elif f"{LSb} {LSa}" in elements:
                    reduced = elements[f"{LSb} {LSa}"]
                    element = f"< {LSb} | {LSa} >"
                else:
                    if not abs(array[i, j]) < 1e-9:
                        success = False
                        print(f"ERROR: unknown SL element < {LSa} | {LSb} >!")
                    continue

                # Calculate matrix element with sign
                # Note: [24] states -1 ^ ((Sb + La + Ja) % 2) for the sign, which is wrong!
                value = factor * reduced
                if (Sb + Lb + Ja) % 2 != 0:
                    value = -value

                # SL diagonal element, sign always +1
                if LSa == LSb:
                    if not abs(value - array[i, j]) < 1e-9:
                        success = False
                        print(f"ERROR: {name}[{i},{j}]: {element} {value} != {array[i, j]}")
                    continue

                # Both SL states with same sign
                if abs(value - array[i, j]) < 1e-9:
                    sl_states.add(LSa)
                    sl_states.add(LSb)
                    phases.append((LSa, LSb, 1))
                    continue

                # Both SL states with opposite sign
                if abs(value + array[i, j]) < 1e-9:
                    sl_states.add(LSa)
                    sl_states.add(LSb)
                    phases.append((LSa, LSb, -1))
                    continue

                # Different magnitude of given and calculated matrix element
                success = False
                print(f"ERROR: {name}[{i},{j}]: {element} {value} != {array[i, j]}")

        # Stop if one element test failed
        assert success

        # Test for state phases which match the detected matrix element signs
        sl_states = tuple(sl_states)
        success = False
        for signs in itertools.product((+1, -1), repeat=len(sl_states)):
            if all(signs[sl_states.index(LSa)] * signs[sl_states.index(LSb)] == sign for LSa, LSb, sign in phases):
                success = True
                break
        assert success
