##########################################################################
# Copyright (c) 2025 Reinhard Caspary                                    #
# <reinhard.caspary@phoenixd.uni-hannover.de>                            #
# This program is free software under the terms of the MIT license.      #
##########################################################################
#
# This module provides an interface supporting the calculation of the
# matrices of tensor operators in the determinantal product space of a
# given electron configuration. It exploits the fact that these matrices
# are necessarily sparse because all matrix elements between states which
# differ in more electrons than the respective tensor operator is acting
# on, must be zero.
#
# The data structure of the interface contains a list of potential
# non-zero matrix elements and a list of binary bra-ket keys for the
# evaluation of one-, two, or three-electron elementary tensor operators.
# Conversion of electron tuples to binary keys reduces the disk space
# significantly and increases the calculation speed.
#
# No quantum numbers are used in this module. Electrons are
# characterised by their index into a list of indistinguishable
# electrons in standard order. In case of the 4f electron configurations
# of the lanthanides this means that the electron index my range
# from 0-13.
#
##########################################################################

from itertools import combinations
from functools import reduce

# Version of the algorithm. If the precomputed lists of indices in the file cache comes with another version number,
# it will be recomputed. Since this is the first element in the chain of dependent cache elements, changing this
# number will effectively render all file caches invalid. These elements are the modules 'unit', 'states', and
# 'matrix', see the init function of the Lanthanide class.
SINGLE_VERSION = 0

##########################################################################
# Binary state keys
##########################################################################

# Number of bits to encode one electron in a binary state key. There are 14 different f electrons, therefore
# four bits are sufficient.
BIT_LENGTH = 4

# Binary mask to extract a single electron from a key
INDEX_MASK = (1 << BIT_LENGTH) - 1

# Binary mask to extract the binary block of initial or final keys from a one-, two- or three-electron binary key
PART_MASK = [(1 << i * BIT_LENGTH) - 1 for i in range(1, 4)]

# Shift steps required to extract a single electron from the initial state encoded in a one-, two-, or
# three-electron bra-ket key
INITIAL_SHIFT = [(1 * BIT_LENGTH,), (3 * BIT_LENGTH, 2 * BIT_LENGTH), (5 * BIT_LENGTH, 4 * BIT_LENGTH, 3 * BIT_LENGTH)]

# Shift steps required to extract a single electron from the final state encoded in a one-, two-, or
# three-electron bra-ket key
FINAL_SHIFT = [(0,), (1 * BIT_LENGTH, 0), (2 * BIT_LENGTH, 1 * BIT_LENGTH, 0)]


def index_key(indices: tuple) -> int:
    """ Build a binary key from a tuple of electron numbers. with BIT_LENGTH bits for each electron. """

    # assert all(index < 1 << BIT_LENGTH for index in indices)
    return reduce(lambda num, index: (num << BIT_LENGTH) | index, indices, 0)


def braket_key(initial: tuple, final: tuple, parity: int) -> int:
    """ Build a binary bra-ket key with initial and final electrons in the higher bits. Parity is stored in bit 3
    and the number of electrons in the initial or final states in bits 0-2. """

    # assert parity in (0, 1)
    # assert len(initial) == len(final)
    # assert 1 <= len(initial) <= 3
    return (index_key(initial + final) << 1 | parity) << 3 | len(initial)


def braket_split_lower(key: int) -> (int, int):
    """ Return the electron number key part and the parity of a binary bra-ket key. In contrast to
    braket_split_upper(), the order of initial and final electrons is kept in the first return value. This
    is the default behaviour, because most tensor operators used in this package are symmetric and only the
    lower triangle is actually calculated. """

    return key >> 4, (key >> 3) & 1


def braket_split_upper(key: int) -> (int, int):
    """ Return the electron number key part and the parity of a binary bra-ket key. In contrast to
    braket_split_lower(), the initial and final electrons are swapped in the first return value to address
    matrix elements in the upper triangle. """

    num = key & 7
    parity = (key >> 3) & 1
    key >>= 4
    shift = num * BIT_LENGTH
    mask = PART_MASK[num - 1]
    initial = (key >> shift) & mask
    final = key & mask
    return (final << shift) | initial, parity


def key_pair(key: int, num: int):
    """ Extract the tuples of initial and final electrons from a binary bra-ket key without the last 4 bits as it
    is delivered in the first return value of the functions braket_split_lower() or braket_split_upper(). The number
    of electrons in each state is expected in the parameter num. """

    initial = tuple((key >> shift) & INDEX_MASK for shift in INITIAL_SHIFT[num - 1])
    final = tuple((key >> shift) & INDEX_MASK for shift in FINAL_SHIFT[num - 1])
    return initial, final


##########################################################################
# Generation of bra-ket lists for elementary tensor operators
##########################################################################

def determinant_one(brakets, electron_initial, electron_final, parity):
    """ Determinant function for one-electron elementary states. The determinant of a one-electron state is just
    the electron itself. The function therefore just appends the bra-ket key of the final and initial electrons. """

    brakets.append(braket_key(electron_initial, electron_final, parity))


def single_one(elements):
    """ Take lists of matrix elements with states differing in 0-3 electrons. Build and return one list of all
    matrix elements, which might be non-zero in case of a one-electron tensor operator and one list of the bra-ket
    keys of all determinantal electron permutations for which an elementary one-electron tensor operator needs to be
    evaluated in order to calculate the actual value of each of the matrix elements. For a one-electron tensor
    operator only states differing in zero (diagonal) or one electron can lead to non-zero matrix elements. All
    others are discarded. """

    # Initialize matrix elements index list and evaluation bra-ket keys list
    single_elements = []
    single_brakets = []

    # For each matrix element with no different electron (diagonal element), the elementary operator needs to be
    # evaluated for every electron from the common final and initial state.
    for state_index, state_electrons in enumerate(elements[0]):
        brakets = []
        for electron in state_electrons:
            determinant_one(brakets, (electron,), (electron,), 0)
        single_elements.append([state_index, state_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with one different electron in the final and initial states, the elementary operator
    # needs to be evaluated only for the different electrons in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[1]:
        brakets = []
        determinant_one(brakets, initial_electrons, final_electrons, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # Return matrix elements index list and evaluation bra-ket keys list
    return single_elements, single_brakets


def determinant_two(brakets, electrons_a_initial, electrons_a_final, parity):
    """ Determinant function for two-electron elementary states. The function appends the bra-ket keys of the four
    determinantal permutations of final and initial electrons. Note that the value of an elementary operator acting
    on the same quantity of two electrons like u1(k)u2(k) will be the same for the key pairs (0, 3) and (1, 2). This
    symmetry is not exploited yet, because it would increase the code complexity significantly and each tensor operator
    is calculated only once anyway and then taken from the file cache. """

    # Build order permutations
    electrons_b_initial = (electrons_a_initial[1], electrons_a_initial[0])
    electrons_b_final = (electrons_a_final[1], electrons_a_final[0])

    # Parity values of even- and odd-order permutations
    pos_parity = parity % 2
    neg_parity = (parity + 1) % 2

    # Append bra-ket keys of all permutations
    brakets.append(braket_key(electrons_a_initial, electrons_a_final, pos_parity))
    brakets.append(braket_key(electrons_a_initial, electrons_b_final, neg_parity))
    brakets.append(braket_key(electrons_b_initial, electrons_a_final, neg_parity))
    brakets.append(braket_key(electrons_b_initial, electrons_b_final, pos_parity))


def single_two(elements):
    """ Take lists of matrix elements with states differing in 0-3 electrons. Build and return one list of all
    matrix elements, which might be non-zero in case of a two-electron tensor operator and one list of the bra-ket
    keys of all determinantal electron permutations for which an elementary two-electron tensor operator needs to be
    evaluated in order to calculate the actual value of each of the matrix elements. For a two-electron tensor
    operator only states differing in zero (diagonal), one, or two electrons can lead to non-zero matrix elements.
    All others are discarded. """

    # Initialize matrix elements index list and evaluation bra-ket keys list
    single_elements = []
    single_brakets = []

    # For each matrix element with no different electron (diagonal element), the elementary operator needs to be
    # evaluated for every combination of two electrons from the common final and initial state.
    for state_index, state_electrons in enumerate(elements[0]):
        brakets = []
        for electrons in combinations(state_electrons, 2):
            determinant_two(brakets, electrons, electrons, 0)
        single_elements.append([state_index, state_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with one different electron in the final and initial states, the elementary operator
    # needs to be evaluated for each pair of one electron from the tuple of same electrons and the different
    # electron in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[1]:
        brakets = []
        for same in same_electrons:
            initial = (same,) + initial_electrons
            final = (same,) + final_electrons
            determinant_two(brakets, initial, final, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with two different electrons in the final and initial states, the elementary operator
    # needs to be evaluated only for the pair of different electrons in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[2]:
        brakets = []
        determinant_two(brakets, initial_electrons, final_electrons, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # Return matrix elements index list and evaluation bra-ket keys list
    return single_elements, single_brakets


def determinant_three(brakets, electrons_initial, electrons_final, element_parity):
    """ Determinant function for three-electron elementary states. The function appends the bra-ket keys of the 36
    determinantal permutations of final and initial electrons. Note that the value of an elementary operator acting
    on the same quantity of two electrons like u1(k)u2(k)u3(k') or three electrons like u1(k)u2(k)u3(k) will be the
    same for certain key triples. This symmetry is not exploited yet, because it would increase the code complexity
    significantly and each tensor operator is calculated only once anyway and then taken from the file cache. """

    # Build order permutations of the initial states
    initial_electrons = [
        electrons_initial,
        (electrons_initial[2], electrons_initial[1], electrons_initial[0]),
        (electrons_initial[1], electrons_initial[2], electrons_initial[0]),
        (electrons_initial[1], electrons_initial[0], electrons_initial[2]),
        (electrons_initial[2], electrons_initial[0], electrons_initial[1]),
        (electrons_initial[0], electrons_initial[2], electrons_initial[1])]

    # Build order permutations of the final states
    final_electrons = [
        electrons_final,
        (electrons_final[2], electrons_final[1], electrons_final[0]),
        (electrons_final[1], electrons_final[2], electrons_final[0]),
        (electrons_final[1], electrons_final[0], electrons_final[2]),
        (electrons_final[2], electrons_final[0], electrons_final[1]),
        (electrons_final[0], electrons_final[2], electrons_final[1])]

    # Append bra-ket keys of all permutations
    for i, initial in enumerate(initial_electrons):
        for j, final in enumerate(final_electrons):
            parity = (element_parity + i + j) % 2
            key = braket_key(initial, final, parity)
            brakets.append(key)


def single_three(elements):
    """ Take lists of matrix elements with states differing in 0-3 electrons. Build and return one list of all
    matrix elements, which might be non-zero in case of a three-electron tensor operator and one list of the bra-ket
    keys of all determinantal electron permutations for which an elementary three-electron tensor operator needs to be
    evaluated in order to calculate the actual value of each of the matrix elements. For a three-electron tensor
    operator only states differing in zero (diagonal), one, two, or three  electrons can lead to non-zero matrix
    elements. All others are discarded. """

    # Initialize matrix elements index list and evaluation bra-ket keys list
    single_elements = []
    single_brakets = []

    # For each matrix element with no different electron (diagonal element), the elementary operator needs to be
    # evaluated for every combination of three electrons from the common final and initial state.
    for state_index, state_electrons in enumerate(elements[0]):
        brakets = []
        for electrons in combinations(state_electrons, 3):
            determinant_three(brakets, electrons, electrons, 0)
        single_elements.append([state_index, state_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with one different electron in the final and initial states, the elementary operator
    # needs to be evaluated for the triple of each combination of two electrons from the tuple of same electrons
    # together with the different electron in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[1]:
        brakets = []
        for same in combinations(same_electrons, 2):
            initial = same + initial_electrons
            final = same + final_electrons
            determinant_three(brakets, initial, final, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with two different electrons in the final and initial states, the elementary operator
    # needs to be evaluated for each triple of one electron from the tuple of same electrons together with the
    # different electrons in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[2]:
        brakets = []
        for same in same_electrons:
            initial = (same,) + initial_electrons
            final = (same,) + final_electrons
            determinant_three(brakets, initial, final, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # For each matrix element with three different electrons in the final and initial states, the elementary operator
    # needs to be evaluated only for the triple of different electrons in the final and initial states.
    for initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity in elements[3]:
        brakets = []
        determinant_three(brakets, initial_electrons, final_electrons, parity)
        single_elements.append([initial_index, final_index, len(brakets)])
        single_brakets += brakets

    # Return matrix elements index list and evaluation bra-ket keys list
    return single_elements, single_brakets


##########################################################################
# Data structure generation
##########################################################################

def matrix_elements(states):
    """ Generate all potentially non-zero matrix elements. Return number of different electrons, indices of the
    matrix element, same electrons, different electrons in the initial and final states, and parity. """

    # Check every matrix element if it might be non-zero
    for final_index, (final_state, final_diff, final_keys) in enumerate(states):
        for initial_index, (initial_state, initial_diff, initial_keys) in enumerate(states[:final_index]):

            # Initial and final state of diagonal elements always share all electrons. This case is treated separately.
            if initial_index == final_index:
                continue

            # Find the minimum number of different electrons between initial and final state of the matrix element
            # in the range 1-3.
            for num in range(3):

                # There might be at most one match of the same electrons
                match = initial_keys[num] & final_keys[num]
                if len(match) > 1:
                    raise RuntimeError("Multiple state matches!")

                # Initial and final state differ in num+1 electrons
                elif len(match) == 1:

                    # Get tuples of same electrons, different initial and final electrons and swap numbers
                    key = match.pop()
                    same_electrons, initial_electrons, initial_swaps = initial_diff[num][key]
                    _, final_electrons, final_swaps = final_diff[num][key]

                    # Calculate parity based on the sum of swaps
                    parity = (initial_swaps + final_swaps) % 2

                    # Return indices of the matrix element, same and different electrons and parity
                    data = initial_index, final_index, same_electrons, initial_electrons, final_electrons, parity
                    yield num + 1, data

                    # Deliver only the match with the smallest number of different electrons
                    break


def diff_electrons(diff, state):
    """ Build and return a diff-dictionary for a given state. Each electron in the state tuple is represented by an
    index to one of 2*(2l+1) different electrons. Take the tuple of num ordered electrons of the state and determine
    every possible splitting into the ordered tuples 'same' containing num-diff electrons and 'other' with diff
    electrons. Store the two tuples and the number of swap operations required for the splitting in the returned
    diff-dictionary. Use a binary key which identifies the tuple 'same'. The set of keys of the diff-dictionary
    provides an easy way to determine if two states differ in exactly diff electrons or not. """

    # Number of electrons
    num_electrons = len(state)

    # Build dictionary with all combinations of diff 'other' electrons
    diff_dict = {}
    if diff <= num_electrons:

        # Go through all combinations of diff electrons picked from the initial tuple
        indices = range(num_electrons)
        for selected in combinations(indices, diff):
            # num-diff electrons remaining in the 'same' tuple and keeping their order
            same = tuple(state[i] for i in indices if i not in selected)

            # Build a binary key from the 'same' tuple
            key = index_key(same)
            # assert key not in diff_dict

            # diff electrons in the 'other' tuple, also keeping their order
            other = tuple(state[i] for i in selected)

            # Number of neighbor swaps required to move the 'other' electrons to the end of the initial tuple
            swaps = diff * num_electrons - sum(range(diff + 1)) - sum(selected)

            # Store split tuples and the swap number in the dictionary
            diff_dict[key] = (same, other, swaps % 2)

    # Return the dictionary
    return diff_dict


def single_elements(states):
    """ Determine all matrix elements which might be non-zero for a given set of states. Return lists of
    parameters for elementary one-, two-, or three-electron tensor operators to be evaluated for the respective
    matrix elements. """

    # Collect comparison data for each product state
    diff_states = []
    for electrons in states:
        # Get the diff-dictionaries for 1, 2, or 3 'other' electrons
        diff_dicts = [diff_electrons(diff, electrons) for diff in range(1, 4)]

        # Build sets of the 'same' key numbers
        diff_keys = [set(diff.keys()) for diff in diff_dicts]

        # Store state, three diff-dictionaries and the three respective sets of keys
        diff_states.append((electrons, diff_dicts, diff_keys))

    # Initialize a list of length 4, containing lists of all matrix elements (state pairs) which differ in exactly
    # 0-3 electrons. Diagonal matrix elements connect identical states, which differ in 0 electrons. Thus, the first
    # element contains the list of states
    elements = [states, [], [], []]

    # Collect all matrix elements with a minimum diff of 1, 2, or 3 electrons. All other elements are discarded.
    for diff_num, data in matrix_elements(diff_states):
        elements[diff_num].append(data)

    # Number of electrons in the configuration
    num_electrons = len(states[0])

    # Determine parameters for all elementary one-electron tensor operators to be evaluated for each
    # potential matrix element (diff is 1)
    one = single_one(elements) if num_electrons >= 1 else ([], [])

    # Determine parameters for all elementary two-electron tensor operators to be evaluated for each
    # potential matrix element (diff is 1 or 2)
    two = single_two(elements) if num_electrons >= 2 else ([], [])

    # Determine parameters for all elementary three-electron tensor operators to be evaluated for each
    # potential matrix element (diff is 1, 2, or 3)
    three = single_three(elements) if num_electrons >= 3 else ([], [])

    # Return evaluation parameters for elementary one-, two-, or three-electron tensor operators
    return one, two, three


##########################################################################
# HDF5 cache interface
##########################################################################

class SingleElements:
    """ This class uses the precalculated lists of non-zero matrix elements and the respective lists of binary
    bra-ket keys from single_elements() to provide support for the calculation of the matrix of a tensor operator
    in the determinantal product space. """

    def __init__(self, group, states):

        # HDF5 group which contains the precalculated data. It contains the three groups 'one', 'two', and 'three'
        # with data for one-, two-, or three-electron tensor operators. Each of these groups contains the two
        # datasets 'indices' with the indices of potentially non-zero matrix elements and 'elements' with the
        # list of bra-ket keys for the elementary operators required for the calculation of these matrix elements.
        self.group = group

        # Electron states of the current configuration
        self.states = states

        # Number of electrons in the current configuration
        self.num_electrons = len(states[0])

    def elements(self, num: int):
        """ Generate the indices of the initial and the final state of each potentially non-zero matrix element
        of a num-electron tensor operator together with a slice into the list of binary bra-ket keys for the
        evaluation of elementary tensor operators required for this matrix element. """

        index = 0
        for initial_index, final_index, size in self[num]["indices"]:
            yield initial_index, final_index, slice(index, index + size)
            index += size

    def lower_keys(self, key_slice: slice, num: int):
        """ Generate the electron number key part and the parity of each binary bra-ket key in the range given by
        the slice for the matrix element of a num-electron tensor operator in the lower triangle. """

        for key in self[num]["elements"][key_slice]:
            yield braket_split_lower(int(key))

    def upper_keys(self, key_slice: slice, num: int):
        """ Generate the electron number key part and the parity of each binary bra-ket key in the range given by
        the slice for the matrix element of a num-electron tensor operator in the upper triangle. """

        for key in self[num]["elements"][key_slice]:
            yield braket_split_upper(int(key))

    def index_pair(self, key: int, num: int) -> (tuple, tuple):
        """ Return tuples of the num initial and the num final electrons from a binary bra-ket key without the
        last 4 bits as it is generated as the first return value by the generators lower_keys() or upper_keys()."""

        return key_pair(key, num)

    def __len__(self):
        """ Number of electrons in the configuration. """

        return self.num_electrons

    def __getitem__(self, num):
        """ Return the list of potentially non-zero matrix elements of a num-electron tensor operator and the
        list of binary bra-ket keys for the evaluation of elementary tensor operators required for these
        matrix elements. """

        key = ["one", "two", "three"][num - 1]
        return self.group[key]

    def __iter__(self):
        """ Generate all pairs of matrix element lists and bra-ket key lists which are applicable to the current
        configuration based on its number of electrons. """

        if self.num_electrons >= 1:
            yield self.group["one"]
        if self.num_electrons >= 2:
            yield self.group["two"]
        if self.num_electrons >= 3:
            yield self.group["three"]

    def __str__(self):
        return f"SingleElements support object for a {self.num_electrons} electron configuration"


def init_single(vault, group_name, states):
    """ Initialize the cache for the support structure for the calculation of tensor operator matrix elements for
    the given determinantal product states in the HDF5 group with given name in the given HDF5 file vault. Return
    the SingleElements support object. """

    # No file cache
    if not vault:
        one, two, three = single_elements(states)
        group = {"one": {"indices": one[0], "elements": one[1]},
                 "two": {"indices": two[0], "elements": two[1]},
                 "three": {"indices": three[0], "elements": three[1]}}

    # Use file cache
    else:

        # Delete the group in the HDF5 file, if the cache is marked as invalid or its version number does not match
        if group_name in vault:
            if not vault.attrs["valid"] or "version" not in vault[group_name].attrs or \
                    vault[group_name].attrs["version"] != SINGLE_VERSION:
                del vault[group_name]
                vault.flush()

        # Generate all data, if the HDF5 group is missing
        if group_name not in vault:
            print("Creating vault single ...")

            # Render all cache structures following in the dependence chain as invalid
            vault.attrs["valid"] = False

            # Create new HDF5 group and store the current version number
            vault.create_group(group_name)
            vault[group_name].attrs["version"] = SINGLE_VERSION

            # Calculate the lists of matrix indices and bra-ket keys for one-, two-, and three-electron operators
            one, two, three = single_elements(states)

            # Store the lists of matrix indices and bra-ket keys for one-electron operators as HDF5 datasets
            group = vault[group_name].create_group("one")
            group.create_dataset("indices", data=one[0], compression="gzip", compression_opts=9)
            group.create_dataset("elements", data=one[1], compression="gzip", compression_opts=9)

            # Store the lists of matrix indices and bra-ket keys for two-electron operators as HDF5 datasets
            group = vault[group_name].create_group("two")
            group.create_dataset("indices", data=two[0], compression="gzip", compression_opts=9)
            group.create_dataset("elements", data=two[1], compression="gzip", compression_opts=9)

            # Store the lists of matrix indices and bra-ket keys for three-electron operators as HDF5 datasets
            group = vault[group_name].create_group("three")
            group.create_dataset("indices", data=three[0], compression="gzip", compression_opts=9)
            group.create_dataset("elements", data=three[1], compression="gzip", compression_opts=9)

            # Flush the cache file
            vault.flush()
            print("Vault single done.")

        # HDF5 single group
        group = vault[group_name]

    # Return an interface object supporting the calculation of the matrices of tensor operators in the determinantal
    # product space of the given configuration
    return SingleElements(group, states)
