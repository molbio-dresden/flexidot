###############################
#     Analysis Functions      #
###############################

import logging

import numpy as np


def wobble_replacement(sequence, general_ambiguity_code):
    """
    get all degenerated sequences for sequence with ambiguous residues
    (only residues considered that are keys in wobble_dictionary)
    """

    # get positions of ambiguous residues
    wobble_pos = []
    for idx in range(len(sequence)):
        letter = sequence[idx]
        if letter in list(general_ambiguity_code.keys()):
            wobble_pos.append(idx)

    text = '\t%d wobbles' % len(wobble_pos)
    logging.debug(text)

    # replace one wobble through each iteration by all possible residues
    # repeat if still wobbles in new kmers
    kmer_variants = [sequence]
    while True:
        text = '\t\t%d kmer variants' % len(kmer_variants)
        logging.debug(text)
        temp_kmers = set()
        for kmer in kmer_variants:
            for idx in wobble_pos:
                letter = kmer[idx]
                if letter in list(general_ambiguity_code.keys()):
                    for base in general_ambiguity_code[kmer[idx]]:
                        newkmer = kmer[:idx] + base + kmer[idx + 1 :]
                        temp_kmers.add(newkmer)
        wobble = False
        for kmer in temp_kmers:
            for idx in range(len(kmer)):
                letter = kmer[idx]
                if letter in list(general_ambiguity_code.keys()):
                    wobble = True
                    break
            if wobble:
                break
        kmer_variants = set(list(temp_kmers)[:])
        if not wobble:
            break

    return kmer_variants


def split_diagonals(data, stepsize=1):
    """
    split array if point difference exceeds stepsize
    data = sorted list of numbers
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, _x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    _x_longest = x
            else:
                m[x][y] = 0
    return longest


def lcs_from_x_values(x_values):
    """
    calculate length of longest common substring based on nested list of numbers
    """
    if len(x_values) == 0:
        return 0
    # get lengths of each subarray data
    lengths = np.array([len(i) for i in x_values])
    return max(lengths)
