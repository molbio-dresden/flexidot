###############################
#     Matching Functions      #
###############################

import logging
import time

import regex
import numpy as np

from flexidot.utils.utils import time_track
from flexidot.utils.alphabets import alphabets
from flexidot.utils.analysis import (
    wobble_replacement,
    split_diagonals,
    lcs_from_x_values,
)


def find_match_pos_diag(
    seq1,
    seq2,
    wordsize,
    report_lcs=False,
    rc_option=True,
    convert_wobbles=False,
    max_N_percentage=49,
    type_nuc=True,
):
    """
    find all matching positions with matches >= wordsize
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """
    t1 = time.time()  # timer

    # look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc == True:
        any_residue = "N"
    else:
        any_residue = "X"

    # read sequences
    seq_one = seq1.upper()
    len_one = len(seq_one)
    seq_two = seq2.upper()
    len_two = len(seq_two)

    # set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[
        2
    ]  # nucleotide_ambiguity_code or aminoacid_ambiguity_code

    # forward
    #################################
    kmer_pos_dict_one = {}
    kmer_pos_dict_two = {}  # dictionaries for both sequences

    # reverse complement
    #################################
    kmer_pos_dict_three = {}
    kmer_pos_dict_four = {}  # dictionaries for both sequences

    # create dictionaries with kmers (wordsize) and there position(s) in the sequence
    if rc_option:
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
            (str(seq_two), kmer_pos_dict_two),
            (str(seq_one), kmer_pos_dict_three),
            (str(seq_two.reverse_complement()), kmer_pos_dict_four),
        ]
    else:
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
            (str(seq_two), kmer_pos_dict_two),
        ]
    for seq, kmer_pos_dict in data_list:
        for i in range(len(seq) - wordsize + 1):
            kmer = seq[i : i + wordsize]
            # discard kmer, if too many Ns included
            if kmer.count(any_residue) * 100.0 / wordsize <= max_N_percentage:
                if not convert_wobbles:
                    try:
                        kmer_pos_dict[kmer].append(i)
                    except KeyError:
                        kmer_pos_dict[kmer] = [i]
                else:
                    wobbles = False
                    for item in list(general_ambiguity_code.keys()):
                        if item in kmer:
                            wobbles = True
                            break
                    if not wobbles:
                        try:
                            kmer_pos_dict[kmer].append(i)
                        except KeyError:
                            kmer_pos_dict[kmer] = [i]
                    else:
                        kmer_variants = wobble_replacement(kmer, general_ambiguity_code)
                        for new_kmer in kmer_variants:
                            # print "\t", new_kmer
                            try:
                                kmer_pos_dict[new_kmer].append(i)
                            except KeyError:
                                kmer_pos_dict[new_kmer] = [i]

    # find kmers shared between both sequences
    matches_for = set(kmer_pos_dict_one).intersection(kmer_pos_dict_two)  # forward
    matches_rc = set(kmer_pos_dict_three).intersection(
        kmer_pos_dict_four
    )  # reverse complement

    text = "[matches: %i for; %.i rc]" % (len(matches_for), len(matches_rc))
    logging.debug(text)

    # create lists of x and y co-ordinates for scatter plot
    # keep all coordinates of all shared kmers (may match multiple times)
    diag_dict_for = {}
    diag_dict_rc = {}
    for match_list, pos_dict1, pos_dict2, diag_dict in [
        (matches_for, kmer_pos_dict_one, kmer_pos_dict_two, diag_dict_for),
        (matches_rc, kmer_pos_dict_three, kmer_pos_dict_four, diag_dict_rc),
    ]:
        for kmer in match_list:
            for i in pos_dict1[kmer]:
                for j in pos_dict2[kmer]:
                    diag = i - j
                    points = set(range(i + 1, i + wordsize + 1))
                    if diag not in list(diag_dict.keys()):
                        diag_dict[diag] = points
                    else:
                        diag_dict[diag].update(points)

    # convert coordinate points to line start and stop positions
    x1 = []  # x values reverse
    y1 = []  # y values forward
    for diag in list(diag_dict_for.keys()):
        x_values = np.array(sorted(diag_dict_for[diag]))
        x1.extend(split_diagonals(x_values))
        y_values = split_diagonals(x_values - diag)
        y1.extend(y_values)

    x2 = []  # x values rc
    y2 = []  # y values rc
    if rc_option:
        for diag in list(diag_dict_rc.keys()):
            factor = len_two + diag + 1
            x_values = np.array(sorted(diag_dict_rc[diag]))
            x2.extend(split_diagonals(x_values))
            y_values = split_diagonals(factor - x_values, -1)
            y2.extend(y_values)

    t1 = time_track(t1)

    if not report_lcs:
        return (
            np.array([np.array(x) for x in x1], dtype=object),
            np.array([np.array(y) for y in y1], dtype=object),
            np.array([np.array(x) for x in x2], dtype=object),
            np.array([np.array(y) for y in y2], dtype=object),
        )
    else:
        # get length of longest common substring based on match lengths
        lcs_for = lcs_from_x_values(x1)
        lcs_rev = lcs_from_x_values(x2)
        return (
            np.array([np.array(x) for x in x1], dtype=object),
            np.array([np.array(y) for y in y1], dtype=object),
            np.array([np.array(x) for x in x2], dtype=object),
            np.array([np.array(y) for y in y2], dtype=object),
            lcs_for,
            lcs_rev,
        )


def find_match_pos_regex(
    seq1,
    seq2,
    wordsize,
    substitution_count=0,
    report_lcs=False,
    rc_option=True,
    convert_wobbles=False,
    max_N_percentage=49,
    type_nuc=True,
):
    """
    find all matching positions with matches >= wordsize via regular expression search
    fuzzy matching - allow up to substitution_count substitutions
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """
    global t1  # timer

    # read sequences
    seq_one = seq1.upper()
    len_one = len(seq_one)
    seq_two = seq2.upper()
    len_two = len(seq_two)

    # set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[
        2
    ]  # nucleotide_ambiguity_code or aminoacid_ambiguity_code
    ambiguity_match_dict = alphabets(type_nuc)[3]

    ambiq_residues = "[%s]" % "".join(list(general_ambiguity_code.keys()))

    # look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc == True:
        any_residue = "N"
    else:
        any_residue = "X"

    # check for wobble presence
    if not (
        regex.search(ambiq_residues, str(seq_one)) == None
        and regex.search(ambiq_residues, str(seq_two)) == None
    ):
        wobble_found = True
    else:
        wobble_found = False

    # dictionary for matches
    diag_dict_for = {}
    diag_dict_rc = {}
    counter = [0, 0]

    # one-way matching
    if rc_option:
        data_list = [
            (str(seq_one), str(seq_two), diag_dict_for, 0),
            (str(seq_one), str(seq_two.reverse_complement()), diag_dict_rc, 1),
        ]
    else:
        data_list = [(str(seq_one), str(seq_two), diag_dict_for, 0)]

    for seq_query, seq_target, diag_dict, counter_pos in data_list:
        # split query sequence into kmers
        if not rc_option and counter_pos == 1:
            break

        for idx in range(len(str(seq_query)) - wordsize + 1):
            kmer = str(seq_query)[idx : idx + wordsize]

            # skip excessive N/X stretches (big black areas)
            if kmer.count(any_residue) * 100.0 / wordsize <= max_N_percentage:
                #  convert kmer to regular expression for wobble_matching
                if convert_wobbles and wobble_found:
                    kmer_string = ""
                    # replace each residue with matching residues or wobbles
                    for jdx in range(len(kmer)):
                        kmer_string += ambiguity_match_dict[kmer[jdx]]
                else:
                    kmer_string = kmer

                # convert to regular expression tolerating substitution errors
                if type(substitution_count) is int and substitution_count != 0:
                    kmer_string = "(%s){s<=%d}" % (kmer_string, substitution_count)

                # search for regular expression in target sequence
                kdx = 0
                start = True
                if regex.search(kmer_string, seq_target[kdx:]) != None:
                    counter[counter_pos] += 1
                    while regex.search(kmer_string, seq_target[kdx:]) != None:
                        # search for regular expression pattern in target sequence
                        result = regex.search(kmer_string, seq_target[kdx:])

                        kmer2 = seq_target[kdx:][result.start() : result.end()]

                        # skip excessive N/X stretches (big black areas)
                        if (
                            kmer2.count(any_residue) * 100.0 / wordsize
                            <= max_N_percentage
                        ):
                            diag = idx - (kdx + result.start())
                            points = set(range(idx + 1, idx + wordsize + 1))
                            if diag not in list(diag_dict.keys()):
                                diag_dict[diag] = points
                            else:
                                diag_dict[diag].update(points)

                        kdx += result.start() + 1
                        if kdx >= len(seq_target):
                            break
                        elif regex.search(kmer_string, seq_target[kdx:]) != None:
                            counter[counter_pos] += 1

    text = "%5.i \tforward matches" % counter[0]
    text += "\n%5.i \treverse complementary matches" % counter[1]
    logging.debug(text)

    # convert coordinate points to line start and stop positions
    x1 = []  # x values reverse
    y1 = []  # y values forward
    for diag in list(diag_dict_for.keys()):
        x_values = np.array(sorted(diag_dict_for[diag]))
        x1.extend(split_diagonals(x_values))
        y_values = split_diagonals(x_values - diag)
        y1.extend(y_values)

    x2 = []  # x values rc
    y2 = []  # y values rc
    if rc_option:
        for diag in list(diag_dict_rc.keys()):
            factor = len_two + diag + 1
            x_values = np.array(sorted(diag_dict_rc[diag]))
            x2.extend(split_diagonals(x_values))
            y_values = split_diagonals(factor - x_values, -1)
            y2.extend(y_values)

    t1 = time_track(t1)

    if not report_lcs:
        return (
            np.array([np.array(x) for x in x1], dtype=object),
            np.array([np.array(y) for y in y1], dtype=object),
            np.array([np.array(x) for x in x2], dtype=object),
            np.array([np.array(y) for y in y2], dtype=object),
        )
    else:
        # get length of longest common substring based on match lengths
        lcs_for = lcs_from_x_values(x1)
        lcs_rev = lcs_from_x_values(x2)
        return (
            np.array([np.array(x) for x in x1], dtype=object),
            np.array([np.array(y) for y in y1], dtype=object),
            np.array([np.array(x) for x in x2], dtype=object),
            np.array([np.array(y) for y in y2], dtype=object),
            lcs_for,
            lcs_rev,
        )
