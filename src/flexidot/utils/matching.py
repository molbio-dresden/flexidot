###############################
#     Matching Functions      #
###############################

from collections import defaultdict
import logging
from typing import Tuple, Union

from Bio.Seq import Seq
import numpy as np
import regex

# from flexidot.utils.utils import time_track
from flexidot.utils.alphabets import alphabets
from flexidot.utils.analysis import (
    lcs_from_x_values,
    split_diagonals,
    wobble_replacement,
)


def find_match_pos_diag(
    seq1: Union[str, Seq],
    seq2: Union[str, Seq],
    wordsize: int,
    report_lcs: bool = False,
    rc_option: bool = True,
    convert_wobbles: bool = False,
    max_N_percentage: float = 0,
    type_nuc: bool = True,
) -> Union[
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int],
]:
    """
    Find all matching positions with matches >= wordsize
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """
    # TODO: check if seq1 and se2 are the same sequence (e.g. in case of self-alignment)
    # If so, then can skip counting the fwd kmers for seq2 and recycle the results from seq1
    if seq1 == seq2:
        logging.debug('Self-alignment detected. Recycling results from seq1 for seq2.')
        self_alignment = True
    else:
        self_alignment = False

    # Look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc:
        unknown_residue = 'N'
    else:
        unknown_residue = 'X'

    # Calculate the maximum number of Ns allowed in a kmer
    max_N_count = (max_N_percentage / 100.0) * wordsize

    # read sequences
    seq_one = seq1.upper()
    # len_one = len(seq_one)
    seq_two = seq2.upper()
    len_two = len(seq_two)

    # Set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[
        2
    ]  # nucleotide_ambiguity_code or aminoacid_ambiguity_code

    # Forward
    #################################
    kmer_pos_dict_one = defaultdict(list)  # Seq1
    kmer_pos_dict_two = defaultdict(list)  # Seq2

    # Reverse complement
    #################################
    kmer_pos_dict_three = defaultdict(list)  # Seq1
    kmer_pos_dict_four = defaultdict(list)  # Seq2

    # Create dictionaries to index kmer (wordsize) positions in the sequence
    if self_alignment and rc_option:
        # Compare seq1 forward to self and self reverse
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
            (str(seq_one.reverse_complement()), kmer_pos_dict_three),
        ]
    elif self_alignment and not rc_option:
        # Compare seq1 forward to self only
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
        ]
    elif rc_option:
        # Compare seq1 forward to seq2 forward and reverse
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
            (str(seq_two), kmer_pos_dict_two),
            # (str(seq_one), kmer_pos_dict_three), #TODO: Check if this is needed
            (str(seq_two.reverse_complement()), kmer_pos_dict_four),
        ]
    else:
        # Compare seq1 and seq2 forward only
        data_list = [
            (str(seq_one), kmer_pos_dict_one),
            (str(seq_two), kmer_pos_dict_two),
        ]

    # Step through each sequence and add kmers to dictionary
    for seq, kmer_pos_dict in data_list:
        # Track number of kmers skipped due to Ns > max_N_count
        skipped_Ns = 0
        # Step through sequence and add kmer positions to dictionary
        for i in range(len(seq) - wordsize + 1):
            # Extract kmer
            kmer = seq[i : i + wordsize]
            # Count Ns in kmer
            Ns_in_kmer = kmer.count(unknown_residue)
            # Discard kmer, if too many Ns included
            if Ns_in_kmer <= max_N_count:
                if not convert_wobbles:
                    if Ns_in_kmer == 0:
                        # Add kmer to dictionary if no Ns.
                        kmer_pos_dict[kmer].append(i)
                    else:
                        # Skip kmers with Ns
                        skipped_Ns += 1
                else:  # Deal with ambiguous characters
                    # Set as True if any wobble characters are present in the kmer
                    has_wobbles = any(
                        item in kmer for item in general_ambiguity_code.keys()
                    )
                    if not has_wobbles:
                        kmer_pos_dict[kmer].append(i)
                    else:
                        # Replace wobble characters with all possible variants
                        kmer_variants = wobble_replacement(kmer, general_ambiguity_code)
                        for new_kmer in kmer_variants:
                            kmer_pos_dict[new_kmer].append(i)
            else:
                skipped_Ns += 1

        # Log number of skipped kmers
        if skipped_Ns > 0:
            if convert_wobbles:
                logging.debug(
                    'Skipped %i kmers due to {unknown_residue}s > %i'
                    % (skipped_Ns, max_N_count)
                )
            else:
                logging.debug(
                    'Skipped %i kmers containing {unknown_residue}s' % (skipped_Ns)
                )

    # If self alignment, duplicate self fwd and rev dictionaries
    if self_alignment:
        # If self alignment, copy kmer_pos_dict_one to kmer_pos_dict_two
        kmer_pos_dict_two = kmer_pos_dict_one.copy()
        # Copy kmer_pos_dict_three to kmer_pos_dict_four
        kmer_pos_dict_four = kmer_pos_dict_three.copy()

    # Find kmers shared between both sequences in forward orientation
    matches_for = set(kmer_pos_dict_one).intersection(kmer_pos_dict_two)

    # Find kmers shared between Seq1 forward and Seq2 reverse orientation
    matches_rc = set(kmer_pos_dict_one).intersection(kmer_pos_dict_four)

    # TODO: Check if above is correct. Do we gain any extra kmers over set(matches_for, matches_rc) if we also check the seq1_rc?
    # print("matches_for: ", type(matches_for) ,matches_for)
    # print("matches_rc: ", matches_rc)

    logging.debug(
        '[matches: %i forward; %i reverse]' % (len(matches_for), len(matches_rc))
    )

    # Create lists of x and y coordinates for scatter plot
    # Keep all coordinates of all shared kmers (may match multiple times)
    diag_dict_for = defaultdict(set)
    diag_dict_rc = defaultdict(set)

    for match_list, pos_dict1, pos_dict2, diag_dict in [
        (matches_for, kmer_pos_dict_one, kmer_pos_dict_two, diag_dict_for),
        (matches_rc, kmer_pos_dict_one, kmer_pos_dict_four, diag_dict_rc),
    ]:
        for kmer in match_list:
            for i in pos_dict1[kmer]:
                for j in pos_dict2[kmer]:
                    diag = i - j
                    points = set(range(i + 1, i + wordsize + 1))
                    diag_dict[diag].update(points)

    # Convert coordinate points to line start and stop positions
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

    if not report_lcs:
        return (
            # x values forward matches
            np.array([np.array(x) for x in x1], dtype=object),
            # y values forward matches
            np.array([np.array(y) for y in y1], dtype=object),
            # x values rc (ascending)
            np.array([np.array(x) for x in x2], dtype=object),
            # y values rc (descending)
            np.array([np.array(y) for y in y2], dtype=object),
        )
    else:
        # Get length of longest common substring based on match lengths
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
    max_N_percentage=0,
    type_nuc=True,
):
    """
    find all matching positions with matches >= wordsize via regular expression search
    fuzzy matching - allow up to substitution_count substitutions
    convert matching points into lines of the length of the match
    (+ optional handling of ambiguities)
    """

    # read sequences
    seq_one = seq1.upper()
    # len_one = len(seq_one)
    seq_two = seq2.upper()
    len_two = len(seq_two)

    # set ambiguity code for wobble replacement
    general_ambiguity_code = alphabets(type_nuc)[
        2
    ]  # nucleotide_ambiguity_code or aminoacid_ambiguity_code
    ambiguity_match_dict = alphabets(type_nuc)[3]

    ambiq_residues = '[%s]' % ''.join(list(general_ambiguity_code.keys()))

    # look for Ns in DNA or Xs in proeins (minimum word size)
    if type_nuc:
        any_residue = 'N'
    else:
        any_residue = 'X'

    # check for wobble presence
    if not (
        regex.search(ambiq_residues, str(seq_one)) is None
        and regex.search(ambiq_residues, str(seq_two)) is None
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
                    kmer_string = ''
                    # replace each residue with matching residues or wobbles
                    for jdx in range(len(kmer)):
                        kmer_string += ambiguity_match_dict[kmer[jdx]]
                else:
                    kmer_string = kmer

                # convert to regular expression tolerating substitution errors
                if type(substitution_count) is int and substitution_count != 0:
                    kmer_string = '(%s){s<=%d}' % (kmer_string, substitution_count)

                # search for regular expression in target sequence
                kdx = 0
                # start = True
                if regex.search(kmer_string, seq_target[kdx:]) is not None:
                    counter[counter_pos] += 1
                    while regex.search(kmer_string, seq_target[kdx:]) is not None:
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
                        elif regex.search(kmer_string, seq_target[kdx:]) is not None:
                            counter[counter_pos] += 1

    text = '%5.i \tforward matches' % counter[0]
    text += '\n%5.i \treverse complementary matches' % counter[1]
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
