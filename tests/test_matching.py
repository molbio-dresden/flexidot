from Bio.Seq import Seq

# Import the function to be tested
from flexidot.utils.matching import find_match_pos_diag


def test_basic_functionality_nucleotide():
    seq1 = Seq('ATCGATCG')
    seq2 = Seq('GATCGATC')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize, rc_option=False)
    assert len(x1) > 0
    assert len(y1) > 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_basic_functionality_protein():
    seq1 = Seq('MKVLY')
    seq2 = Seq('KVLYM')
    wordsize = 3
    x1, y1, x2, y2 = find_match_pos_diag(
        seq1, seq2, wordsize, type_nuc=False, rc_option=False
    )
    assert len(x1) > 0
    assert len(y1) > 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_reverse_complement_enabled():
    seq1 = Seq('ATCGATCG')
    seq2 = Seq('CGATCGAT')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize, rc_option=True)
    assert len(x1) > 0
    assert len(y1) > 0
    assert len(x2) > 0
    assert len(y2) > 0


def test_reverse_complement_disabled():
    seq1 = Seq('ATCGATCG')
    seq2 = Seq('CGATCGAT')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize, rc_option=False)
    assert len(x1) > 0
    assert len(y1) > 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_ambiguous_residues_no_conversion():
    seq1 = Seq('ATCGNNNN')
    seq2 = Seq('NNNNGATC')
    wordsize = 4
    # When convert_wobbles is False, the function should skip any kmers with ambiguous residues
    x1, y1, x2, y2 = find_match_pos_diag(
        seq1, seq2, wordsize, convert_wobbles=False, max_N_percentage=60
    )
    assert len(x1) == 0
    assert len(y1) == 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_ambiguous_residues_with_conversion():
    seq1 = Seq('ATCGNNNN')
    seq2 = Seq('NNNNGATC')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(
        seq1, seq2, wordsize, convert_wobbles=True, max_N_percentage=26
    )
    # One rev match ATCG -> NGAT
    assert len(x1) == 0
    assert len(y1) == 0
    assert len(x2) == 1
    assert len(y2) == 1


def test_report_lcs_nucleotide():
    seq1 = Seq('ATCGATCG')
    seq2 = Seq('GATCGATC')
    wordsize = 4
    x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(
        seq1, seq2, wordsize, report_lcs=True
    )
    assert lcs_for == 7
    assert lcs_rev == 7


def test_report_lcs_protein():
    seq1 = Seq('MKVLY')
    seq2 = Seq('KVLYM')
    wordsize = 3
    x1, y1, x2, y2, lcs_for, lcs_rev = find_match_pos_diag(
        seq1, seq2, wordsize, report_lcs=True, type_nuc=False
    )
    assert lcs_for > 0
    assert lcs_rev == 0


def test_short_sequences():
    seq1 = Seq('ATC')
    seq2 = Seq('GAT')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize)
    assert len(x1) == 0
    assert len(y1) == 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_all_ambiguous_residues():
    seq1 = Seq('NNNN')
    seq2 = Seq('NNNN')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize)
    assert len(x1) == 0
    assert len(y1) == 0
    assert len(x2) == 0
    assert len(y2) == 0


def test_empty_sequences():
    seq1 = Seq('')
    seq2 = Seq('')
    wordsize = 4
    x1, y1, x2, y2 = find_match_pos_diag(seq1, seq2, wordsize)
    assert len(x1) == 0
    assert len(y1) == 0
    assert len(x2) == 0
    assert len(y2) == 0
