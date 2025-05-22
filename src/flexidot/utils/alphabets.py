def alphabets(type_nuc=True):
    """
    provide ambiguity code for sequences
    """

    nucleotide_alphabet = ['A', 'C', 'G', 'T']

    nucleotide_alphabet_full = [
        'A',
        'C',
        'G',
        'T',
        'N',
        'B',
        'D',
        'H',
        'V',
        'Y',
        'R',
        'W',
        'S',
        'K',
        'M',
    ]

    nucleotide_ambiguity_code = {
        'N': ['A', 'C', 'G', 'T'],  # any
        'B': ['C', 'G', 'T'],  # not A
        'D': ['A', 'G', 'T'],  # not C
        'H': ['A', 'C', 'T'],  # not G
        'V': ['A', 'C', 'G'],  # not T
        'Y': ['C', 'T'],  # pyrimidine
        'R': ['A', 'G'],  # purine
        'W': ['A', 'T'],  # weak
        'S': ['C', 'G'],  # strong
        'K': ['G', 'T'],  # keto
        'M': ['A', 'C'],
    }  # amino

    nucleotide_match_dict = {
        'N': '[ACGTNBDHVYRWSKM]',  # any
        'B': '[CGTNBDHVYRWSKM]',  # not A
        'D': '[AGTNBDHVYRWSKM]',  # not C
        'H': '[ACTNBDHVYRWSKM]',  # not G
        'V': '[ACGNBDHVYRWSKM]',  # not T
        'K': '[GTNBDHVYRWSK]',  # keto - not A,C,M
        'M': '[ACNBDHVYRWSM]',  # amino - not G,T,K
        'W': '[ATNBDHVYRWKM]',  # weak - not C,G,S
        'S': '[CGNBDHVYRSKM]',  # strong - not A,G,W
        'Y': '[CTNBDHVYWSKM]',  # pyrimidine - not A,G,R
        'R': '[AGNBDHVRWSKM]',  # purine - not C,T,Y
        'A': '[ANDHVRWM]',
        'C': '[CNBHVYSM]',
        'G': '[GNBDVRSK]',
        'T': '[TNBDHYWK]',
    }

    aminoacid_alphabet = [
        'A',
        'R',
        'N',
        'D',
        'C',
        'E',
        'Q',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V',
        'U',
        'O',
        '*',
    ]

    aminoacid_alphabet_full = [
        'A',
        'R',
        'N',
        'D',
        'C',
        'E',
        'Q',
        'G',
        'H',
        'I',
        'L',
        'K',
        'M',
        'F',
        'P',
        'S',
        'T',
        'W',
        'Y',
        'V',
        'U',
        'O',
        '*',
        'J',
        'Z',
        'B',
        'X',
    ]

    aminoacid_ambiguity_code = {
        'J': ['I', 'L'],
        'Z': ['Q', 'E'],
        'B': ['N', 'D'],
        'X': [
            'A',
            'R',
            'N',
            'D',
            'C',
            'E',
            'Q',
            'G',
            'H',
            'I',
            'L',
            'K',
            'M',
            'F',
            'P',
            'S',
            'T',
            'W',
            'Y',
            'V',
            'U',
            'O',
            '*',
        ],
    }  # any

    aminoacid_match_dict = {
        'J': '[ILJ]',
        'Z': '[QEZ]',
        'B': '[NDB]',
        # "X": ".",
        'X': '[ARNDCEQGHILKMFPSTWYVUO*XBZJ]',
        'A': '[AX]',
        'R': '[RX]',
        'N': '[NXB]',
        'D': '[DXB]',
        'C': '[CX]',
        'E': '[EXZ]',
        'Q': '[QXZ]',
        'G': '[GX]',
        'H': '[HX]',
        'I': '[IXJ]',
        'L': '[LXJ]',
        'K': '[KX]',
        'M': '[MX]',
        'F': '[FX]',
        'P': '[PX]',
        'S': '[SX]',
        'T': '[TX]',
        'W': '[WX]',
        'Y': '[YX]',
        'V': '[VX]',
        'U': '[UX]',
        'O': '[OX]',
        '*': '[*X]',
    }

    # aa_only = {'E', 'F', 'I', 'J', 'L', 'O', 'Q', 'P', 'U', 'X', 'Z', '*'}
    # return nucleotide_alphabet, nucleotide_alphabet_full, nucleotide_ambiguity_code, aminoacid_alphabet, aminoacid_alphabet_full, aminoacid_ambiguity_code, aa_only

    if type_nuc:
        return (
            nucleotide_alphabet,
            nucleotide_alphabet_full,
            nucleotide_ambiguity_code,
            nucleotide_match_dict,
        )
    else:
        return (
            aminoacid_alphabet,
            aminoacid_alphabet_full,
            aminoacid_ambiguity_code,
            aminoacid_match_dict,
        )
