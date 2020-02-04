#!/usr/bin/env python3
"""
Author:         Sam Overduin

Description:    This script is an implementation of the Needleman-Wunsch algorithm.
                It takes 2 protein sequences and aligns them
                Under __main__ examples are implemented according to the stated
                exercises.
"""

blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""


def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int, parts[1:])))
    return order, blosum_matrix


BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()


def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues

    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]


# write your own functions below here

def initialize_matrices(seq1, seq2, end_gap_penalty):
    """Initializes score_matrix and arrow_matrix by filling the first column and row.

    :param seq1: string of amino acids.
    :param seq2: string of amino acids.
    :param end_gap_penalty: int of end gap cost outside alignment

    :return: returns score_matrix, arrow_matrix
    """
    # Create score_matrix with seq2 horizontal and seq1 vertical
    score_matrix = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
    # Initialise arrow_matrix which shows where the score came from.
    arrow_matrix = [['' for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]

    # Fill first column and row with end gaps (start is also the end)
    for i in range(1, len(score_matrix)):
        score_matrix[i][0] = -end_gap_penalty * i
        arrow_matrix[i][0] = 'V'
    for j in range(1, len(score_matrix[0])):
        score_matrix[0][j] = -end_gap_penalty * j
        arrow_matrix[0][j] = 'H'

    return score_matrix, arrow_matrix

def make_align_matrices(seq1, seq2, gap_penalty, end_gap_penalty):
    """Calculates the score alignment matrix and keeps track of path with arrow_matrix

    :param seq1: string of amino acids.
    :param seq2: string of amino acids.
    :param gap_penalty: int of gap cost within alignment
    :param end_gap_penalty: int of end gap cost outside alignment

    :return: returns score_matrix, arrow_matrix
    """
    # Create score_matrix and arrow_matrix with seq2 horizontal and seq1 vertical
    score_matrix, arrow_matrix = initialize_matrices(seq1, seq2, end_gap_penalty)

    # i is Vertical, j is horizontal
    for i in range(1, len(score_matrix)):
        for j in range(1, len(score_matrix[i])):
            if i == len(score_matrix) - 1:
                hor_score = score_matrix[i][j - 1] - end_gap_penalty
            else:
                hor_score = score_matrix[i][j - 1] - gap_penalty
            if j == len(score_matrix[0]) - 1:
                vert_score = score_matrix[i - 1][j] - end_gap_penalty
            else:
                vert_score = score_matrix[i - 1][j] - gap_penalty
            diag_score = score_matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])

            # decide which score to use
            score_matrix[i][j] = max(hor_score, vert_score, diag_score)

            # update arrows, prefers diagonal over horizontal over vertical
            if score_matrix[i][j] == diag_score:
                arrow_matrix[i][j] += 'D'
            if score_matrix[i][j] == hor_score:
                arrow_matrix[i][j] += 'H'
            if score_matrix[i][j] == vert_score:
                arrow_matrix[i][j] += 'V'

    return score_matrix, arrow_matrix


def traceback_align(score_matrix, arrow_matrix, seq1, seq2):
    """Produces the alignment from score_matrix and arrow_matrix

    :param score_matrix: list of lists as produced by make_align_matrix with alignment scores.
    :param arrow_matrix: list of lists as produced by make_align_matrix with source direction ('D', 'H' or 'V')
    :param seq1: string of amino acids.
    :param seq2: string of amino acids.

    :return: align_seq1 and align_seq2. The aligned sequences with gaps inserted
    """
    assert len(score_matrix) == len(arrow_matrix), ValueError("score_matrix and arrow_matrix not same height")
    assert len(score_matrix[0]) == len(arrow_matrix[0]), ValueError("score_matrix and arrow_matrix not same width")
    align_lst1 = []
    align_lst2 = []
    y_pos = len(score_matrix) - 1
    x_pos = len(score_matrix[0]) - 1

    # I first had 'and' here but it did not accept both conditions. With 'or' it works. Don't know exactly why.
    # Now I do. You want to run until both are false. 'and' is only true when both are true. when one condition is
    # false, it evaluates to false. 'or' needs both to be false to return false.
    while y_pos != 0 or x_pos != 0:
        if 'D' in arrow_matrix[y_pos][x_pos]:
            align_lst1 += seq1[y_pos - 1]
            align_lst2 += seq2[x_pos - 1]
            x_pos = x_pos - 1
            y_pos = y_pos - 1
        elif 'H' in arrow_matrix[y_pos][x_pos]:
            align_lst1 += '-'
            align_lst2 += seq2[x_pos - 1]
            x_pos = x_pos - 1
        elif 'V' in arrow_matrix[y_pos][x_pos]:
            align_lst1 += seq1[y_pos - 1]
            align_lst2 += '-'
            y_pos = y_pos - 1
        else:
            break
    align_seq1 = ''.join(align_lst1[::-1])
    align_seq2 = ''.join(align_lst2[::-1])

    return align_seq1, align_seq2


def clean_seq(seq):
    """Does a sanity check on amino acid sequences to see if it uses amino acid letters. Raises exception if wrong.

    :param seq: string of amino acids
    :return: uppercase string of amino acids
    """
    assert type(seq) == str, "seq should be strings"
    seq = seq.upper()

    allowed_chars = list(BLOSUM62_ORDER.keys())[:-1]  # Don't allow * in chars.
    for ch in set(seq):
        if ch not in allowed_chars:
            raise ValueError("Please use valid amino acid letters in seq1:\n"
                             "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X")
    return seq


def align_seqs(seq1, seq2, gap_penalty=4, end_gap_penalty=None):
    """Returns the alignment of 2 protein sequences.

    :param seq1: string of amino acids.
    :param seq2: string of amino acids.
    :param gap_penalty: int of gap penalty to use.
    :param end_gap_penalty: int of beginning or end gap penalty to use. If none given same as gap_penalty

    :return: 2 Strings in alignment and int of total score.
    """
    # Sanity checks
    assert type(gap_penalty) == int, "gap_penalty should be int"

    # Check if end_gap_penalty used
    if end_gap_penalty is None:
        end_gap_penalty = gap_penalty
    else:
        assert type(end_gap_penalty) == int, "end_gap_penalty should be int"

    # make sure seqs are correctly formatted
    seq1 = clean_seq(seq1)
    seq2 = clean_seq(seq2)
    # Input is good & ready!

    # Populate score and arrow matrices with alignments of seq1 and seq2
    score_matrix, arrow_matrix = make_align_matrices(seq1, seq2, gap_penalty, end_gap_penalty)

    # for debugging:
    # for line in score_matrix:
    #    print(line)

    # Find alignment path
    align1, align2 = traceback_align(score_matrix, arrow_matrix, seq1, seq2)
    total_score = score_matrix[-1][-1]

    return align1, align2, total_score


def perc_identity(align_seq1, align_seq2):
    """Calculates percentage identity between two aligned sequences (gaps inserted)

    :param align_seq1:
    :param align_seq2:
    :return: percentage of identity.
    """
    assert len(align_seq2) == len(align_seq1), ValueError("sequences in perc_identity should have same length")
    num_identical_residues = 0
    for i in range(len(align_seq1)):
        if align_seq1[i] == align_seq2[i]:
            num_identical_residues += 1
    perc_iden = num_identical_residues / len(align_seq1) * 100
    return perc_iden


def print_alignment(align_seq1, align_seq2, score, print_score_identity=True):
    """Fancy printing of aligned sequences

    :param align_seq1: str of aligned sequence #1
    :param align_seq2: str of aligned sequence #2
    :param score: int of alignment score
    :param print_score_identity: Boolean flag to print score and identity

    :return: nothing, this prints
    """
    assert len(align_seq2) == len(align_seq1), ValueError("sequences should have same length, did you align them?")
    bars = []
    for i in range(len(align_seq1)):
        if align_seq1[i] == align_seq2[i]:
            bars += '|'
        else:
            bars += ' '
    print(align_seq1)
    print(''.join(bars))
    print(align_seq2)
    if print_score_identity:
        identity_report = 'Alignment score: {}.  Percentage identity {:.2f}%\n'.format(score, perc_identity(align_seq1,
                                                                                                            align_seq2))
        print(identity_report)
    return None


def main():
    """Main function to showcase usage of functions"""
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"

    seq1_aligned, seq2_aligned, total_score = align_seqs(seq1, seq2, gap_penalty=4)

    print_alignment(seq1_aligned, seq2_aligned, total_score)


if __name__ == "__main__":
    main()
