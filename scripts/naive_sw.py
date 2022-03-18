"""Naive Python implementation of the Smith-Waterman algorithm with affine gap scoring."""

import numpy as np
import argparse
from Bio import SeqIO

def sw(v, w, scoring_matrix, aa_dict, gap_opening, gap_extension):
    n = len(v)
    m = len(w)
    s_middle = np.zeros((n+1, m+1), dtype=float)
    b_middle = np.full((n+1, m+1), 3, dtype=int)
    s_lower = np.zeros((n+1, m+1), dtype=float)
    b_lower = np.full((n+1, m+1), 3, dtype=int)
    s_upper = np.zeros((n+1, m+1), dtype=float)
    b_upper = np.full((n+1, m+1), 3, dtype=int)

    # Begin filling:
    for i in range(1, n+1):
        for j in range(1, m+1):
            # lower
            if i > 1:
                b_lower[i, j], s_lower[i, j] = max(zip([3, 2, 0, 1],[0, s_middle[i-1, j] - gap_opening, s_lower[i-1, j] - gap_extension, float("-inf")]), key=lambda x: x[1])
            # upper
            if j > 1:
                b_upper[i, j], s_upper[i, j] = max(zip([3, 2, 0, 1],[0, s_middle[i, j-1] - gap_opening, float("-inf"), s_upper[i, j-1] - gap_extension]), key=lambda x: x[1])
            # middle
            match_score = scoring_matrix[aa_dict[v[i-1]], aa_dict[w[j-1]]] + s_middle[i-1, j-1]
            b_middle[i, j], s_middle[i, j] = max(zip([3, 2, 0, 1],[0, match_score, s_lower[i, j], s_upper[i, j]]), key=lambda x: x[1])
    
    # get (one) local alignment with maximal score
    max_score = 0
    sink = ()
    for i in range(n+1):
        for j in range(m+1):
            if s_middle[i, j] > max_score:
                max_score = s_middle[i, j]
                sink = (i, j)

    return s_middle[sink[0], sink[1]], b_middle, b_lower, b_upper, sink


def OutputLCS_AffineGaps(b_middle, b_lower, b_upper, v, w, i, j):
    al_v = ""
    al_w = ""
    matrix = "b_middle"
    while i + j != 0:
        if matrix == "b_middle":
            direction = b_middle[i, j]
            if direction == 0:
                matrix = "b_lower"
            elif direction == 1:
                matrix = "b_upper"
            elif direction == 2:
                al_v += v[i-1]
                al_w += w[j-1]
                i -= 1
                j -= 1
            else:
                break
        elif matrix == "b_lower":
            direction = b_lower[i, j]
            al_v += v[i-1]
            al_w += "-"
            i -= 1
            if direction == 2:
                matrix = "b_middle"
            elif direction == 3:
                break
        elif matrix == "b_upper":
            direction = b_upper[i, j]
            al_v += "-"
            al_w += w[j-1]
            j -= 1
            if direction == 2:
                matrix = "b_middle"
            elif direction == 3:
                break
    return al_v[::-1], al_w[::-1]


def substitution_matrix_parser(matrix_name):
    with open(f"../data/matrices/{matrix_name}") as f:
        matrix_string = ""
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()
        aas = line.split()
        while line:
            matrix_string += line.replace("\n", " ")
            line = f.readline()
    aa_dict = {}
    matrix_scores = []
    for el in matrix_string.split():
        try:
            matrix_scores.append(int(el))
        except ValueError:
            continue
    i = 0
    for aa in aas:
        aa_dict[aa] = i
        i += 1

    matrix_out = np.array(matrix_scores, dtype=int).reshape(len(aas), len(aas))

    return matrix_out, aa_dict


def main():
    """ 
    Smith-Waterman algorithm with affine gap scoring.
    """
    parser = argparse.ArgumentParser(description = 'Smith-Waterman algorithm with affine gap scoring.')

    parser.add_argument("query", 
        help="Query sequence(s)")
    parser.add_argument("db", 
        help="db sequence(s))")
    parser.add_argument("-m", default="EBLOSUM62",
        help="Substitution matrix (default: EBLOSUM62)")
    parser.add_argument("-o", type=float, default=10.0,
        help="Gap opening penalty")
    parser.add_argument("-e", type=float, default=0.5,
        help="Gap extension penalty. Penalty for a gap of n positions: gap opening penalty + (n - 1) * gap extension penalty")
    args = parser.parse_args()

    # parse substitution matrix
    substitution_matrix, aa_dict = substitution_matrix_parser(args.m)
    # loop over all query/db pairs
    for query in SeqIO.parse(args.query, "fasta"):
        #SeqIO.write(seq_record, f"data/{input_file.split('/')[-1]}_separate/{seq_record.name.split('|')[-1]}.fasta", "fasta")
        for db in SeqIO.parse(args.db, "fasta"):
            score, b_middle, b_lower, b_upper, sink = sw(query.seq, db.seq, substitution_matrix, aa_dict, args.o, args.e)
            al1, al2 = OutputLCS_AffineGaps(b_middle, b_lower, b_upper, query.seq, db.seq, sink[0], sink[1])
            print(f"{query.name}\t{db.name}\t{score}")
            print(al1)
            print(al2)
            print()
 
if __name__ == '__main__':
    main()
