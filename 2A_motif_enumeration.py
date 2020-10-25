# Implementing code to find the (k, d) motifs in a collection of strings, where k is the length of the motif and d is
# allowed Hamming distance between variants of the motif.
# Problem 2A in the BALA textbook/Rosalind

# Input is a text file with first line as k and d, separated by a space, and subsequent lines corresponding to dna
# 3 1
# ATTTGGC
# TGCCTTA
# CGGTATC
# GAAAATT

# Output are the (k, d) motifs, separated by a space
# ATA ATT GTT TTT

# Usage: python3 2A_motif_enumeration.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_hamming_distance(str1, str2):
    """
    str1, str2 as strings of equal length
    Returns an int as a count of the positional mismatches between str1 and str2
    """
    # return sum([pattern1[i] != pattern2[i] for i in range(0, len(pattern1))])
    str1 = list(str1)
    str2 = list(str2)
    assert len(str1) == len(str2)
    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist


def neighbours(pattern, d):
    """
    text as a string containing at least one of [A,C,G,T]
    d as an int >= 0
    Returns a set of strings that are allowed up to d mismatches to pattern
   """
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighbourhood = set()
    suffix_neighbours = neighbours(pattern[1:], d)
    for text in suffix_neighbours:
        if get_hamming_distance(pattern[1:], text) < d:
            for nucleotide in ["A", "C", "G", "T"]:
                neighbourhood.add(nucleotide + text)
        else:
            neighbourhood.add(pattern[0] + text)
    return neighbourhood


def motif_enumeration(dna, k, d):
    """
    dna as a list of strings
    d as an int >= 0
    k as an int >= 0
    Returns a list of the (k, d) motifs in dna
    """
    all_neighbourhoods = [[] for i in range(len(dna))]
    for i in range(len(dna)):
        for j in range(len(dna[i]) - k+1):
            kmer = dna[i][j:j+k]
            neighbourhood = neighbours(kmer, d)
            for pattern in neighbourhood:
                all_neighbourhoods[i].append(pattern)
    result = set(all_neighbourhoods[0])
    for s in all_neighbourhoods[1:]:
        result.intersection_update(s)
    return " ".join(result)


def main():
    """
    Read the input file and parse the arguments. Return the (k, d) motifs of input dna
    """
    argv = list(sys.argv)
    input_text = []
    input_param = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_param is None:
                input_param = line.replace('\n', '').split()
                input_k = int(input_param[0])
                input_d = int(input_param[1])
            else:
                input_text.append(line.replace('\n', ''))

    output = motif_enumeration(input_text, input_k, input_d)

    print(output)


if __name__ == "__main__":
    main()

