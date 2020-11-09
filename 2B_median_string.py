# Implementing code to find the median string - the kmer that minimizes the sum of Hamming distances between a
# kmer/pattern and the possible motifs within dna, a collection of strings.
# Problem 2B in the BALA textbook/Rosalind

# Input is a text file with first line as k, and subsequent lines corresponding to dna
# 3
# AAATTGACGCAT
# GACGACCACGTT
# CGTCAGCGCCTG
# GCTGAGCACCGG
# AGTACGGGACAG

# Output is the median string
# GAC

# Usage: python3 2B_median_string.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import math


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


def number_to_symbol(number):
    """
    number as an integer in [0, 1, 2, 3]
    returns one of "A", "C", "G", "T"
    """
    assert number in [0, 1, 2, 3], "Invalid index"
    if number == 0:
        return "A"
    elif number == 1:
        return "C"
    elif number == 2:
        return "G"
    elif number == 3:
        return "T"


def number_to_pattern(index, k):
    """
       index as an integer between 0 and (4 ** k) - 1
       k as an integer (practically >= 3, <= 12)
       returns a pattern of length k corresponding to the lexicographically sorted kmers of length k
       """
    if k == 1:
        return number_to_symbol(index)
    prefix_index = index // 4
    remainder = index % 4
    symbol = number_to_symbol(remainder)
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    return prefix_pattern + symbol


def get_pattern_distance_from_dna(pattern, dna):
    """
       pattern as a string of dna nucleotides
       dna as a list of equal length strings
       returns an integer corresponding to the sum of the minimized distances between pattern and its closest motif
       match in each string of dna
       """
    k = len(pattern)
    distance = 0
    for string in dna:
        hdist = math.inf
        for i in range(len(string) - k + 1):
            kmer = string[i:i+k]
            if hdist > get_hamming_distance(pattern, kmer):
                hdist = get_hamming_distance(pattern, kmer)
        distance += hdist
    return distance


def median_string(dna, k):
    """
    dna as a list of equal length strings
    k as an integer (practically >= 3, <= 12)
    returns a string that represents the motif that minimizes the distance of d(pattern, dna)
    """
    distance = math.inf
    for i in range(4**k):
        pattern = number_to_pattern(i, k)
        pattern_dist = get_pattern_distance_from_dna(pattern, dna)
        if distance > pattern_dist:
            distance = pattern_dist
            median = pattern
    return median


def main():
    """
    Read the input file and parse the arguments. Return the (k, d) motifs of input dna
    """
    argv = list(sys.argv)
    input_text = []
    input_k = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_k is None:
                input_k = int(line.replace('\n', ''))
            else:
                input_text.append(line.replace('\n', ''))

    output = median_string(input_text, input_k)

    print(output)


if __name__ == "__main__":
    main()
