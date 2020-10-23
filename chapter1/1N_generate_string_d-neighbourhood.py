# Implementing code to find the all strings that have a Hamming distance <= d to an input string
# Problem 1N in the BALA textbook/Rosalind

# Input is a text file with first line as the string, and the second line the d (max distance allowed)
# ACG
# 1

# Output are the neighbours, separated by new line
# ATGT ACAT

# Usage: python3 1N_generate_string_d-neighbourhood.py input.txt > output.txt
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


def main():
    """
    Read the input file and parse the arguments. Return the most frequent kmers with at most d mismatch in input text.
    """
    argv = list(sys.argv)
    input_text = None
    input_param = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_text is None:
                input_text = line.replace('\n', '')
            elif input_param is None:
                input_param = line.replace('\n', '').split()
                input_d = int(input_param[0])

    output = neighbours(input_text, input_d)

    print("\n".join(output))


if __name__ == "__main__":
    main()
