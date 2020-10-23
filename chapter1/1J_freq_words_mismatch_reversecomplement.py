# Implementing code to find the most frequent kmer in a string of DNA nucleotides, allowing for up to d mismatches.
# Problem 1J in the BALA textbook/Rosalind

# Input is a text file formatted as string, k then d over two lines:
# ACGTTGCATGTCGCATGATGCATGAGAGCT
# 4 1

# Output is kmers separated by a space
# ATGT ACAT

# Usage: IJ_freq_words_mismatch_reversecomplement.py input.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_reverse_complement(text):
    """
    text as a string containing at least one of [A,C,G,T]
    Returns a string of the reverse complement of text (ACTG -> CAGT)
    """
    complement = []
    for letter in text:
        if letter == 'A':
            complement.append('T')
        elif letter == 'C':
            complement.append('G')
        elif letter == 'T':
            complement.append('A')
        elif letter == 'G':
            complement.append('C')
    return "".join(complement[::-1])


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


def freq_kmers_mismatch_sort_rev(text, k, d):
    """
    text as a string containing at least one of [A,C,G,T]
    k as an int >= 2
    d as an int >= 0
    Returns a list of strings that represent the most frequent kmers in text, allowing up to d mismatches and c
    considering the forward and reverse direction
    """
    freq_patterns = []
    neighbourhood = {}

    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        rev_pattern = get_reverse_complement(pattern)
        for j in neighbours(pattern, d):
            if j in neighbourhood:
                neighbourhood[j] += 1
            else:
                neighbourhood[j] = 1
        for j in neighbours(rev_pattern, d):
            if j in neighbourhood:
                neighbourhood[j] += 1
            else:
                neighbourhood[j] = 1

    max_count = max(neighbourhood.values())

    for i in neighbourhood:
        if neighbourhood[i] == max_count:
            freq_patterns.append(i)

    return freq_patterns


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
                input_k = int(input_param[0])
                input_d = int(input_param[1])

    output = freq_kmers_mismatch_sort_rev(input_text, input_k, input_d)

    print(output)


if __name__ == "__main__":
    main()
