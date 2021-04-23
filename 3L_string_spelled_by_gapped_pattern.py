# Implementing code to spell the path obtained from a paired De Bruijn graph, where the edges are (k,d)mers of form
# (kmer1|kmer2). Not every Eulerian path from a paired de bruijn graph will result in a fully reconstructed string,
# so this process can both filter a set of paths and spell the actual full string.

# Problem 3L in the BALA textbook/Rosalind

# Input is a text file where the first line is k (size of kmer) and d (gap) and following lines are the reads/kmers
# 4 2
# GACC|GCGC
# ACCG|CGCC
# CCGA|GCCG
# CGAG|CCGG
# GAGC|CGGA

# Output is a string containing the Euler path
# GACCGAGCGCCGGA

# Usage: python3 string_spelled_by_gapped_pattern.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def spell_string_by_path(text):
    """
    text: list of DNA strings assumed to be overlapping by len(string)-1
    returns the complete glued string
    """
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


def spell_string_by_gapped_pattern(patterns, k, d):
    """
    patterns: a list of paired read kmers of form "Kmer1|Kmer2"
    k: the size of the kmers
    d: the size of the gap between the paired read
    returns the path spelled by the paired edges
    """
    first_patterns = [kmer.split("|")[0] for kmer in patterns]
    second_patterns = [kmer.split("|")[1] for kmer in patterns]
    prefix_string = spell_string_by_path(first_patterns)
    suffix_string = spell_string_by_path(second_patterns)
    # imagine lining up rows of the reads, where each column position must match
    for i in range(k+d+1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i-k-d]:
            return "There is no consensus pattern"
    return prefix_string + suffix_string[len(suffix_string) - (k+d):]


def main():
    """
    Read the input file, parse the arguments, and returns the reconstructed string
    """
    argv = list(sys.argv)
    input_param = None
    input_k = None
    input_d = None
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_param is None:
                input_param = line.replace('\n', '').split()
                input_k = int(input_param[0])
                input_d = int(input_param[1])
            else:
                input_text.append(line.replace('\n', ''))

    print(spell_string_by_gapped_pattern(input_text, input_k, input_d))


if __name__ == "__main__":
    main()
