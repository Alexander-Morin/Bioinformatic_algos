# Implementing code to print out the adjacency list/dict of an overlap graph
# Problem 3C in the BALA textbook/Rosalind

# Input is a text file with each line as a pattern in text
# ATGCG
# GCATG
# CATGC
# AGGCA
# GGCAT

# Output is the printed adjacency list
# AGGCA -> GGCAT
# CATGC -> ATGCG
# GCATG -> CATGC
# GGCAT -> GCATG

# Usage: python3 3C_construct_overlap_graph.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_adjacency_dict(patterns):
    """
    patterns: list of DNA strings
    returns an adjacency dictionary corresponding to the overlap graph
    """
    adjacency_dict = dict.fromkeys(patterns)
    for node1 in adjacency_dict:
        suffix = node1[1:]
        for node2 in patterns:
            prefix = node2[:-1]
            if suffix == prefix:
                adjacency_dict[node1] = node2
    return adjacency_dict


def main():
    """
    Read the input file, parse the arguments, build the adjacency list/dict, then print it out line by line
    """
    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    output = get_adjacency_dict(input_text)

    for node1 in output.keys():
        if output[node1] is not None:
            print(node1 + " -> " + output[node1])


if __name__ == "__main__":
    main()
