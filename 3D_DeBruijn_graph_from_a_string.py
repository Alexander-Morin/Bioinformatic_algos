# Implementing code to print out the adjacency list/dict of a De Bruijn graph (DGB) for kmers in input text. The DBG
# represents the kmers as edges, and so a node represents the k-1 overlap between two kmers. Identically labelled nodes
# are 'glued' together - so if ATG was repeated 3 times in a text, where would be 3 edges between AT and TG

# Problem 3D in the BALA textbook/Rosalind

# Input is a text file with first line as k, and the following line as the DNA string
# 4
# AAGATTCTCTAC

# Output is the printed adjacency list
# AAG -> AGA
# AGA -> GAT
# ATT -> TTC
# CTA -> TAC
# CTC -> TCT
# GAT -> ATT
# TCT -> CTA,CTC
# TTC -> TCT

# Usage: python3 3D_DeBruijn_graph_from_a_string.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def de_graph(text, k):
    """
    text: string of DNA characters
    k: length of the kmer (which is represented by the edges, so nodes are kmers of length k-1
    returns an adjacency list (as a dict) for the DeBruijn graph for kmers in text
    """
    adjacency_dict = dict()
    for i in range(len(text) - k+1):
        kmer = text[i:i+k-1]
        next_kmer = text[i+1:i+k]
        if kmer not in adjacency_dict:
            adjacency_dict[kmer] = [next_kmer]
        else:
            adjacency_dict[kmer].append(next_kmer)
    return adjacency_dict


def print_sorted_de_graph(de_dict):
    """
    de_dict: the adjacency dict for a de bruijn graph.
    prints the nodes as node1 -> node2,node3... in lexico. order
    """
    for node in sorted(de_dict):
        print(node + " -> " + ",".join(sorted(de_dict[node])))


def main():
    """
    Read the input file, parse the arguments, build the adjacency list/dict, then print it out line by line
    """

    argv = list(sys.argv)
    input_k = None
    input_text = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_k is None:
                input_k = int(line.replace('\n', ''))
            else:
                input_text = str(line.replace('\n', ''))

    output = de_graph(input_text, input_k)

    print_sorted_de_graph(output)


if __name__ == "__main__":
    main()
