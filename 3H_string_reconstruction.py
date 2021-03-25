# Implementing code to solve the string reconstruction problem by finding a Eulerian path in the De Bruijn graph of
# sequenced reads.

# Problem 3H in the BALA textbook/Rosalind

# Input is a text file where the first line is k (size of kmer) and following lines are the reads/kmers
# 4
# CTTA
# ACCA
# TACC
# GGCT
# GCTT
# TTAC

# Output is a string containing the Euler path
# GGCTTACCA

# Usage: python3 3H_string_reconstruction.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import pandas as pd
from collections import defaultdict


def spell_string_by_path(text):
    """
    text: list of DNA strings assumed to be overlapping by len(string)-1
    returns the complete glued string
    """
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


def debruin_graph(patterns):
    """
    patterns: a list of equal length strings (kmers/reads)
    return the de bruijn graph as an adjacency list/dict
    """
    adjacency_dict = defaultdict(list)
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_dict[prefix].append(suffix)
        adjacency_dict[suffix]  # init suffix in case it has no outgoing edges
    return adjacency_dict


def count_edges(graph):
    """
    graph: adjacency list (as a dict) of the graph
    returns a pandas DF with nodes as index and columns for count of edges in and out
    """
    nodes = set(graph)
    edge_df = pd.DataFrame(
        {"Edges_in": [0] * len(nodes),
         "Edges_out": [0] * len(nodes)},
        index=nodes)
    for node in graph:
        edge_df.at[node, 'Edges_out'] = len(graph[node])
        for target in graph[node]:
            edge_df.at[target, 'Edges_in'] += 1
    return edge_df


def get_start_node(graph):
    """
    graph: adjacency list (as a dict) of the graph
    returns the start node as a str
    """
    count_df = count_edges(graph)
    start = None
    for node in count_df.index:
        if count_df.at[node, "Edges_out"] - count_df.at[node, "Edges_in"] == 1:
            start = node
    if start is None:
        raise ValueError("No start node was found")
    return start


def eulerian_path(graph):
    """
    graph: adjacency list/dictionary
    returns a string of the eulerian path
    """
    current_path = [get_start_node(graph)]  # initialize a stack with the starting node
    path = []  # list to store final path
    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining cycle
            path.append(current_path.pop())  # Remove the current node and put it in the cycle

    path = "->".join(path[::-1])  # print in reverse
    return path


def string_reconstruction(kmer_list):
    """
    kmer_list: list of equal length kmers
    returns a string that uses each kmer exactly once, with kmers overlapping by len(kmer) - 1
    """
    graph = debruin_graph(kmer_list)
    path = eulerian_path(graph).split("->")
    return spell_string_by_path(path)


def main():
    """
    Read the input file, parse the arguments, and returns the reconstructed string
    """
    argv = list(sys.argv)
    input_k = None  # don't actually need to use this!
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_k is None:
                input_k = int(line.replace('\n', ''))
            else:
                input_text.append(line.replace('\n', ''))

    print(string_reconstruction(input_text))


if __name__ == "__main__":
    main()
