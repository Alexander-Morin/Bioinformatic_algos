# Implementing code to find all of the contigs from patterns (a collection of kmers). Contigs are contiguous stretches
# of DNA in the genome, and are used as a compromise in assembly as repeats in genomes make it impossible to find unique
# Eulerian paths. Contigs are represented by non-branching paths in the De Bruijn graph of patterns - paths whose
# intermediate nodes have in and out degrees equal to one. Look to find all max length non-branching paths.

# Problem 3K in the BALA textbook/Rosalind

# Input is a text file where each line corresponds to a kmer in patterns
# ATG
# ATG
# TGT
# TGG
# CAT
# GGA
# GAT
# AGA

# Output is all contigs in the De Bruijn graph of patterns
# AGA ATG ATG CAT GAT TGGA TGT

# Usage: python3 3K_generate_all_contigs.py input.txt > output.txt
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


def is_1in1out(count_df, node):
    """
    count_df: pandas DF that has the count of edges in and edges out
    node: an index node of count_df
    returns True/False whether the input node has edges in and edges out equal to one
    """
    return (count_df.loc[node, "Edges_in"] == 1) and (count_df.loc[node, "Edges_out"] == 1)


def get_contigs(graph):
    """
    graph: adjacency list (as a dict) of the graph
    returns a list of all the contigs in the graph
    """
    contigs = []
    count_df = count_edges(graph)
    for start_node in graph.keys():
        if not is_1in1out(count_df, start_node):  # iterate and start from all non 1-in-1-out nodes
            if count_df.loc[start_node, "Edges_out"] > 0:
                for next_node in graph[start_node]:
                    nb_path = [start_node, next_node]
                    while is_1in1out(count_df, next_node).all():  # all intermediate nodes are 1-in-1-out
                        next_node = graph[next_node][0]
                        nb_path.append(next_node)
                    contigs.append(nb_path)
    contigs = [spell_string_by_path(i) for i in contigs]  # glue series of nodes into paths
    return contigs


def main():
    """
    Read the input file, parse the arguments, and prints all maximal non branching paths/contigs
    """
    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    graph = debruin_graph(input_text)
    contigs = get_contigs(graph)
    for contig in contigs:
        print(contig)


if __name__ == "__main__":
    main()
