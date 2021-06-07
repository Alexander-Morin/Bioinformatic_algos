# Implementing code to find the all of the maximal non-branching paths of a graph - the longest paths whose intermediate
# nodes all have in/out degree equal to 1. This is used to assemble contigs in a genome, which are typically required
# as repeats in genomes make it impossible to find a unique Eulerian path, even with perfect coverage. The strings
# spelled by these contigs/maximal non branching paths are useful as they will be present in any given assembly with a
# given kmer composition. Unlike 3K, which finds all contigs in a given set of nucleotide kmers, this implementation
# also finds the isolated cycles of a graph (eg, 6 -> 7 -> 6) and expects a set of numeric graph edges as input.

# Problem 3M in the BALA textbook/Rosalind

# Input is a text file where each line corresponds to a set of edges in an adjacency list for a graph
# 1 -> 2
# 2 -> 3
# 3 -> 4,5
# 6 -> 7
# 7 -> 6

# Output is the set of maximum non branching paths in the graph
# 1 -> 2 -> 3
# 3 -> 4
# 3 -> 5
# 7 -> 6 -> 7

# Usage: python3 3M_all_max_nonbranching_paths.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import pandas as pd
from collections import defaultdict


def parse_graph(list_input):
    """
    list_input: list containing nodes/edges of form '2 -> 1,6'
    return an adjacency list/dictionary
    """
    graph = defaultdict(list)
    for line in list_input:
        split_input = line.split(" -> ")
        in_node = split_input[0].strip()
        out_nodes = split_input[1].split(",")
        for node in out_nodes:
            graph[in_node].append(node)
            graph[node]  # initialize a node in case it only has edges coming in
    return graph


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


def get_isolated_cycle(graph):
    """
    graph: adjacency list (as a dict) of the graph
    returns a list of the isolated cycles in the graph - eg, [7 -> 6 -> 7]
    """
    all_cycles = []
    count_df = count_edges(graph)
    count_df["Visited"] = False  # init column tracking if node has been visited
    for start_node in graph.keys():
        if count_df.at[start_node, "Visited"]:  # skip if seen to prevent adding both n1->n2->n1 and n2->n1->n1
            continue
        count_df.at[start_node, "Visited"] = True
        if is_1in1out(count_df, start_node):
            next_node = graph[start_node][0]
            cycle = [start_node, next_node]
            while is_1in1out(count_df, next_node).all():
                count_df.at[next_node, "Visited"] = True
                next_node = graph[next_node][0]
                cycle.append(next_node)
                if next_node is start_node:
                    all_cycles.append(cycle)
                    break
    return all_cycles


def maximal_nonbranching_paths(graph):
    """
    graph: adjacency list (as a dict) of the graph
    returns a list of all of the longest paths in graph that are non branching (all intermediate nodes have in and out
    degree equal to one)
    """
    paths = []
    count_df = count_edges(graph)
    for start_node in graph.keys():
        if not is_1in1out(count_df, start_node):  # iterate and start from all non 1-in-1-out nodes
            if count_df.loc[start_node, "Edges_out"] > 0:
                for next_node in graph[start_node]:
                    nb_path = [start_node, next_node]
                    while is_1in1out(count_df, next_node).all():  # all intermediate nodes are 1-in-1-out
                        next_node = graph[next_node][0]
                        nb_path.append(next_node)
                    paths.append(nb_path)
    for cycle in get_isolated_cycle(graph):  # need to include isolated cycles, eg 6 -> 7 -> 6
        paths.append(cycle)
    return paths


def main():
    """
    Read the input file, parse the arguments, and prints all maximal non branching paths
    """
    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    graph = parse_graph(input_text)
    paths = maximal_nonbranching_paths(graph)
    for path in paths:
        print(" -> ".join(path))


if __name__ == "__main__":
    main()
