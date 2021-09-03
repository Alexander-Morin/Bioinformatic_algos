# Implementing topological ordering, where an input DAG's nodes are ordered in such a manner that a node is only visited
# after its predecessors have been visited. A DAG has at least one node with no incoming edges, which is used to seed
# the traversal.

# Problem 5N in the BALA textbook/Rosalind

# Input is a text file where each line represents edges of a DAG
# 1 -> 2
# 2 -> 3
# 4 -> 2
# 5 -> 3

# Output is the ordered nodes
# 1, 4, 5, 2, 3

# Usage: python3 5N_topological_ordering.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import pandas as pd
from collections import defaultdict
import copy


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
            graph[int(in_node)].append(int(node))
            graph[int(node)]  # initialize a node in case it only has edges coming in
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


def topological_ordering(graph):
    """
    graph: adjacency list (as a dict) of the graph
    candidates are nodes with no incoming edges. the algo goes through each candidate, removing the candidate and
    its outgoing edges, and then checking the next node with incoming edges in the resulting graph.
    returns a list of the ordered nodes
    """
    order = []
    count = count_edges(graph)
    candidates = list(count[count["Edges_in"] == 0].index.values)
    graph_cp = copy.deepcopy(graph)
    while candidates:
        node_a = candidates[0]
        order.append(node_a)
        candidates.remove(node_a)
        for node_b in graph[node_a]:
            graph_cp[node_a].remove(node_b)
            count.at[node_b, "Edges_in"] -= 1
            if count.at[node_b, "Edges_in"] == 0:
                candidates.append(node_b)
    if all(count["Edges_in"]) != 0:
        return "Input graph is not a DAG"
    return order


def main():
    """
    Read the input file, parse the arguments, and find the topologically ordered nodes
    """
    argv = list(sys.argv)
    input_graph_list = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_graph_list.append(line.replace('\n', ''))

    input_graph = parse_graph(input_graph_list)
    output = topological_ordering(input_graph)
    print(output)


if __name__ == "__main__":
    main()
