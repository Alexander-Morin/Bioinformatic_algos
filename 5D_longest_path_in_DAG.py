# Finding the longest path/highest score in a weighted DAG. Implement a dynamic programming approach that finds the
# max path score following a topological ordering, along with a backtracking approach to actually construct the path
# from source to sink.

# Problem 5D in the BALA textbook/Rosalind

# Input is a text file the first line is the graph source, the second line is the graph sink, and all following lines
# are the graph edges in form node_a->node_b:edge_weight
# 0
# 4
# 0->1:7
# 0->2:4
# 2->3:2
# 1->4:1
# 3->4:3

# Output is the score of the longest path and the corresponding path
# (9, '0->2->3->4')

# Usage: python3 5D_lcs_backtrack.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from collections import defaultdict
import copy
import pandas as pd
import math


def parse_weighted_graph(line_input):
    """
    list_input: list containing nodes/edges of form '2->1:6'
    return an adjacency list/dictionary
    """
    graph = defaultdict(list)
    for line in line_input:
        edge = line.split("->")
        node_a = int(edge[0])
        target = [int(x) for x in edge[1].split(":")]
        graph[node_a].append(target)
        graph[target[0]]  # init in case it only has edges coming in
    return graph


def count_weighted_edges(graph):
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
            target_node = target[0]  # just refer to node (first element) for key in count df
            edge_df.at[target_node, 'Edges_in'] += 1
    return edge_df


def topological_ordering_weighted_graph(graph):
    """
    graph: adjacency list (as a dict) of the graph
    candidates are nodes with no incoming edges. the algo goes through each candidate, removing the candidate and
    its outgoing edges, and then checking the next node with incoming edges in the resulting graph.
    returns a list of the ordered nodes
    """
    order = []
    count = count_weighted_edges(graph)
    candidates = list(count[count["Edges_in"] == 0].index.values)
    graph_cp = copy.deepcopy(graph)
    while candidates:
        node_a = candidates[0]
        order.append(node_a)
        candidates.remove(node_a)
        for node_b in graph[node_a]:
            graph_cp[node_a].remove(node_b)
            count.at[node_b[0], "Edges_in"] -= 1
            if count.at[node_b[0], "Edges_in"] == 0:
                candidates.append(node_b[0])
    if all(count["Edges_in"]) != 0:
        return "Input graph is not a DAG"
    return order


def recursive_backtrack(backtrack, source, node, path=[]):
    """
    backtrack: array of nodes used to reach the index-matched current node
    source: integer of the source node of the graph
    node: integer of the current node to start backtracking
    path: list of the traversed node path
    returns a character of the path travelled from node to source of form '0->2->3->4'
    """
    path.append(node)
    node = backtrack[node]
    if node == source:
        path.append(source)
        path.reverse()
        path = "->".join([str(x) for x in path])
        return path
    else:
        return recursive_backtrack(backtrack, source, node, path)


def longest_path(graph, source, sink):
    """
    graph: adjacency list (as a dict) of the graph
    source: integer of the source node of the graph
    sink: integer of the destination node of the graph
    returns the score and corresponding longest path from source to sink of form (9, '0->2->3->4')
    """
    order = topological_ordering_weighted_graph(graph)
    scores = [-math.inf] * (sink+1)
    scores[source] = 0
    backtrack = [None] * len(scores)
    for node in order:
        for target in graph[node]:
            if len(target) > 0:
                if scores[node] + target[1] > scores[target[0]]:
                    scores[target[0]] = scores[node] + target[1]
                    backtrack[target[0]] = node
    return scores[sink], recursive_backtrack(backtrack, source, sink)


def main():
    """
    Read the input file, parse the arguments, and returns the score and corresponding longest path
    """
    argv = list(sys.argv)
    input_source = None
    input_sink = None
    input_graph_list = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_source is None:
                input_source = int(line.replace('\n', ''))
            elif input_sink is None:
                input_sink = int(line.replace('\n', ''))
            else:
                input_graph_list.append(line.replace('\n', ''))

    input_graph = parse_weighted_graph(input_graph_list)
    print(longest_path(input_graph, input_source, input_sink))


if __name__ == "__main__":
    main()
