# Implementing code to find a Eulerian path in a graph: the path where every edge is visited exactly once, but unlike
# Eulerian circuit, the start and end nodes are different. Requires that exactly one node has out edges - in edges == 1,
# which is the start node, and exactly one node has in edges - out edges == 1, which is the final node. Once the start
# node is found, use the same Hierholzer's algo as the Eulerian circuit problem to find the path in O(E) time.
# https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
# https://www.youtube.com/watch?v=8MpoO2zA2l4

# Problem 3G in the BALA textbook/Rosalind

# Input is a text file where each line corresponds to the nodes/edges of the graph
# 0 -> 2
# 1 -> 3
# 2 -> 1
# 3 -> 0,4
# 6 -> 3,7
# 7 -> 8
# 8 -> 9
# 9 -> 6

# Output is a string containing the Euler path
# 6->7->8->9->6->3->0->2->1->3->4

# Usage: python3 3G_eulerian_path_in_graph.py input.txt > output.txt
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
    print(path)


def main():
    """
    Read the input file, parse the arguments, and finds eulerian path
    """
    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    graph = parse_graph(input_text)
    eulerian_path(graph)


if __name__ == "__main__":
    main()
