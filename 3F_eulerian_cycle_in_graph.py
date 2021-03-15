# Implementing code to find a Euler circuit (every edge visited once and start and end on the same node) in a strongly
# connected and balanced graph using Hierholzer's algorithm, which is O(E) (linear where E is the count of edges).
# As the graph is balanced, we can only get stuck on the starting node. Once stuck, the algo backtracks to a node in the
# current circuit with unused edges, adding this new circuit to the previous, until all edges have been visited
# https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
# https://www.youtube.com/watch?v=8MpoO2zA2l4

# Problem 3F in the BALA textbook/Rosalind

# Input is a text file where each line corresponds to the nodes/edges of the graph
# 0 -> 3
# 1 -> 0
# 2 -> 1,6
# 3 -> 2
# 4 -> 2
# 5 -> 4
# 6 -> 5,8
# 7 -> 9
# 8 -> 7
# 9 -> 6

# Output is a string containing the Euler path
# 6->8->7->9->6->5->4->2->1->0->3->2->6

# Usage: python3 3F_eulerian_cycle_in_graph.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from random import choice


def parse_graph(list_input):
    """
    list_input: list containing nodes/edges of form '2 -> 1,6'
    return an adjacency list/dictionary
    """

    adjacency_dict = dict()
    for line in list_input:
        split_input = line.split("->")
        in_node = split_input[0].strip()
        out_nodes = split_input[1].split(",")
        if in_node not in adjacency_dict:
            adjacency_dict[in_node] = []
        for node in out_nodes:
            adjacency_dict[in_node].append(node.strip())
    return adjacency_dict


def eulerian_cycle(graph):
    """
    graph: adjacency list/dictionary
    returns a string of the eulerian circuit
    assumes: graph is strongly connected and balanced
    """

    current_path = list(choice(list(graph.keys())))  # initialize a stack with a random node
    cycle = []  # list to store final cycle

    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining cycle
            cycle.append(current_path.pop())  # Remove the current node and put it in the cycle

    cycle = "->".join(cycle[::-1])  # print in reverse
    print(cycle)


def main():
    """
    Read the input file, parse the arguments, and finds eulerian cycle
    """

    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    graph = parse_graph(input_text)
    eulerian_cycle(graph)


if __name__ == "__main__":
    main()
