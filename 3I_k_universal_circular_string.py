# Implementing code to find a k-universal circular binary string: find a string of 0s and 1s that contains every binary
# kmer for an input k. This task reduces to finding a Eulerian cycle in the De Bruijn graph, and models the string
# reconstruction problem for a circular genome, like a bacterial chromosome.

# Problem 3I in the BALA textbook/Rosalind

# Input is a text file where the first/only line is an int
# 4

# Output
# 0000110010111101

# Usage: python3 3I_k_universal_circular_string.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from random import choice
from collections import defaultdict


def debruijn_graph(patterns):
    """
    patterns: list of binary patterns
    returns an adjacency list/dict of the de bruin graph of patterns
    """
    adjacency_dict = defaultdict(list)
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_dict[prefix].append(suffix)
        adjacency_dict[suffix]  # init suffix in case it has no outgoing edges
    return adjacency_dict


def eulerian_cycle(graph):
    """
    graph: adjacency list/dictionary
    returns a string of the eulerian circuit
    assumes: graph is strongly connected and balanced
    """

    current_path = [choice(list(graph.keys()))]  # initialize a stack with a random node
    cycle = []  # list to store final cycle

    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining cycle
            cycle.append(current_path.pop())  # Remove the current node and put it in the cycle

    cycle = "->".join(cycle[::-1])  # print in reverse
    return cycle


def binary_kmers(k):
    """
    k: an int, recommended 1 < k < 20
    returns a list of all binary kmers
    https://stackoverflow.com/questions/64890117/
    """
    def recursion(k, kmer="", kmer_list=[]):
        if len(kmer) == k:
            kmer_list.append(kmer)
        else:
            recursion(k, kmer + "0", kmer_list)
            recursion(k, kmer + "1", kmer_list)
        return kmer_list
    return recursion(k)


def spell_string_by_path(text):
    """
    text: list of (k-1)mers
    returns the complete glued string
    """
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


def k_circular_universal_string(k):
    """
    k: an int, recommended 1 < k < 20
    returns the k universal circular string
    """
    kmer_list = binary_kmers(k)
    dbg = debruijn_graph(kmer_list)
    cycle = eulerian_cycle(dbg).split("->")
    path = spell_string_by_path(cycle)
    path = path[:len(path) - int(k)+1]
    return path


def main():
    """
    Read the input file, parse the arguments, and prints the k universal circular string
    """
    argv = list(sys.argv)
    input_k = None

    for line in fileinput.input(argv[1]):
        if input_k is None:
            input_k = int(line)

    print(k_circular_universal_string(input_k))


if __name__ == "__main__":
    main()
