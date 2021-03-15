# Chapter 3 is focused on genome assembly

from random import choice
# import numpy as np
from collections import defaultdict

# First, focus on the more simplistic case where are all reads are kmers of length k, all reads come from the same
# strand, there are no errors, and they exhibit perfect coverage (every kmer substring of the genome has a read).
# Consider the kmer composition of a string (lexicographically sorted kmers of a string) - to model genome assembly
# must perform the inverse problem, where kmers are assembled into a string. The chapter starts by motivating the
# (k-1) overlap of kmers to spell out the string. This task is complicated by repeats and constructing the string's
# path requires knowing the genome in advance
# ----------------------------------------------------------------------------------------------------------------------


def string_composition(text, k):
    kmers = [text[i:i+k] for i in range(0, len(text) - k+1)]
    kmers.sort()
    return kmers


def spell_string_by_path(text):
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


# The problem is then phrased in terms of graphs: Kmers are connected if the suffix of one is equal to the prefix of
# another, resulting in a directed overlap graph where the kmer nodes are sorted lexicographically. Spelling out the
# string means visiting each node exactly once. Introduces the overlap graph problem, with (k-1) overlap and returning
# an adjacency list (which is actually better represented as a mapping => dict)
# ----------------------------------------------------------------------------------------------------------------------


def get_adjacency_dict(patterns):
    adjacency_dict = dict.fromkeys(patterns)
    for node1 in adjacency_dict:
        suffix = node1[1:]
        for node2 in patterns:
            prefix = node2[:-1]
            if suffix == prefix:
                adjacency_dict[node1] = node2
    return adjacency_dict


# To solve the string reconstruction problem, looking for a Hamiltonian path: a path in the graph that visits every
# node once. Note that there may be more than one Hamiltonian path in a graph. States that an alternative way to find
# Hamiltonian paths in large graphs is to follow de Bruijn's approach: represent a kmer composition as a graph by
# treating the kmers as edges, rather than nodes. Pairs of consecutive edges represent kmers that overlap in (k-1)
# nucleotides, so the node with the incoming and outgoing edges are labeled as the (k-1)mers. This results in the
# merging of identical nodes into one - decreases count of nodes while preserving the count of edges. For example, if
# text has 3 repeats of ATG, there would be 3 edges connecting node AT to node TG
# ----------------------------------------------------------------------------------------------------------------------


def de_graph_from_string(text, k):
    adjacency_dict = dict()
    for i in range(len(text) - k+1):
        kmer = text[i:i+k-1]
        next_kmer = text[i+1:i+k]
        if kmer not in adjacency_dict:
            adjacency_dict[kmer] = [next_kmer]
        else:
            adjacency_dict[kmer].append(next_kmer)
    return adjacency_dict


# Solving the string reconstruction problem reduces to finding a path that visits every edge of the De Bruijn graph
# exactly once == Eulerian path. The text then talks about the composition graph of text, which has all of the kmers
# as isolated edges. Gluing nodes wit the same label in the composition graph produces DeBruijn(text). The text then
# suggests a way of constructing the De Bruijn graph that does involve gluing: given a collection of kmer patterns, the
# nodes of DeBruijn(patterns) are all unique (k-1)mers occurring as a prefix or suffix. DeBruijn(patterns) forms a
# directed graph by connecting prefix to suffix. However, I am not certain how this is different from constructing
# from a text and k...
# ----------------------------------------------------------------------------------------------------------------------


def de_graph_from_pattern(patterns):
    adjacency_dict = dict()
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in adjacency_dict:
            adjacency_dict[prefix] = [suffix]
        else:
            adjacency_dict[prefix].append(suffix)
    return adjacency_dict


# Anywho, states that we can now solve the string reconstruction problem by either finding a Hamiltonian path (ever node
# once) in the Overlap graph, or a Eulerian path (every edge once) in the DBG. Their provided visualization shows that,
# for small graphs, the DBG is more attractive as it is smaller. But the key is to find an efficient algorithm that
# solves large graphs. Hamiltonian path is NP complete. Euler's theorem: every balanced (every node has equal in
# degree to out degree) and strongly connected (possible to reach any node from any node) graph is Eulerian, meaning
# you can visit every edge once while starting/finishing on the same node. Because the graph is balanced, you can only
# get stuck on the starting node. The idea is to backtrack to a node in the current cycle that has unused edges and
# form a new circuit, adding this circuit to the previous until all edges have been visited. O(E) - linear over the
# number of edges in the graph.
# ----------------------------------------------------------------------------------------------------------------------


def parse_graph(list_input):
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



graph_input = [
    "0 -> 3",
    "1 -> 0",
    "2 -> 1,6",
    "3 -> 2",
    "4 -> 2",
    "5 -> 4",
    "6 -> 5,8",
    "7 -> 9",
    "8 -> 7",
    "9 -> 6"
]



graph = parse_graph(graph_input)
eulerian_cycle(graph)
