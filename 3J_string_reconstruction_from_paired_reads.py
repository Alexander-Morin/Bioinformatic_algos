# Implementing code to reconstruct a string from its paired composition.

# Problem 3J in the BALA textbook/Rosalind

# Input is a text file where the first line is k (size of kmer) and d (gap) and following lines are the reads
# 4 2
# GAGA|TTGA
# TCGT|GATG
# CGTG|ATGT
# TGGT|TGAG
# GTGA|TGTT
# GTGG|GTGA
# TGAG|GTTG
# GGTC|GAGA
# GTCG|AGAT

# Output is a string containing the Euler path
# GTGGTCGTGAGATGTTGA

# Usage: python3 3J_string_reconstruction_from_paired_reads.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import pandas as pd
from random import choice
from collections import defaultdict, deque


def split_read_patterns(paired_list):
    """
    paired_list: list containing paired reads of form "Kmer1|Kmer2"
    returns two lists containing all reads from the respective pairings
    """
    read1_patterns = []
    read2_patterns = []
    for pattern in paired_list:
        read_split = pattern.split("|")
        read1_patterns.append(read_split[0])
        read2_patterns.append(read_split[1])
    return read1_patterns, read2_patterns


def paired_dbg_from_composition(paired_list):
    """
    paired_list: list containing paired reads of form "Kmer1|Kmer2"
    returns an adjacency list/dict of the De Bruijn graph of the paired reads
    """
    read1_patterns, read2_patterns = split_read_patterns(paired_list)
    graph = defaultdict(list)
    for i in range(len(read1_patterns)):
        read_pair_prefix = "|".join([read1_patterns[i][:-1], read2_patterns[i][:-1]])
        read_pair_suffix = "|".join([read1_patterns[i][1:], read2_patterns[i][1:]])
        graph[read_pair_prefix].append(read_pair_suffix)
        graph[read_pair_suffix]  # init suffix in case it has no outgoing edges
    return graph


def count_edges(graph):
    """
    graph: adjacency list/dict of the graph
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
    Used for finding a Eulerian path, where the first node in the path has node degree out - in == 1
    graph: adjacency list/dict of the graph
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
    Spells out the path that visits every edge of the graph exactly once
    graph: adjacency list/dict of the graph
    returns a string of the eulerian path
    """
    current_path = [get_start_node(graph)]  # initialize a stack with the starting node
    path = []  # list to store final path
    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining path
            path.append(current_path.pop())  # Remove the current node and put it in the cycle
    path = "->".join(path[::-1])  # print in reverse
    return path


def bfs(graph):
    """
    Breadth first search to see if a graph is connected
    graph: adjacency list/dict of the graph
    Returns a logical vector of whether the nodes in graph are able to be visited or not
    """
    s = choice(list(graph.keys()))  # start with random node
    visited = [False] * len(graph.keys())  # array keeps track of what nodes we have already seen
    queue = deque([s])  # use a queue to track nodes that have yet to be visited, beginning with start
    visited[int(s)] = True
    while queue:
        s = queue.popleft()
        for i in graph[s]:
            if not visited[int(i)]:
                queue.append(i)
                visited[int(i)] = True
    return visited


def get_bypass_graphs(graph):
    """
    If graph has nodes with in degrees greater than 1, return a list of simple graphs where these nodes are bypassed
    graph: adjacency list/dict of the graph
    returns a list of graphs as adjacency lists/dicts
    """
    all_graphs = [graph]
    for g in all_graphs:
        edges_in = count_edges(g)["Edges_in"]
        if all(edges_in < 2):
            continue
        in_gt1 = choice(edges_in[edges_in > 1].index.values)
        edges_in = [edge for edge in g if in_gt1 in g[edge]]
        edges_out = g[in_gt1]
        for node1 in edges_in:
            for node2 in edges_out:
                g_cp = g.copy()
                g_cp[node1] = [node2]
                if all(bfs(g_cp)):  # the more simple graph must be connected
                    all_graphs.append(g_cp)
        all_graphs.remove(g)
    return all_graphs


def spell_string_by_path(text):
    """
    text: list of DNA strings assumed to be overlapping by len(string)-1
    returns the complete glued string
    """
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


def spell_string_by_gapped_pattern(patterns, k, d):
    """
    patterns: a list of paired read kmers of form "Kmer1|Kmer2"
    # k: the size of the kmers
    # d: the size of the gap between the paired read
    returns the path spelled by the paired edges
    """
    first_patterns = [kmer.split("|")[0] for kmer in patterns]
    second_patterns = [kmer.split("|")[1] for kmer in patterns]
    prefix_string = spell_string_by_path(first_patterns)
    suffix_string = spell_string_by_path(second_patterns)
    # imagine lining up rows of the reads, where each column position must match
    for i in range(k+d+1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i-k-d]:
            return "There is no consensus pattern"
    return prefix_string + suffix_string[len(suffix_string) - (k+d):]


def spell_all_eulerian_paths(graph_list, k, d):
    """
    graph_list: list of graphs as adjacency lists/dicts
    # k: the size of the kmers
    # d: the size of the gap between the paired read
    returns
    """
    for g in graph_list:
        cycle = eulerian_path(g).split("->")
        path = spell_string_by_gapped_pattern(cycle, k, d)
        print(path)
    return None


def main():
    """
    Read the input file, parse the arguments, and returns the reconstructed string
    """
    argv = list(sys.argv)
    input_param = None
    input_k = None
    input_d = None
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_param is None:
                input_param = line.replace('\n', '').split()
                input_k = int(input_param[0])
                input_d = int(input_param[1])
            else:
                input_text.append(line.replace('\n', ''))

    graph = paired_dbg_from_composition(input_text)
    bpg = get_bypass_graphs(graph)
    spell_all_eulerian_paths(bpg, input_k, input_d)


if __name__ == "__main__":
    main()
