# This chapter explores how we compare DNA sequences using dynamic programming

import numpy as np
import pandas as pd
from collections import defaultdict
import copy
import math

# Text begins with the motivation of cracking the non-ribosomal code - non-ribosomal peptides (NRPs) produced by NRP
# synthetase. Specifically, looking at the adenylation domains (A-domains) of NRP synthetase responsible for adding
# a single amino acid. Researchers wanted to compare the sequences of these A domains, and began to appreciate that
# while there was a conserved core, it was obscured by spaces in the pairwise alignment of the sequences. The goal is
# to then align these sequences while being aware that insertions and deletions are possible. Finding the longest
# common subsequence can be phrased as the Manhattan tourist problem, wherein the longest path in a direct acyclic graph
# is found. This is accomplished using dynamic programming - the text motives using dynamic programming by looking at
# the minimum coin change problem. Start with greedy implementation (doesn't guarantee minimum), recursive solution
# (very inefficient b/c same values are re-calculated), then finally DP implementation
# ----------------------------------------------------------------------------------------------------------------------


def greedy_change(money, coins):
    change = []
    coins = np.array(coins)
    while money > 0:
        coin = max(coins[coins <= money])
        change.append(coin)
        money -= coin
    return change


def recursive_change(money, coins):
    if money == 0:
        return 0
    min_num_coins = np.inf
    for i in range(len(coins)):
        if money >= coins[i]:
            num_coins = recursive_change(money - coins[i], coins)
            if num_coins + 1 < min_num_coins:
                min_num_coins = num_coins + 1
    return min_num_coins


# O(money * len(coins))
def dp_change(money, coins):
    min_num_coins = [0] * (money + 1)
    for m in range(1, money + 1):
        min_num_coins[m] = np.inf
        for i in range(len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coins[i]] + 1
    return min_num_coins[money]


# The goal then becomes to implement a DP approach to the Manhattan tourist problem. Imagine an ixj table starting at
# top left (seed: 0,0) and you can only go right from (i, j-1) or down from (i-1, j). Each path has a weight, and want
# to get to the bottom right point (sink, n, m) via a path that has the highest sum weights. The first solution is
# recursive and re-calculates steps many times. The DP approach calculates every sub-step to find the ultimate max.
# ----------------------------------------------------------------------------------------------------------------------


# recursive
def south_or_east(i, j, v_edge_weight, h_edge_weight):
    if i == 0 and j == 0:
        return 0
    x, j = -np.inf
    if i > 0:
        x = south_or_east(i-1, j) + h_edge_weight(i, j)
    if j > 0:
        y = south_or_east(i, j-1) + v_edge_weight(i, j)
    return max(x, y)


def manhattan_tourist(i, j, v_edge_weight, h_edge_weight):
    scores = np.zeros([i+1, j+1])
    for x in range(1, i+1):
        scores[x, 0] = scores[x-1, 0] + v_edge_weight[x-1, 0]
    for y in range(1, j+1):
        scores[0, y] = scores[0, y-1] + h_edge_weight[0, y-1]
    for x in range(1, i+1):
        for y in range(1, j+1):
            scores[x, y] = max(scores[x-1, y] + v_edge_weight[x-1, y], scores[x, y-1] + h_edge_weight[x, y-1])
    return int(scores[i, j])


# From here, the task becomes to adapt the Manhattan tourist to alignment graphs that have diagonal edges, which
# represent matches. The text notes that we can find the longest common subsequence (LCS) without building the "city
# grid." However, "the arguments to do so are tedious." It also notes that many applications of alignment are much
# more complex than the LCS problem, and require building a DAG with specifically learned edge weights to reflect
# biology. Goal is to just then learn a generic DP approach that will find the longest path in any DAG, which itself
# is a tool that can be applied to many other problems. Note that before we were working with a rectangular graph with
# inherit order to the columns/rows. This ordering is important as DP requires previous problems to be solved/nodes
# have been visited. The goal of using DP approach for finding longest paths is to choose a path that maximizes the
# recurrence: sb = max of all predecessors nodes a of b { sa + weight of edge a to b }. Therefore must topologically
# order the graph such that every edge (ai, aj) of the DAG connects with a node with a smaller index to a larger index.
# This ordering can be done in time proportional to the number of edges in the graph. Works off of the fact that every
# DAG has at least one node with no incoming edges - remove this node and its outgoing edges. The following DAG will
# also have a node with no incoming edges. Proceed until all nodes removed.
# ----------------------------------------------------------------------------------------------------------------------


def parse_graph(list_input):
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


# Once we have a topological ordering, find the length of the longest path from source to sink by visiting the nodes
# in their ordered form. For simplicity assume that the only the source node is the only node with in degree == 0.
# Note modification to graph helpers to account for node_a->node_b:weight structure of the input - each edge now
# represented as a list where the first element is node a and the second element is the weight
# ----------------------------------------------------------------------------------------------------------------------


def parse_weighted_graph(line_input):
    graph = defaultdict(list)
    for line in line_input:
        edge = line.split("->")
        node_a = int(edge[0])
        target = [int(x) for x in edge[1].split(":")]
        graph[node_a].append(target)
        graph[target[0]]  # init in case it only has edges coming in
    return graph


def count_weighted_edges(graph):
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


# While we have the longest path, we still need to reconstruct the LCS. Note that tracing a path from source to sink
# may reach a dead end, but all paths starting from sink will lead to the source. Use a backtrack algorithm with
# pointers that keeps track of the path (down=0, right=0, diagonal=1) in a matrix representing the match between two
# strings.
# ----------------------------------------------------------------------------------------------------------------------


def lcs_backtrack(string1, string2):

    scores = np.zeros([len(string1)+1, len(string2)+1])
    backtrack = np.empty([len(string1)+1, len(string2)+1], dtype="object")

    for i in range(1, len(string1)+1):
        for j in range(1, len(string2)+1):

            if string1[i-1] == string2[j-1]:
                scores[i, j] = max(scores[i-1, j-1] + 1, scores[i-1, j], scores[i, j-1])
            else:
                scores[i, j] = max(scores[i-1, j], scores[i, j-1])

            if scores[i, j] == scores[i-1, j]:
                backtrack[i, j] = "Down"
            elif scores[i, j] == scores[i, j-1]:
                backtrack[i, j] = "Right"
            elif (scores[i, j] == scores[i-1, j-1] + 1) & (string1[i-1] == string2[j-1]):
                backtrack[i, j] = "Diagonal"

    return backtrack


def output_lcs(backtrack, string1, i, j):

    output = []

    def recursive_output(backtrack, string1, i, j, output):
        if i == 0 or j == 0:
            return None
        if backtrack[i, j] == "Down":
            recursive_output(backtrack, string1, i-1, j, output)
        elif backtrack[i, j] == "Right":
            recursive_output(backtrack, string1, i, j-1, output)
        else:
            recursive_output(backtrack, string1, i-1, j-1, output)
            output.append(string1[i-1])

    recursive_output(backtrack, string1, i, j, output)
    return "".join(output)


#


# Example inputs
# ----------------------------------------------------------------------------------------------------------------------

# print(greedy_change(48, [120, 40, 30, 24, 20, 10, 5, 4, 1]))
# print(recursive_change(40, [5, 4, 1]))
# print(dp_change(40, [1, 5, 10, 20, 25, 50]))

# i, j = 4, 4
# v_weights = np.array([[1, 0, 2, 4, 3],
#                       [4, 6, 5, 2, 1],
#                       [4, 4, 5, 2, 1],
#                       [5, 6, 8, 5, 3]])

# h_weights = np.array([[3, 2, 4, 0],
#                       [3, 2, 4, 2],
#                       [0, 7, 3, 3],
#                       [3, 3, 0, 2],
#                       [1, 3, 2, 2]])
# print(manhattan_tourist(i, j, v_weights, h_weights))

# graph_input = ["1 -> 2",
#                "2 -> 3",
#                "4 -> 2",
#                "5 -> 3"]
# graph = parse_graph(graph_input)
# print(topological_ordering(graph))

weighted_graph_input = ["0->1:7",
                        "0->2:4",
                        "2->3:2",
                        "1->4:1",
                        "3->4:3"]
graph = parse_weighted_graph(weighted_graph_input)
# string1 = "AACCTTGG"
# string2 = "ACACTGTGA"
print(graph)
# print(topological_ordering_weighted_graph(graph))
# print(count_weighted_edges(graph))
print(longest_path(graph, 0, 4))
# backtrack = lcs_backtrack(string1, string2)
# print(backtrack)
# print(output_lcs(backtrack, string1, len(string1), len(string2)))
