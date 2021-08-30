# This chaper explores how we compare DNA sequences using dynamic programming

import numpy as np
import pandas as pd

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
# recrusive and re-calculates steps many times. The DP approach calculates every sub-step to find the ultimate max.
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


# DP
def manhattan_tourist(i, j, v_edge_weight, h_edge_weight):
    scores = np.zeros([i+1, j+1])
    for x in range(1, i+1):
        scores[x][0] = scores[x-1][0] + v_edge_weight[x-1][0]
    for y in range(1, j+1):
        scores[0][y] = scores[0][y-1] + h_edge_weight[0][y-1]
    for x in range(1, i+1):
        for y in range(1, j+1):
            scores[x][y] = max(scores[x-1][y] + v_edge_weight[x-1][y], scores[x][y-1] + h_edge_weight[x][y-1])
    return int(scores[i][j])


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
#
# h_weights = np.array([[3, 2, 4, 0],
#                       [3, 2, 4, 2],
#                       [0, 7, 3, 3],
#                       [3, 3, 0, 2],
#                       [1, 3, 2, 2]])
# print(manhattan_tourist(i, j, v_weights, h_weights))
