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



# Example inputs
# ----------------------------------------------------------------------------------------------------------------------

# print(greedy_change(48, [120, 40, 30, 24, 20, 10, 5, 4, 1]))
# print(recursive_change(40, [5, 4, 1]))
print(dp_change(40, [1, 5, 10, 20, 25, 50]))
