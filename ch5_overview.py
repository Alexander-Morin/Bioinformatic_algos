# This chaper explores how we compare DNA sequences using dynamic programming

import numpy as np

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



# Example inputs
# ----------------------------------------------------------------------------------------------------------------------


coins = [120, 40, 30, 24, 20, 10, 5, 4, 1]
print(greedy_change(48, coins))