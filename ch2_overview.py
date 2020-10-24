# Chapter 2 is focused on finding DNA motifs

import sys
from ch1_overview import *

# Review: There are 4^k possible kmers. There are L - k + 1 kmers in a string of length L.
# So, the expected count of occurrences of a kmer appearing in 500x 1000bp sequences is:
# (1 / 4^9) * 500 * (1000 - 9 + 1) = 1.89

# Brute force solution for finding all (k, d) motifs in a collection of strings
# ----------------------------------------------------------------------------------------------------------------------


def motif_enumeration(dna, k, d):
    all_neighbourhoods = [[] for i in range(len(dna))]
    for i in range(len(dna)):
        for j in range(len(dna[i]) - k+1):
            kmer = dna[i][j:j+k]
            neighbourhood = neighbours(kmer, d)
            for pattern in neighbourhood:
                all_neighbourhoods[i].append(pattern)
    result = set(all_neighbourhoods[0])
    for s in all_neighbourhoods[1:]:
        result.intersection_update(s)
    return " ".join(result)



# print(motif_enumeration(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 3, 1))