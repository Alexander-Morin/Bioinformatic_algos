# Chapter 2 is focused on finding DNA motifs

import ch1_overview as ch1
import numpy as np

# Review: There are 4^k possible kmers. There are L - k + 1 kmers in a string of length L.
# So, the expected count of occurrences of a kmer appearing in 500x 1000bp sequences is:
# (1 / 4^9) * 500 * (1000 - 9 + 1) = 1.89

# Brute force solution for finding all (k, d) motifs in a collection of strings. In this implementation, a list of
# neighbourhoods (which themselves are a list) is generated. Because the motif must occur in the d-neighbourhood of all
# strings in dna, can intersect the sets of neighbourhoods
# ----------------------------------------------------------------------------------------------------------------------


def motif_enumeration(dna, k, d):
    all_neighbourhoods = [[] for i in range(len(dna))]
    for i in range(len(dna)):
        for j in range(len(dna[i]) - k+1):
            kmer = dna[i][j:j+k]
            neighbourhood = ch1.neighbours(kmer, d)
            for pattern in neighbourhood:
                all_neighbourhoods[i].append(pattern)
    result = set(all_neighbourhoods[0])
    for s in all_neighbourhoods[1:]:
        result.intersection_update(s)
    return " ".join(result)


# However, this requires that every DNA string must contain a (k, d) motif, otherwise it will not be selected. To move
# towards a better solution, consider t DNA strings of length n, and select a kmer from each string to form a
# t x k motif matrix. From this matrix, find the most frequent nucleotide per column. The end goal is to find a means
# to find the most conserved motif. First, we must define a score, which we want to minimize. This involves:
# A 4 x k count matrix, count(motifs), giving the column-wise occurrences of the nucleotides
# A profile matrix, profile(motifs), which divides these counts by t, the number of rows in motifs
# A consensus string, consensus(motifs), of the most popular motifs in each column of the motif matrix
