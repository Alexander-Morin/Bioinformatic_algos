# Chapter 2 is focused on finding DNA motifs

import ch1_overview as ch1
import numpy as np
from collections import Counter
import math

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
# to get the most conserved motif. First, we must define a score, which we want to minimize. This involves:

# Representing the list of strings in DNA as an array, where [i, j] is nucleotide j of string i.
# A scoring system, which acts column-wise by summing the amount of nucleotides that differ from the dominant nucleotide
# A 4 x k count matrix, count(motifs), giving the column-wise occurrences of the nucleotides
# A profile matrix, profile(motifs), which divides these counts by t, the number of rows in motifs
# A consensus string, consensus(motifs), of the most popular motifs in each column of the motif matrix
# ----------------------------------------------------------------------------------------------------------------------


def get_motif_array(motifs):
    array = np.empty((len(motifs), len(motifs[0])), dtype=str)  # initialize empty 2d array
    for i in range(len(motifs)):
        for j in range(len(motifs[i])):
            array[i, j] = motifs[i][j]
    return array


def score_motifs(motifs):
    motif_array = get_motif_array(motifs)
    scores = []
    for i in range(motif_array.shape[1]):
        counts = Counter(motif_array[:, i])
        scores.append(len(motif_array[:, 0]) - max(counts.values()))
    return sum(scores)


def count_motifs(motifs):
    motif_array = get_motif_array(motifs)
    a = []
    c = []
    g = []
    t = []
    for column in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, column])
        a.append(nucleotides.count("A"))
        c.append(nucleotides.count("C"))
        g.append(nucleotides.count("G"))
        t.append(nucleotides.count("T"))
    return np.array([a, c, g, t])


def profile_motifs(motifs):
    counts = count_motifs(motifs)
    return counts / len(motifs)


def get_consensus_motifs(motifs):
    consensus = []
    motif_array = get_motif_array(motifs)
    for i in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, i])
        consensus.append(max(nucleotides, key=nucleotides.count))
    return consensus


# However, this consensus approach isn't satisfying: consider a column where there are 6 Cs, 2 As and 2 Ts, versus
# a column with 6 Cs and 4 Ts. Both are contributing scores of 4, yet one is more conserved.
# So switch to using entropy as a measure of uncertainty of the probability of each nucleotide's occurrence.
# H(p1...pN) = -sum(pi * log2(pi))
# In turn, these entropies are used to derive information content (1 - H(p1...pN) which is used in motif logos
# ----------------------------------------------------------------------------------------------------------------------


def get_entropy(prob_vec):
    prob_vec = prob_vec[prob_vec != 0]  # remove 0s b/c of log2  - count these as 0 anyway
    values = map(lambda x: x * math.log(x, 2), prob_vec)
    return round(-sum(list(values)), 3)


def get_motifs_entropy(motifs):
    entropy = []
    profile = profile_motifs(motifs)
    for j in range(profile.shape[1]):
        entropy.append(get_entropy(profile[:, j]))
    return entropy


# The text then introduces the Motif Finding Problem: Given a collection of strings, find a set of k-mers, one from each
# string, that minimizes the score of the resulting motif.
# A brute force solution would consider every kmer in the collection of of strings (each containing n nucleotides).
# As there are n - k + 1 motifs per string, and a total of t strings, there are (n - k + 1)**t ways to form motifs.
# For each choice of motifs the scoring function must get called, which requires k * t steps. Assuming k is smaller than
# n, this approach would have O(n**t * k * t)

# The text then states that rather than trying to derive motifs from each string, and then find consensus among this
# collection of motifs, to instead explore all potential consensus kmers, and then find the best collection of motifs
# for each consensus string. So consensus(motifs) -> motifs

# This change of strategy requires finding a new scoring system. Originally, was considering column-wise to find the
# difference from the consensus. Now consider row wise, where each string gets compared to a consensus string, and the
# row score is the Hamming distance. The final score is the sum of these row scores.
# ---------------------------------------------------------------------------------------------------------------------


def pattern_distance(pattern, dna):
    scores = []
    for string in dna:
        scores.append(ch1.get_hamming_distance(string, pattern))
    return sum(scores)


# The text then notes that this row-wise scoring (when compared against the consensus motif/pattern) yields the same
# score as the column-wise approach (which is also comparing the difference from the most common nucleotide). This
# "gives us an idea" - rather than searching for a collection of motifs minimizing the score, instead search for a
# potential consensus string pattern that minimizes pattern_distance among all possible k-mers 'pattern' and all
# possible choices of kmers in dna. This is equivalent to the motif finding problem.

# Equivalent Motif Finding Problem: Given a collection of strings, find a collection of k-mers (one from each string)
# that minimizes the distance between all possible patterns and all possible collections of k-mers.

# The text notes that on the surface, this equivalent problem seems harder, as we must consider all motifs as well as
# all kmer patterns. The 'key observation' is that if we are given pattern, we don't need to explore all possible
# collections of motif in order to minimize d(pattern, motifs)

# Define Motifs(Pattern, dna) as the set of motifs (one per string in dna) that minimizes the Hamming Distance between
# pattern and all possible kmers in dna. d(pattern, text) refers to this minimum Hamming distance. Note that there can
# be multiple kmers that minimize this distance, but this isn't an issue. d(pattern, dna) is therefore the sum of these
# distances over all text/strings in dna.

# The goal then becomes to find the 'median string' - the kmer pattern that minimizes d(pattern, dna) over all patterns.
# ----------------------------------------------------------------------------------------------------------------------


def get_pattern_distance_from_dna(pattern, dna):
    k = len(pattern)
    distance = 0
    for string in dna:
        hdist = math.inf
        for i in range(len(string) - k + 1):
            kmer = string[i:i+k]
            if hdist > ch1.get_hamming_distance(pattern, kmer):
                hdist = ch1.get_hamming_distance(pattern, kmer)
        distance += hdist
    return distance


def median_string(dna, k):
    distance = math.inf
    for i in range(4**k):
        pattern = ch1.number_to_pattern(i, k)
        pattern_dist = get_pattern_distance_from_dna(pattern, dna)
        if distance > pattern_dist:
            distance = pattern_dist
            median = pattern
    return median



motifs = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "ATGGGGACTTCC",
    "TCGGGGACTTCC",
    "TCGGGGATTCAT",
    "TAGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
]


print(median_string(motifs, 5))