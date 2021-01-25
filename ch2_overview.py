# Chapter 2 is focused on finding DNA motifs

import ch1_overview as ch1
import numpy as np
from collections import Counter
import math
import random

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
    motif_array = get_motif_array(motifs)
    a = []
    c = []
    g = []
    t = []
    for column in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, column])
        a.append(nucleotides.count("A")+1)  # add pseudocounts to stop multiplying probs by 0
        c.append(nucleotides.count("C")+1)
        g.append(nucleotides.count("G")+1)
        t.append(nucleotides.count("T")+1)
    counts_matrix = np.array([a, c, g, t])
    profile_matrix = counts_matrix / sum(counts_matrix[:, 0])
    return profile_matrix


def get_consensus_motifs(motifs):
    consensus = []
    motif_array = get_motif_array(motifs)
    for i in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, i])
        consensus.append(max(nucleotides, key=nucleotides.count))  # faster for small lists
        # consensus.append(Counter(nucleotides).most_common(1))  # faster for larger lists
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
# Note that the following implementation must consider d(pattern, string) for every string in dna, requiring approx
# k * n * t operations for t strings of length n in dna. Therefore O(4**k * n * k * t), which compares favorably
# to the brute force motif search [O(n**t * k * t)] as k is usually under 20, while t is measures in the thousands.
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


# The text then introduces greedy algorithms, noting that they are typically heuristics that sacrifice speed for
# accuracy, but often find good use in biological problems. Introduces the profile most probable kmer problem, where
# given a string text, an integer k, and a 4xk profile matrix, the algorithm will produce the most probable kmer in
# text. Can then use this as the basis of greedy motif search, which looks to find the most probably kmer iteratively,
# reconstructing the profile matrix each time after finding the most prob kmer and adding it to the set of motifs.
# Because some columns may have probability = 0 for a nucleotide, resulting in 0 for the entire motif after multiplying,
# use pseudocounts (Laplace's rule of succession) for correction.
# ----------------------------------------------------------------------------------------------------------------------


def most_probable_kmer(text, k, profile):
    kmer = text[0:k]
    most_prob = 0.0
    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        prob = 1.0
        for j in range(len(pattern)):
            prob *= profile[ch1.base_to_index(pattern[j]), j]
        if prob > most_prob:
            most_prob = prob
            kmer = pattern
    return kmer


def greedy_motif_search(dna, k, t):
    # http://www.mrgraeme.co.uk/greedy-motif-search/
    best_motifs = np.array([list(text[0:k]) for text in dna])
    for i in range(len(dna[0]) - k+1):
        motifs = [dna[0][i: i + k]]
        for j in range(1, t):
            profile = profile_motifs(motifs)
            most_prob = most_probable_kmer(dna[j], k, profile)
            motifs.append(most_prob)
        if score_motifs(motifs) < score_motifs(best_motifs):
                best_motifs = motifs
    return best_motifs


# The book then discusses randomized algorithms: Las Vegas algorithms give an exact answer, while Monte Carlo algos may
# not give an exact answer, but they are quick and can be repeated many times to point to a consensus answer. Introduces
# the randomized motif search, where kmers are randomly chosen to construct a profile, which is used to select the most
# probable kmer. Define a separate function to call this process 1000 times, tracking the best overall score and motif
# ----------------------------------------------------------------------------------------------------------------------


def random_kmer(string, k):
    start = random.randint(0, len(string) - k)
    return string[start:start+k]


def random_motifs(dna, k):
    return [random_kmer(string, k) for string in dna]


def randomized_motif_search(dna, k, t):
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    while True:
        profile = profile_motifs(motifs)
        motifs = [most_probable_kmer(string, k, profile) for string in dna]
        score = score_motifs(motifs)
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs, best_score


def repeated_randomized_motif_search(dna, k, t, r):
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for i in range(0, r):
        motifs, score = randomized_motif_search(dna, k, t)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


# The text notes that this randomized motif search works because the implanted motif ultimately biases the generated
# profile matrix. However, it is 'reckless' as it regenerates the set of motifs at each step, potentially throwing away
# a true motif. It then introduces the Gibbs sampler, which uses a weighted random function to only throw away one of
# the motifs (note that this is random, while the randomized motif search is deterministic after the initial motif
# selection as it selects the profile most probable kmers.

# The weighted random function accepts probabilities that are not constrained by having to sum to 1. It sums the probs
# and weighs each prob by this sum.

# Gibbs sampler also starts by initializing random motifs and finding the score. It then goes through an iterative
# process of (normal) randomly removing one of the strings in dna, as well as the corresponding motif. A profile is then
# made from this subset, and the probabilities from this profile are used in the (weighted) random function to sample a
# motif from the held out dna string. This sampled motif is added back to the collection of motifs, and their score
# compared to the best score. This process must be repeated to avoid getting trapped in local optimums.


def weighted_random_integer(probs):
    weighted_probs = [i / sum(probs) for i in probs]
    return random.choices(range(0, len(probs)), weighted_probs)[0]


def profile_kmer_probs(text, k, profile):
    kmer_probs = dict()
    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        prob = 1.0
        for j in range(len(pattern)):
            prob *= profile[ch1.base_to_index(pattern[j]), j]
        kmer_probs[pattern] = prob
    return kmer_probs


def profile_random_kmer(text, k, profile):
    kmer_probs = profile_kmer_probs(text, k, profile)
    i = weighted_random_integer(kmer_probs.values())
    return list(kmer_probs.keys())[i]


def gibbs_sampler(dna, k, t, n=100):
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for j in range(n):
        test_motifs = best_motifs.copy()
        i = random.randint(0, t-1)
        dna_i = dna[i]
        del test_motifs[i]
        profile = profile_motifs(test_motifs)
        random_motif = profile_random_kmer(dna_i, k, profile)
        test_motifs.insert(i, random_motif)
        score = score_motifs(test_motifs)
        if score < best_score:
            best_motifs = test_motifs
            best_score = score
    return best_motifs, best_score


def repeated_gibbs_sampler(dna, k, t, n=100, r=20):
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for i in range(0, r):
        motifs, score = gibbs_sampler(dna, k, t, n)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs
