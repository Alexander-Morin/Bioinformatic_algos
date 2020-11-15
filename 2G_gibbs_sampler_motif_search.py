# Implementing the Gibbs sampler motif search, a Monte Carlo algorithm which looks to find the set of motifs within a
# collection of DNA strings that minimizes the distance of all possible motifs within dna
# Problem 2G in the BALA textbook/Rosalind

# Input is a text file with first line as k (length of kmer), t (count of strings), n (number of iterations) and
# subsequent lines corresponding to strings in DNA
# 8 5 100
# CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
# GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
# TAGTACCGAGACCGAAAGAAGTATACAGGCGT
# TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
# AATCCACCAGCTCCACGTGCAATGTTGGCCTA

# Output
# TCTCGGGG
# CCAAGGTG
# TACAGGCG
# TTCAGGTG
# TCCACGTG

# Usage: python3 2G_gibbs_sampler_motif_search.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from collections import Counter
import numpy as np
import random


def get_motif_array(motifs):
    """
    motif: list of equal length DNA strings
    returns a t (count of strings) x n (length of string) numpy array of the aligned nucleotides
    """
    array = np.empty((len(motifs), len(motifs[0])), dtype=str)
    for i in range(len(motifs)):
        for j in range(len(motifs[i])):
            array[i, j] = motifs[i][j]
    return array


def score_motifs(motifs):
    """
    motifs: list of equal length DNA strings
    returns an int corresponding to the sum of how many nucleotides differ from the dominant nucleotide across each
    column of the motif array
    """
    motif_array = get_motif_array(motifs)
    scores = []
    for i in range(motif_array.shape[1]):
        counts = Counter(motif_array[:, i])
        # length of the column minus max count gives count of differing nucleotides for that column
        scores.append(len(motif_array[:, 0]) - max(counts.values()))
    return sum(scores)


def profile_motifs(motifs):
    """
    motif: list of equal length DNA strings
    returns a 4 (nucleotides A, C, G, T) by n (length of string) numpy array of the ratio of each aligned nucleotide
    """
    motif_array = get_motif_array(motifs)
    a = []
    c = []
    g = []
    t = []
    for column in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, column])
        a.append(nucleotides.count("A")+1)  # add pseudocounts to stop multiplying by 0
        c.append(nucleotides.count("C")+1)
        g.append(nucleotides.count("G")+1)
        t.append(nucleotides.count("T")+1)
    counts_matrix = np.array([a, c, g, t])
    profile_matrix = counts_matrix / sum(counts_matrix[:, 0])
    return profile_matrix


def random_kmer(string, k):
    """
    returns a random kmer of length k from string
    """
    start = random.randint(0, len(string) - k)
    return string[start:start+k]


def random_motifs(dna, k):
    """"
    returns a list of motifs, taking a random kmer of length k from each string in dna
    """
    return [random_kmer(string, k) for string in dna]


def weighted_random_integer(probs):
    """
    probs: list of non-negative probabilities
    returns a weighted random index corresponding to the given probs, which are weighted by the sum of the probs
    """
    weighted_probs = [i / sum(probs) for i in probs]
    return random.choices(range(0, len(probs)), weighted_probs)[0]


def index_to_base(i):
    """ given an index (assumed 0-3), give the corresponding nucleotide in the order of A/C/G/T"""
    return "ACGT"[i]


def base_to_index(b):
    """ given an A/C/G/T nucleotide, return the corresponding index 0-3"""
    return "ACGT".index(b)


def profile_kmer_probs(text, k, profile):
    """
    text: a DNA string
    k: an int (ideally 3 <= k <= 12) representing length of kmer
    profile: the motif profile matrix
    returns a dict of the profile-driven probabilities of the kmers in text
    """
    kmer_probs = dict()
    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        prob = 1.0
        for j in range(len(pattern)):
            prob *= profile[base_to_index(pattern[j]), j]
        kmer_probs[pattern] = prob
    return kmer_probs


def profile_random_kmer(text, k, profile):
    """
    returns a weighted random kmer in text, based on the profile matrix
    """
    kmer_probs = profile_kmer_probs(text, k, profile)
    i = weighted_random_integer(kmer_probs.values())
    return list(kmer_probs.keys())[i]


def gibbs_sampler(dna, k, t, n=100):
    """
    dna: list of equal length strings
    k: integer (practically >= 3, <= 12)
    t: integer that is the count of strings in dna
    n: integer of the number of iterations
    returns a string (motif) and the motif score

    initializes random motifs from dna and gets the motif score. in each iteration, a random string from dna (and the
    corresponding motif within) are removed, and a motif profile generated from the remaining motifs. the probabilities
    from this profile are used in a weighted random function, which is used to choose a random motif from the held out
    dna string. this random motif is then inserted back into the current set of motifs, and the score compared to the
    best set of motifs, which is updated if the new set of motifs has a lower score
    """
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for j in range(1, n):
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
    """
    dna: list of equal length strings
    k: integer (practically >= 3, <= 12)
    t: integer that is the count of strings in dna
    n: integer of the number of iterations
    r: integer of the number of times gibb sampler should be called
    returns a list of strings (best motifs) equal to length of dna

    calls gibbs sampler r times, tracking the global best motifs
    """
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for i in range(0, r):
        motifs, score = gibbs_sampler(dna, k, t, n)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


def main():
    """
    Read the input file and parse the arguments. Calls repeated gibbs sampler motif search
    """
    argv = list(sys.argv)
    input_text = []
    input_param = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_param is None:
                input_param = line.replace('\n', '').split()
                input_k = int(input_param[0])
                input_t = int(input_param[1])
                input_n = int(input_param[2])
            else:
                input_text.append(line.replace('\n', ''))

    output = repeated_gibbs_sampler(input_text, input_k, input_t, input_n, 20)

    for string in output:
        print(string)


if __name__ == "__main__":
    main()
