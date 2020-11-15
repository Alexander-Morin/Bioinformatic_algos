# Implementing code to implement the Monte Carlo randomized motif search
# Problem 2F in the BALA textbook/Rosalind

# Input is a text file with first line as k (length of kmer) and t (count of strings) and subsequent lines corresponding
# to strings in DNA
# 8 5
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

# Usage: python3 2F_randomized_motif_search.py input.txt > output.txt
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


def index_to_base(i):
    """ given an index (assumed 0-3), give the corresponding nucleotide in the order of A/C/G/T"""
    return "ACGT"[i]


def base_to_index(b):
    """ given an A/C/G/T nucleotide, return the corresponding index 0-3"""
    return "ACGT".index(b)


def most_probable_kmer(text, k, profile):
    """
    text: a DNA string
    k: an int (ideally 3 <= k <= 12) representing length of kmer
    profile: the motif profile matrix
    returns the most probable kmer in text, given the probabilities contained in profile
    """
    kmer = text[0:k]
    most_prob = 0.0
    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        prob = 1.0
        for j in range(len(pattern)):
            prob *= profile[base_to_index(pattern[j]), j]
        if prob > most_prob:
            most_prob = prob
            kmer = pattern
    return kmer


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


def randomized_motif_search(dna, k, t):
    """
    dna: list of equal length strings
    k: integer (practically >= 3, <= 12)
    t: integer that is the count of strings in dna
    returns a string (motif) and the motif score

    initializes random motifs from dna and gets the motif score. creates a motif profile from these random motifs,
    which is used to find the most probable kmer for each string in dna. this collection of motifs are scored and
    compared to the initial score. this process is repeated until the collection of generated motifs do not beat the
    current best score
    """
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


def repeated_randomized_motif_search(dna, k, t, r=1000):
    """
    dna: list of equal length strings
    k: integer (practically >= 3, <= 12)
    t: integer that is the count of strings in dna
    r: integer that is the amount of times the random motif process should be repeated
    returns a list of the motifs corresponding to the best (lowest) score across all runs

    executes randomized motif search r times, tracking the overall best score and motifs
    """
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    best_score = score_motifs(best_motifs)
    for i in range(0, r):
        motifs, score = randomized_motif_search(dna, k, t)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


def main():
    """
    Read the input file and parse the arguments. Calls repeated randomized motif search
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
            else:
                input_text.append(line.replace('\n', ''))

    output = repeated_randomized_motif_search(input_text, input_k, input_t, 1000)

    for string in output:
        print(string)


if __name__ == "__main__":
    main()
