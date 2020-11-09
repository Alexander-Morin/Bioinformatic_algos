# Implementing code to use a greedy algorithm for finding the 'profile-most-probable' kmers in a collection of DNA
# strings. Returns the most likely kmer for each string, building the profile matrix with previous 'best kmers.'
# Problem 2D in the BALA textbook/Rosalind
# Credit to # http://www.mrgraeme.co.uk/greedy-motif-search/ for help on this problem

# Input is a text file with first line as k (length of kmer) and t (count of strings) and subsequent lines corresponding
# to strings in DNA
# 3 5
# GGCGTTCAGGCA
# AAGAATCAGTCA
# CAAGGAGTTCGC
# CACGTCAATCAC
# CAATAATATTCG

# Output is the median string
# CAG
# CAG
# CAA
# CAA
# CAA

# Usage: python3 2D_greedy_motif_search.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from collections import Counter
import numpy as np


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
        a.append(nucleotides.count("A"))
        c.append(nucleotides.count("C"))
        g.append(nucleotides.count("G"))
        t.append(nucleotides.count("T"))
    counts_matrix = np.array([a, c, g, t])
    profile_matrix = counts_matrix / len(motifs)
    return profile_matrix


def most_probable_kmer(text, k, profile):
    """
    text: a DNA string
    k: an int (ideally 3 <= k <= 12) representing length of kmer
    profile: the motif profile matrix
    returns the most probable kmer in text, given the probabilities contained in profile
    """
    kmer = text[0:k]
    best_prob = 0
    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        prob = 1
        for j in range(len(pattern)):
            if pattern[j] == "A":
                prob *= profile[0, j]
            elif pattern[j] == "C":
                prob *= profile[1, j]
            elif pattern[j] == "G":
                prob *= profile[2, j]
            elif pattern[j] == "T":
                prob *= profile[3, j]
        if prob > best_prob:
            kmer = pattern
            best_prob = prob
    return kmer


def greedy_motif_search(dna, k, t):
    """
    dna: list of equal length DNA strings
    k: an int (ideally 3 <= k <= 12) representing length of kmer
    t: count of strings in dna
    returns a list of motifs, one for each string in dna, representing the greedy most probable kmer

    best motifs initialized as the first kmer in each string. for each kmer in the first string, iteratively find the
    profile most probable kmer in the subsequent strings, reconstructing the profile at each step with the newly added
    most probable motifs. then calculate if the score of this new collection of motifs is lower (minimizes motif
    distance) compared to the current best motifs
    """
    best_motifs = np.array([[text[0:k] for text in dna]])
    for i in range(len(dna[0]) - k+1):
        motifs = [dna[0][i:i+k]]
        for j in range(1, t):
            profile = profile_motifs(motifs)
            most_prob = most_probable_kmer(dna[j], k, profile)
            motifs.append(most_prob)
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs


def main():
    """
    Read the input file and parse the arguments. Call greedy motif search.
    """
    argv = list(sys.argv)
    input_text = []
    input_param = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_param is None:
                if input_param is None:
                    input_param = line.replace('\n', '').split()
                    input_k = int(input_param[0])
                    input_t = int(input_param[1])
            else:
                input_text.append(line.replace('\n', ''))

    output = greedy_motif_search(input_text, input_k, input_t)

    for motif in output[0]:
        print(motif)


if __name__ == "__main__":
    main()
