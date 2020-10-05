# Chapter 1 is centered around text matching algorithms, framed by asking how to find the ori in bacteria.

import numpy as np

# Begin with simple string match/count, sliding a window of length pattern down text. Also introduce function to
# retrieve the reverse complement of a DNA string.
# ----------------------------------------------------------------------------------------------------------------------

def pattern_count(text, pattern):
    count = 0
    k = len(pattern)
    for i in range(0, len(text) - k-1):
        if text[i:i+k] == pattern:
            count += 1
    return count


def pattern_match(text, pattern):
    index_array = []
    k = len(pattern)
    for i in range(0, len(text) - k-1):
        if text[i:i + k] == pattern:
            index_array.append(i)
    return index_array


def get_reverse_complement(text):
    complement = []
    for letter in text:
        if letter == 'A':
            complement.append('T')
        elif letter == 'C':
            complement.append('G')
        elif letter == 'T':
            complement.append('A')
        elif letter == 'G':
            complement.append('C')
    return "".join(complement[::-1])


# Then introduces a naive algorithm for finding frequent kmer: which patterns of length k appear most frequently in a
# text string. Store counts in an array, where each element of the count array corresponds to the amount of times that
# text[i:i+k] appears in text.
# ----------------------------------------------------------------------------------------------------------------------

def init_count_array(text):
    text = list(text)
    count = [0] * len(text)
    array = np.array([text, count])
    return array


def get_freq_kmer(text, k):
    freq_kmers = []
    count_array = init_count_array(text)
    for i in range(0, len(text) - k-1):
        pattern = text[i:i+k]
        count_array[1, i] = pattern_count(text, pattern)
    max_count = max(count_array[1])
    for i in range(0, len(text) - k-1):
        if count_array[1, i] == max_count:
            freq_kmers.append(text[i:i+k])
    return(np.unique(freq_kmers))


# This implementation is slow: pattern_count is called for every kmer/window, which search the whole string from the
# start each time. O(|text|^2 * k).
# An alternative approach is provided so that a window must only slide down text once. There are 4^k possible kmers
# ( (length of alphabet) ^ size of k ). Sort these lexicographically, convert each to an integer corresponding to the
# 4 ^ k possible kmers, and create a frequency array. The ith element of the frequency array

