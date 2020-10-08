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
# text[i:i+k] appears in text. This implementation is slow: pattern_count is called for every kmer/window, which search
# the whole string from the start each time. O(|text|^2 * k).
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
    return np.unique(freq_kmers)


# An alternative approach is provided so that a window must only slide down text once. There are 4^k possible kmers
# ( (length of alphabet) ^ size of k ). Sort these lexicographically, convert each to an integer corresponding to the
# 4 ^ k possible kmers, and create a frequency array. The ith element of the frequency array corresponds to the count of
# times the (lexicographically sorted) kmer appears in text.

# This 4-ary encoding system has A=0 C=1 G=2 T=3. If you remove the last symbol of these sorted kmers, they are still
# sorted lexicographically. These new (k-1)mers are now repeated 4 times. Can define a recursive algorithm that strips
# these patterns until it gets to a base symbol, adding up the recursions to get the final index of the sorted array.
# AAAACT ([0*(4^5)] + [0*(4^4)] + [0*(4^3)] + [0*(4^2)] + [1*(4^1)] + [3*(4^0)]) = 7, where 7 is the index

# Reverse: when you divide index = pattern_to_number(prefix) by 4, the remainder will be equal to symbol to number and
# the quotient will be equal to pattern to number (prefix).


def symbol_to_number(symbol):
    assert symbol in ["A", "C", "G", "T"], "Invalid nucleotide"
    if symbol == "A":
        return 0
    elif symbol == "C":
        return 1
    elif symbol == "G":
        return 2
    elif symbol == "T":
        return 3


def pattern_to_number(pattern):
    if len(pattern) == 1:
        return symbol_to_number(pattern)
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)


def number_to_symbol(number):
    assert number in [0, 1, 2, 3], "Invalid index"
    if number == 0:
        return "A"
    elif number == 1:
        return "C"
    elif number == 2:
        return "G"
    elif number == 3:
        return "T"


def number_to_pattern(index, k):
    if k == 1:
        return number_to_symbol(index)
    prefix_index = index // 4
    remainder = index % 4
    symbol = number_to_symbol(remainder)
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    return prefix_pattern + symbol


def get_freq_array(text, k):
    freq_array = np.array([0] * 4**k)
    for i in range(0, len(text) - k+1):
        pattern = text[i:i+k]
        index = pattern_to_number(pattern)
        freq_array[index] += 1
    return freq_array


# Put it all together to implement a faster version of finding the most frequent kmers using a frequency array
# However, this implementation is only faster for small values of k.
# Essentially, FrequentWords (FW) is (O((|Text|^2) * k)) and FasterFrequentWords (FFW) is (O((|Text|*k) + ((4^k) * k))).
# For all values of k such that 4^k  < |Text|^2 - |Text|, FFW will always be faster (lesser steps to execute) than FW.
# For all values of  k such that 4^k  > |Text|^2 - |Text|, FFW will always be slower (more steps to execute) than FW.

# Finally, an implementation that first sorts the frequency is introduced. Given a string, list all of its kmers in the
# order they appear, and convert these to their index from symbol to number. Sorting this array causes all of the same
# kmers (via their indices) to clump together. So the most frequent kmers are those that have the longest run.

# However, the simplest implementation is to take advantage of python's dictionary data structure
# ----------------------------------------------------------------------------------------------------------------------


def faster_freq_words(text, k):
    freq_patterns = []
    freq_array = get_freq_array(text, k)
    max_count = max(freq_array)
    for i in range(0, len(freq_array)):
        if freq_array[i] == max_count:
            freq_patterns.append(number_to_pattern(i, k))
    return freq_patterns


def sorted_freq_words(text, k):
    freq_patterns = []
    count_array = np.array([0] * 4**k)
    index_array = np.array([0] * 4**k)

    for i in range(len(text) - k+1):
        pattern = text[i:i+k]
        index_array[i] = pattern_to_number(pattern)
        count_array[i] = 1

    index_array[::-1].sort()  # decreasing sort

    for i in range(1, len(text) - k+1):
        if index_array[i] == index_array[i-1]:
            count_array[i] += count_array[i-1]

    max_count = max(count_array)

    for i in range(0, len(text) - k+1):
        if count_array[i] == max_count:
            freq_patterns.append(number_to_pattern(index_array[i], k))

    return freq_patterns


def get_freq_dict(text, k):
    freq_dict = {}
    for i in range(0, len(text) - k+1):
        pattern = text[i:i+k]
        if pattern not in freq_dict:
            freq_dict[pattern] = 1
        else:
            freq_dict[pattern] += 1
    return freq_dict


def get_freq_kmer_dict(text, k):

