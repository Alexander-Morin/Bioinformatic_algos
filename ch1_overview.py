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


def base_to_index(b):
    assert b in ["A", "C", "G", "T"], "Invalid base"
    return "ACGT".index(b)


def index_to_base(i):
    assert i in [0, 1, 2, 3], "Invalid index"
    return "ACGT"[i]


def pattern_to_number(pattern):
    if len(pattern) == 1:
        return base_to_index(pattern)
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + base_to_index(symbol)


def number_to_pattern(index, k):
    if k == 1:
        return index_to_base(index)
    prefix_index = index // 4
    remainder = index % 4
    symbol = index_to_base(remainder)
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
    freq_dict = get_freq_dict(text, k)
    freq_pattern = []
    max_count = max(freq_dict.values())
    for kmer in freq_dict:
        if freq_dict[kmer] == max_count:
            freq_pattern.append(kmer)
    return freq_pattern


# From here, the idea of finding clumped kmers within a window is introduced in order to find clusters of DNaA boxes for
# the oriC. A naive solution defines that a kmer of size k forms an (l, t)-clump, where L the size of a window to search
# for clusters of kmers and t is the minimum amount of time the kmer must appear.
# ----------------------------------------------------------------------------------------------------------------------


def find_clump_kmers(genome, l, k, t):
    freq_patterns = []
    for i in range(0, len(genome) - l+1):
        window = genome[i:i+l]
        freq_dict = get_freq_dict(window, k)
        for kmer in freq_dict:
            if freq_dict[kmer] >= t and kmer not in freq_patterns:
                freq_patterns.append(kmer)
    return freq_patterns


# However, this is inefficient as the windows overlap, such that only the first and last kmers of adjacent windows
# differ. Leverage this fact: only build the frequency array once, then update the counts as we slide down genome
# ----------------------------------------------------------------------------------------------------------------------


def faster_clump_kmers(genome, l, k, t):
    freq_patterns = []
    clump_array = np.array([0] * 4**k)
    text = genome[0:l]
    freq_array = get_freq_array(text, k)

    for i in range(0, len(freq_array)):
        if freq_array[i] >= t:
            clump_array[i] = 1

    for i in range(1, len(genome) - l+1):
        first_pattern = genome[i-1:i-1+k]
        index = pattern_to_number(first_pattern)
        freq_array[index] -= 1
        last_pattern = genome[i+l-k:i+l]
        index = pattern_to_number(last_pattern)
        freq_array[index] += 1
        if freq_array[index] >= t:
            clump_array[index] = 1

    for i in range(len(clump_array)):
        if clump_array[i] == 1:
            freq_patterns.append(number_to_pattern(i, k))

    return freq_patterns


# To this point of the chapter, the focus was finding oriC in a bacterial genome. The authors then note that you can
# apply the clump finding approach and still find a wealth of results. They then introduce the fact that DNA replication
# in prokaryotes is done via replication forks, and the asynchronous nature of the forward and reverse strands results
# in differential deamination of C->T (which results in a T-G pairing that gets repaired to T-A). You can use this
# skew to further pinpoint sites of replication.
# ----------------------------------------------------------------------------------------------------------------------


def get_skew_array(genome):
    skew_array = np.array([0] * (len(genome) + 1))
    for i in range(len(genome)):
        if genome[i] == "C":
            skew_array[i+1] = skew_array[i] - 1
        elif genome[i] == "G":
            skew_array[i + 1] = skew_array[i] + 1
        else:
            skew_array[i + 1] = skew_array[i]
    return skew_array


def get_min_skew(genome):
    skew_array = get_skew_array(genome)
    return np.where(skew_array == skew_array.min())


# Until now, the pattern must have been an exact match. The book mentions that an DNaA box may have a slight mismatch
# to the queried pattern. Introduces Hamming distance as a means of comparing the differences between two strings. Use
# this distance to threshold how much the strings are allowed to differ, and find the most frequent kmer while allowing
# mismatches. Note that this means that the most frequent kmer may have no exact matches.
# ----------------------------------------------------------------------------------------------------------------------


def get_hamming_distance(str1, str2):
    # return sum([pattern1[i] != pattern2[i] for i in range(0, len(pattern1))])
    str1 = list(str1)
    str2 = list(str2)
    assert len(str1) == len(str2)
    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist


def approx_pattern_match(pattern, text, d):
    ix = []
    k = len(pattern)
    for i in range(0, len(text) - k+1):
        kmer = text[i:i+k]
        if get_hamming_distance(kmer, pattern) <= d:
            ix.append(i)
    return ix


def approx_pattern_count(pattern, text, d):
    count = 0
    k = len(pattern)
    for i in range(0, len(text) - k + 1):
        kmer = text[i:i + k]
        if get_hamming_distance(kmer, pattern) <= d:
            count += 1
    return count


# A naive approach would be to generate all of the 4**k kmers and see how many approximately match the pattern, but this
# is highly inefficient, as many of these kmers, as well as their d-mismatches, have no match to the pattern. Instead,
# consider the 'd-neighbourhood' of kmers within pattern, which permutes the allowed (up to d mismatches) approx
# matches. Can implement a recursive function that strips to the 'suffix' of the pattern. If we consider a (k-1)mer
# pattern' that belongs to the neighbours(suffix(pattern), d), we know that the Hamming distance between pattern' and
# suffix(pattern) is equal to or less than d. If it's equal, we can add the first symbol of pattern to pattern' to
# obtain a kmer belonging to neighbours(pattern, d). If the distance is less than d, we can add any symbol to the start
# of pattern' and obtain a kmer belonging to neighbours(pattern, d).
# ----------------------------------------------------------------------------------------------------------------------


def immediate_neighbours(pattern):
    neighbourhood = []
    pattern = list(pattern)
    for i in range(0, len(pattern)):
        nucleotides = ["A", "C", "G", "T"]
        symbol = pattern[i]
        nucleotides.remove(symbol)
        for j in nucleotides:
            neighbour = pattern
            neighbour[i] = j
            neighbourhood.append("".join(neighbour))
    return neighbourhood


def neighbours(pattern, d):
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighbourhood = set()
    suffix_neighbours = neighbours(pattern[1:], d)
    for text in suffix_neighbours:
        if get_hamming_distance(pattern[1:], text) < d:
            for nucleotide in ["A", "C", "G", "T"]:
                neighbourhood.add(nucleotide + text)
        else:
            neighbourhood.add(pattern[0] + text)
    return neighbourhood


# Can now create an implementation of finding the most frequent kmers, while allowing for mismatches. Here, instead of
# generating all kmers, we only consider those that are 'close' to kmer in text (has up to d mismatches). Use an array
# of size 4^k to keep track of patterns that are within d to some kmer in text. This allows only having to call approx
# pattern count to close kmers instead of all kmers.
# ----------------------------------------------------------------------------------------------------------------------


def freq_kmers_mismatches(text, k, d):
    freq_patterns = []
    count_array = np.array([0] * 4**k)
    close_array = np.array([0] * 4**k)

    for i in range(0, len(text) - k+1):
        neighbourhood = neighbours(text[i:i+k], d)
        for pattern in neighbourhood:
            index = pattern_to_number(pattern)
            close_array[index] = 1

    for i in range(len(close_array)):
        if close_array[i] == 1:
            pattern = number_to_pattern(i, k)
            count_array[i] = approx_pattern_count(pattern, text, d)

    max_count = max(count_array)

    for i in range(len(count_array)):
        if count_array[i] == max_count:
            freq_patterns.append(number_to_pattern(i, k))

    return freq_patterns, max_count


# Can further speed up the mismatch problem using sorting. Create a neighbourhood array that contains every d-neighbour
# of each kmer in text. Build an index and count array of length equal to all of the generated kmers. Populate the
# index array with the pattern to number indices, and as before, sorting this array causes runs of the same kmer to
# clump together, which can be tracked using the count array.
# ----------------------------------------------------------------------------------------------------------------------


def freq_kmers_mismatch_sort(text, k, d):
    freq_patterns = []
    neighbourhood = []
    for i in range(0, len(text) - k+1):
        pattern = text[i:i+k]
        for j in neighbours(pattern, d):
            neighbourhood.append(j)

    index_array = np.array([0] * len(neighbourhood))
    count_array = np.array([0] * len(neighbourhood))

    for i in range(len(neighbourhood)):
        pattern = neighbourhood[i]
        index_array[i] = pattern_to_number(pattern)
        count_array[i] = 1

    index_array[::-1].sort()  # decreasing sort

    for i in range(1, len(neighbourhood)):
        if index_array[i] == index_array[i-1]:
            count_array[i] += count_array[i-1]

    max_count = max(count_array)

    for i in range(len(neighbourhood)):
        if count_array[i] == max_count:
            freq_patterns.append(number_to_pattern(index_array[i], k))

    return freq_patterns, max_count


# Finally, consider the reverse complement (using the dict implementation)
# ----------------------------------------------------------------------------------------------------------------------


def freq_kmers_mismatch_rev(text, k, d):
    freq_patterns = []
    neighbourhood = {}

    for i in range(0, len(text) - k + 1):
        pattern = text[i:i+k]
        rev_pattern = get_reverse_complement(pattern)
        for j in neighbours(pattern, d):
            if j in neighbourhood:
                neighbourhood[j] += 1
            else:
                neighbourhood[j] = 1
        for j in neighbours(rev_pattern, d):
            if j in neighbourhood:
                neighbourhood[j] += 1
            else:
                neighbourhood[j] = 1

    max_count = max(neighbourhood.values())

    for i in neighbourhood:
        if neighbourhood[i] == max_count:
            freq_patterns.append(i)

    return freq_patterns
