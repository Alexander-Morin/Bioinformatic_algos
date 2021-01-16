# Inspired by two stack posts regarding finding the most common element in a list in python. Testing implementations as
# apparently a theoretically slower algorithm that is written in CPython can outperform a 'faster' algo written in pure
# python for lists of certain sizes
# https://stackoverflow.com/a/47845960
# https://stackoverflow.com/a/6988979
# https://www.geeksforgeeks.org/timeit-python-examples/

# Conclusion: max(nucleotides, key=nucleotides.count) implementation quickest for small list (list of 5 3mers), but is
# massively outperformed by statistics.mode() and Counter.most_common(1) when the list is 50 3mers

import timeit
import statistics
from collections import Counter
import numpy as np


def get_motif_array(motifs):
    array = np.empty((len(motifs), len(motifs[0])), dtype=str)  # initialize empty 2d array
    for i in range(len(motifs)):
        for j in range(len(motifs[i])):
            array[i, j] = motifs[i][j]
    return array


# O(n^2) as must compare every element of the list, but apparently written in C
def consensus_max_count(motifs):
    consensus = []
    motif_array = get_motif_array(motifs)
    for i in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, i])
        consensus.append(max(nucleotides, key=nucleotides.count))
    return consensus


# Uncertain how statistics.mode() is implemented
# Note in current version (3.7) throws error if there is no mode - in 3.8 apparently takes first
def consensus_stats_mode(motifs):
    consensus = []
    motif_array = get_motif_array(motifs)
    for i in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, i])
        consensus.append(statistics.mode(nucleotides))
    return consensus


# Uncertain how the Counter object is implemented, but it is specialized for this usage
# Note in current version (3.7) throws error if there is no mode - in 3.8 apparently takes first
def consensus_counter_mostcommon(motifs):
    consensus = []
    motif_array = get_motif_array(motifs)
    for i in range(motif_array.shape[1]):
        nucleotides = list(motif_array[:, i])
        consensus.append(Counter(nucleotides).most_common(1))
    return consensus


# define timeit call functions

def max_count_short_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_max_count'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"]
consensus_max_count(motif_list)'''

    # timeit.repeat statement
    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Max_count short motif list search time: {}'.format(min(times)))


def max_count_long_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_max_count'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"] * 10
consensus_max_count(motif_list)'''

    # timeit.repeat statement
    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Max_count short motif list search time: {}'.format(min(times)))


def stats_mode_short_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_stats_mode
import statistics'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"]
consensus_stats_mode(motif_list)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Statistics.mode() short motif list search time: {}'.format(min(times)))


def stats_mode_long_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_stats_mode
import statistics'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"] * 10
consensus_stats_mode(motif_list)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Statistics.mode() long motif list search time: {}'.format(min(times)))



def counter_mostcommon_short_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_counter_mostcommon
from collections import Counter'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"]
consensus_counter_mostcommon(motif_list)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Counter most common short motif list search time: {}'.format(min(times)))


def counter_mostcommon_long_time():
    SETUP_CODE = ''' 
from __main__ import get_motif_array, consensus_counter_mostcommon
from collections import Counter'''

    TEST_CODE = ''' 
motif_list = ["ACT", "ACT", "AAT", "GGT", "CCT"] * 10
consensus_counter_mostcommon(motif_list)'''

    times = timeit.repeat(setup=SETUP_CODE,
                          stmt=TEST_CODE,
                          repeat=3,
                          number=10000)

    print('Counter most common long list search time: {}'.format(min(times)))


if __name__ == "__main__":
    max_count_short_time()
    max_count_long_time()
    stats_mode_short_time()
    stats_mode_long_time()
    counter_mostcommon_short_time()
    counter_mostcommon_long_time()
