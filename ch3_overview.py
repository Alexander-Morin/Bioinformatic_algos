# Chapter 3 is focused on genome assembly

import ch1_overview as ch1
import numpy as np

# First, focus on the more simplistic case where are all reads are kmers of length k, all reads come from the same
# strand, there are no errors, and they exhibit perfect coverage (every kmer substring of the genome has a read).
# Consider the kmer composition of a string (lexicographically sorted kmers of a string) - to model genome assembly
# must perform the inverse problem, where kmers are assembled into a string. The chapter starts by motivating the
# (k-1) overlap of kmers to spell out the string. This task is complicated by repeats and constructing the string's
# path requires knowing the genome in advance
# ----------------------------------------------------------------------------------------------------------------------


def string_composition(text, k):
    kmers = [text[i:i+k] for i in range(0, len(text) - k+1)]
    kmers.sort()
    return kmers


def spell_string_by_path(text):
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


# The problem is then phrased in terms of graphs: Kmers are connected if the suffix of one is equal to the prefix of
# another, resulting in a directed overlap graph where the kmer nodes are sorted lexicographically. Spelling out the
# string means visiting each node exactly once. Introduces the overlap graph problem, with (k-1) overlap and returning
# an adjacency list (which is actually better represented as a mapping => dict)
# ----------------------------------------------------------------------------------------------------------------------


def get_adjacency_dict(patterns):
    adjacency_dict = dict.fromkeys(patterns)
    for node1 in adjacency_dict:
        suffix = node1[1:]
        for node2 in patterns:
            prefix = node2[:-1]
            if suffix == prefix:
                adjacency_dict[node1] = node2
    return adjacency_dict


# To solve the string reconstruction problem, looking for a Hamiltonian path: a path in the graph that visits every
# node. Note that there may be more than one Hamiltonian path in a graph. States that an alternative way to find
# Hamiltonian paths in large graphs is to follow de Bruijn's approach: represent a kmer composition as a graph by
# treating the kmers as edges, rather than nodes. Pairs of consecutive edges represent kmers that overlap in (k-1)
# nucleotides, so the node with the incoming and outgoing edges are labeled as the (k-1)mers. This results in the
# merging of identical nodes into one - decreases count of nodes while preserving the count of edges. For example, if
# text has 3 repeats of ATG, there would be 3 edges connecting node AT to node TG
# ----------------------------------------------------------------------------------------------------------------------


def de_graph(text, k):
    adjacency_dict = dict()
    for i in range(len(text) - k+1):
        kmer = text[i:i+k-1]
        next_kmer = text[i+1:i+k]
        if kmer not in adjacency_dict:
            adjacency_dict[kmer] = [next_kmer]
        else:
            adjacency_dict[kmer].append(next_kmer)
    return adjacency_dict




patterns = [
    "ATGCG",
    "GCATG",
    "CATGC",
    "AGGCA",
    "GGCAT"
    ]




