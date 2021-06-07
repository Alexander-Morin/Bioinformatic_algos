# Chapter 3 is focused on genome assembly
# https://doi.org/10.1038/nbt.2023
# https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf

from random import choice
import pandas as pd
from collections import defaultdict, deque

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


def overlap_graph(patterns):
    adjacency_dict = dict.fromkeys(patterns)
    for node1 in adjacency_dict:
        suffix = node1[1:]
        for node2 in patterns:
            prefix = node2[:-1]
            if suffix == prefix:
                adjacency_dict[node1] = node2
    return adjacency_dict


# To solve the string reconstruction problem, looking for a Hamiltonian path: a path in the graph that visits every
# node once. Note that there may be more than one Hamiltonian path in a graph. States that an alternative way to find
# Hamiltonian paths in large graphs is to follow de Bruijn's approach: represent a kmer composition as a graph by
# treating the kmers as edges, rather than nodes. Pairs of consecutive edges represent kmers that overlap in (k-1)
# nucleotides, so the node with the incoming and outgoing edges are labeled as the (k-1)mers. This results in the
# merging of identical nodes into one - decreases count of nodes while preserving the count of edges. For example, if
# text has 3 repeats of ATG, there would be 3 edges connecting node AT to node TG
# ----------------------------------------------------------------------------------------------------------------------


def dbg_from_string(text, k):
    patterns = [text[i:i+k] for i in range(len(text))]
    adjacency_dict = defaultdict(list)
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_dict[prefix].append(suffix)
        adjacency_dict[suffix]  # init suffix in case it has no outgoing edges
    return adjacency_dict


# Solving the string reconstruction problem reduces to finding a path that visits every edge of the De Bruijn graph
# exactly once == Eulerian path. The text then talks about the composition graph of text, which has all of the kmers
# as isolated edges. Gluing nodes with the same label in the composition graph produces DeBruijn(text). The text then
# suggests a way of constructing the De Bruijn graph that does involve gluing: given a collection of kmer patterns, the
# nodes of DeBruijn(patterns) are all unique (k-1)mers occurring as a prefix or suffix. DeBruijn(patterns) forms a
# directed graph by connecting prefix to suffix. However, I am not certain how this is different from constructing
# from a text and k...
# ----------------------------------------------------------------------------------------------------------------------


def debruijn_graph(patterns):
    adjacency_dict = defaultdict(list)
    for kmer in patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_dict[prefix].append(suffix)
        adjacency_dict[suffix]  # init suffix in case it has no outgoing edges
    return adjacency_dict


# Anywho, states that we can now solve the string reconstruction problem by either finding a Hamiltonian path (ever node
# once) in the Overlap graph, or a Eulerian path (every edge once) in the DBG. Their provided visualization shows that,
# for small graphs, the DBG is more attractive as it is smaller. But the key is to find an efficient algorithm that
# solves large graphs. Hamiltonian path is NP complete. Euler's theorem: every balanced (every node has equal in
# degree to out degree) and strongly connected (possible to reach any node from any node) graph is Eulerian, meaning
# you can visit every edge once while starting/finishing on the same node. Because the graph is balanced, you can only
# get stuck on the starting node. The idea is to backtrack to a node in the current cycle that has unused edges and
# form a new circuit, adding this circuit to the previous until all edges have been visited. O(E) - linear over the
# number of edges in the graph.
# ----------------------------------------------------------------------------------------------------------------------


def parse_graph(list_input):
    graph = defaultdict(list)
    for line in list_input:
        split_input = line.split(" -> ")
        in_node = split_input[0].strip()
        out_nodes = split_input[1].split(",")
        for node in out_nodes:
            graph[in_node].append(node)
            graph[node]  # initialize a node in case it only has edges coming in
    return graph


def eulerian_cycle(graph):
    current_path = [choice(list(graph.keys()))]  # initialize a stack with a random node
    cycle = []  # list to store final cycle
    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining cycle
            cycle.append(current_path.pop())  # Remove the current node and put it in the cycle
    cycle = "->".join(cycle[::-1])  # print in reverse
    return cycle


# If the graph is not balanced, you may still be able find a Eulerian path - every edge is visited once, but the start
# and end nodes will be different. Eulerian path requires that at most one node has out degree - in degree == 1 (start),
# and at most one node has in degree - out degree == 1 (end), with all other nodes being balanced.
# ----------------------------------------------------------------------------------------------------------------------


def count_edges(graph):
    nodes = set(graph)
    edge_df = pd.DataFrame(
        {"Edges_in": [0] * len(nodes),
         "Edges_out": [0] * len(nodes)},
        index=nodes)
    for node in graph:
        edge_df.at[node, 'Edges_out'] = len(graph[node])
        for target in graph[node]:
            edge_df.at[target, 'Edges_in'] += 1
    return edge_df


def get_start_node(graph):
    count_df = count_edges(graph)
    start = None
    for node in count_df.index:
        if count_df.at[node, "Edges_out"] - count_df.at[node, "Edges_in"] == 1:
            start = node
    if start is None:
        raise ValueError("No start node was found")
    return start


def eulerian_path(graph):
    current_path = [get_start_node(graph)]  # initialize a stack with the starting node
    path = []  # list to store final path
    while current_path:
        current_node = current_path[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()  # Find and remove next node adjacent to current node
            current_path.append(next_node)  # Push the new node to the stack
        else:  # back-track to find remaining cycle
            path.append(current_path.pop())  # Remove the current node and put it in the cycle

    path = "->".join(path[::-1])  # print in reverse
    return path


# The text then states we now have a method to assemble a genome, as the string reconstruction problem reduces to
# finding a Eulerian path in the De Bruijn graph generated from reads
# ----------------------------------------------------------------------------------------------------------------------


def string_reconstruction(kmer_list):
    graph = debruijn_graph(kmer_list)
    path = eulerian_path(graph).split("->")
    return spell_string_by_path(path)


# Then notes that we are capable of constructing the k-universal string for any k, as this reduces to finding a Eulerian
# cycle in the De Bruijn Graph. In particular, interested in a circular string like a bacterial chromosome. Uses binary
# kmers as an example, where a k-universal string contains every binary kmer exactly once. Eg, for k=3 the string
# 01110100 contains all 8 binary kmers, mindful that you can loop back (011, 111, 110, 101, 010, 100, 000, 001)
# ----------------------------------------------------------------------------------------------------------------------


def binary_kmers(k):
    # https://stackoverflow.com/questions/64890117/
    def recursion(k, kmer="", kmer_list=[]):
        if len(kmer) == k:
            kmer_list.append(kmer)
        else:
            recursion(k, kmer + "0", kmer_list)
            recursion(k, kmer + "1", kmer_list)
        return kmer_list
    return recursion(k)


def k_circular_universal_string(k):
    kmer_list = binary_kmers(k)
    dbg = debruijn_graph(kmer_list)
    cycle = eulerian_cycle(dbg).split("->")
    path = spell_string_by_path(cycle)
    return path[:len(path) - k+1]


# Notes that assembling a random text is trivial, as random strings are not expected to have long repeats. And as read
# length becomes longer and longer, the DBG becomes less tangled - when reads become longer than repeats, it turns into
# a path. However, difficult to generate long and accurate reads. Read pairs then used as an improvement - two reads of
# length k, separated by unknown distance, d. A read pair of Read1 and Read2 corresponds to two edges in DBG(reads) -
# these two reads are separated by d, there must be a path of length k + d + 1 connecting the node at the beginning of
# Read1 with the node at the end of edge for Read2. **If** there is only one path connecting these nodes, or if all such
# paths spell out the same string, the read pair can be transformed into a virtual read of length 2 * k + d. In reality,
# repetitive regions complicate this scenario. Instead, analyze using paired composition: given a string text, a
# (k,d)-mer is a pair of kmers separated by d. (ATG|GGG) is a (3,4)-mer in TAATGCCATGGGATGTT.
# ----------------------------------------------------------------------------------------------------------------------


def paired_composition(text, k, d):
    kmer_list = []
    for i in range(len(text) - (2*k-d+1)):
        kmer1 = text[i:i+k]
        kmer2 = text[(i+d+k):(i+2*k+d)]
        kmer_list.append(kmer1 + "|" + kmer2)
    kmer_list.sort()
    return kmer_list


# It is noted that while the example text TAATGCCATGGGATGTT has 3 repeated 3mers, it has no repeated (3,1)mers in
# its paired composition. Additionally, the similar text TAAAGCCATGGGATGTT has the same 3mer composition, but different
# (3,1)mer paired compositions, allowing us to distinguish between them. The goal then becomes to adapt the DBG to
# reconstruct a string from its (k,d)mer paired composition. Given a (k,d)-mer (a1...ak | b1...bk) we have the following
# (k-1, d+1)mers: PREFIX = (a1...ak-1 | b1...bk-1) and SUFFIX = (a2...ak | b2...bk). For consecutive (k,d)mers in text,
# the suffix of the first (k,d)mer is equal to the prefix of the second (k,d)mer. Eg, for the consecutive (3,1)mers
# (TAA|GCC) and (AAT|CCA) in TAATGCCATGGGATGTT: SUFFIX((TAA|GCC)) = PREFIX((AAT|CCA)) = (AA|CC)). We can form a path
# graph of text that represents a path formed by len(text) - (k+d+k) + 1 edges corresponding to all (k,d)mers in text.
# Edges in this path are labeled by (k,d)mers
# ----------------------------------------------------------------------------------------------------------------------


def spell_edge(node_list):
    path = []
    for i in range(len(node_list) - 1):  # combine the (k-1) prefix/suffix to form the kmer edge
        path.append(node_list[i] + node_list[i+1][-1])
    return path


def path_graph(text, k, d):
    node_pair1 = []
    node_pair2 = []
    for i in range(len(text) - (2*k-d)):  # creates paired (k-1)mers prefix/suffix nodes
        node_pair1.append(text[i:i+k-d])  # read 1
        node_pair2.append(text[(i+d+k):(i+2*k+d-1)])  # read 2
    path1 = spell_edge(node_pair1)
    path2 = spell_edge(node_pair2)
    paired_path = list(map((lambda x, y: x + "|" + y), path1, path2))
    return " -> ".join(paired_path)


# The paired DBG is formed by gluing identical nodes in path graph. This results in a less tangled graph than DBG
# constructed from the individual reads of the paired reads.
# ----------------------------------------------------------------------------------------------------------------------


def paired_dbg_from_text(text, k, d):
    read1_patterns = [text[i:i+k] for i in range(len(text) - (2*k+d-1))]
    read2_patterns = [text[i:i+k] for i in range(k+d, len(text) - k+1)]
    graph = defaultdict(list)
    for i in range(len(read1_patterns)):
        read_pair_prefix = "|".join([read1_patterns[i][:-1], read2_patterns[i][:-1]])
        read_pair_suffix = "|".join([read1_patterns[i][1:], read2_patterns[i][1:]])
        graph[read_pair_prefix].append(read_pair_suffix)
        graph[read_pair_suffix]  # init suffix in case it has no outgoing edges
    return graph


# The text then defines composition graph(k, d, text) as the graph consisting of the len(text) - (k + d + k) + 1
# isolated edges labeled by the (k,d)mers in text - gluing identical nodes results in the same DBG as gluing identical
# nodes in the path graph. In practice we don't know text, but we can form the composition graph directly from the
# (k,d)mer composition of text, and from here construct a DBG and find a Eulerian path.
# ----------------------------------------------------------------------------------------------------------------------


def split_read_patterns(paired_list):
    read1_patterns = []
    read2_patterns = []
    for pattern in paired_list:
        read_split = pattern.split("|")
        read1_patterns.append(read_split[0])
        read2_patterns.append(read_split[1])
    return read1_patterns, read2_patterns


def paired_dbg_from_composition(paired_list):
    read1_patterns, read2_patterns = split_read_patterns(paired_list)
    graph = defaultdict(list)
    for i in range(len(read1_patterns)):
        read_pair_prefix = "|".join([read1_patterns[i][:-1], read2_patterns[i][:-1]])
        read_pair_suffix = "|".join([read1_patterns[i][1:], read2_patterns[i][1:]])
        graph[read_pair_prefix].append(read_pair_suffix)
        graph[read_pair_suffix]  # init suffix in case it has no outgoing edges
    return graph


# However, not every Eulerian path in the paired DBG spells out a solution to the String Reconstruction from Read-Pairs
# problem. Therefore, need to be able to generate all Eulerian cycles: The idea is to transform a single labeled
# directed graph containing n >= 1 cylces to n different simple directed graphs, each with a single Eulerian cycle.
# This transformation is easily invertible: given a unique Eulerian cycle in one of the simple graphs, we can find it in
# the original. Given a node v in Graph with indegree > 1 and with incoming edge (u, v) and outgoing edge (v, w),
# we construct a simpler (u, v, w) bypass graph where edges (u, v) and (v, w) are removed, and add a new node x with
# edges (u, x) and (v, x). These new edges inherit the labels for the removed edges: given an incoming edge (u, v) into
# v along with k outgoing edges (v, w1) ... (v, wk) from v, we can construct k different bypass graphs, and each will
# have a unique Eulerian cycle. Iteratively construct every bypass graph, keeping only those that are connected, which
# can be checked using breadth first search
# ----------------------------------------------------------------------------------------------------------------------


def bfs(graph):
    s = choice(list(graph.keys()))  # start with random node
    visited = [False] * len(graph.keys())  # array keeps track of what nodes we have already seen
    queue = deque([s])  # use a queue to track nodes that have yet to be visited, beginning with start
    visited[int(s)] = True
    while queue:
        s = queue.popleft()
        for i in graph[s]:
            if not visited[int(i)]:
                queue.append(i)
                visited[int(i)] = True
    return visited


def get_bypass_graphs(graph):
    all_graphs = [graph]
    for g in all_graphs:
        edges_in = count_edges(g)["Edges_in"]
        if all(edges_in < 2):
            continue
        in_gt1 = choice(edges_in[edges_in > 1].index.values)
        edges_in = [edge for edge in g if in_gt1 in g[edge]]
        edges_out = g[in_gt1]
        for node1 in edges_in:
            for node2 in edges_out:
                g_cp = g.copy()
                g_cp[node1] = [node2]
                if all(bfs(g_cp)):  # the more simple graph must be connected
                    all_graphs.append(g_cp)
        all_graphs.remove(g)
    return all_graphs


# Once we have all of the simple bypass graphs, we can find a Eulerian path for each, then spell the path from the
# ordered gapped paired reads (this can be imagined by lining up the reads as rows with gaps, and matching the columns
# to get the consensus base). However, not all Eulerian paths in the paired DBG will give a consensus string
# ----------------------------------------------------------------------------------------------------------------------


def spell_string_by_gapped_pattern(patterns, k, d):
    first_patterns = [kmer.split("|")[0] for kmer in patterns]
    second_patterns = [kmer.split("|")[1] for kmer in patterns]
    prefix_string = spell_string_by_path(first_patterns)
    suffix_string = spell_string_by_path(second_patterns)
    for i in range(k+d+1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i-k-d]:
            return "There is no consensus pattern"
    return prefix_string + suffix_string[len(suffix_string) - (k+d):]


def spell_all_eulerian_paths(graph_list, k, d):
    for g in graph_list:
        cycle = eulerian_path(g).split("->")
        path = spell_string_by_gapped_pattern(cycle, k, d)
        print(path)
    return None


def paired_read_string_reconstruction(input_text, k, d):
    graph = paired_dbg_from_composition(input_text)
    bpg = get_bypass_graphs(graph)
    spell_all_eulerian_paths(bpg, k, d)
    return None


# Defines coverage: if we have a kmer within the genome, coverage is the number of reads to which this kmer belongs. So
# far we have assumed perfect kmer coverage, which just doesn't happen. Eg, most Illumina reads are ~300nt, but still
# misses many 300mers in the genome, and nearly all reads have errors. Then shows example of read breaking: gives a
# 20mer and 4 10mers, which don't cover every 10mer in the 20mer. However, fragmenting the reads into 5mers provides
# perfect coverage of the 20mer. Trade-off: smaller k increases chance of perfect coverage, but results in a more
# tangled DBG which may be harder to infer Euler path. Even after read breaking, most assemblies still have gaps in kmer
# coverage, resulting in a DBG with missing edges, causing Euler path to fail. So rather than go for entire genome,
# settle for contigs, which can be found from the DBG. A path is non branching if it has node in and out = 1 for
# intermediate nodes of the path. Want the maximal of such paths: these strings of nucleotides must be present in any
# given assembly with a given kmer composition. Contigs are strings spelled by the maximal non branching paths in a DBG.
# Even if you have perfect coverage, rely on contigs because repeats prevent inferring a unique Eulerian path. Note
# that the text uses a collection of kmers as the input to the get_contigs problem
# ----------------------------------------------------------------------------------------------------------------------


def is_1in1out(count_df, node):
    return (count_df.loc[node, "Edges_in"] == 1) and (count_df.loc[node, "Edges_out"] == 1)


def get_contigs(graph):
    contigs = []
    count_df = count_edges(graph)
    for start_node in graph.keys():
        if not is_1in1out(count_df, start_node):  # iterate and start from all non 1-in-1-out nodes
            if count_df.loc[start_node, "Edges_out"] > 0:
                for next_node in graph[start_node]:
                    nb_path = [start_node, next_node]
                    while is_1in1out(count_df, next_node).all():  # all intermediate nodes are 1-in-1-out
                        next_node = graph[next_node][0]
                        nb_path.append(next_node)
                    contigs.append(nb_path)
    contigs = [spell_string_by_path(i) for i in contigs]  # glue series of nodes into paths
    return contigs


# The text notes that this implementation fails to account for self isolated cycles: cycles where the nodes also have
# in and out degree equal to one (eg, 6 -> 7 -> 6). The following implementation also includes these isolated cycles.
# Note that unlike the get_contigs problem, whose input is nucleotide kmers, the text phrases this problem as using
# sets of numeric edges of a graph as input
# ----------------------------------------------------------------------------------------------------------------------


def get_isolated_cycle(graph):
    all_cycles = []
    count_df = count_edges(graph)
    count_df["Visited"] = False  # init column tracking if node has been visited
    for start_node in graph.keys():
        if count_df.at[start_node, "Visited"]:  # skip if seen to prevent adding both n1->n2->n1 and n2->n1->n1
            continue
        count_df.at[start_node, "Visited"] = True
        if is_1in1out(count_df, start_node):
            next_node = graph[start_node][0]
            cycle = [start_node, next_node]
            while is_1in1out(count_df, next_node).all():
                count_df.at[next_node, "Visited"] = True
                next_node = graph[next_node][0]
                cycle.append(next_node)
                if next_node is start_node:
                    all_cycles.append(cycle)
                    break
    return all_cycles


def maximal_nonbranching_paths(graph):
    paths = []
    count_df = count_edges(graph)
    for start_node in graph.keys():
        if not is_1in1out(count_df, start_node):  # iterate and start from all non 1-in-1-out nodes
            if count_df.loc[start_node, "Edges_out"] > 0:
                for next_node in graph[start_node]:
                    nb_path = [start_node, next_node]
                    while is_1in1out(count_df, next_node).all():  # all intermediate nodes are 1-in-1-out
                        next_node = graph[next_node][0]
                        nb_path.append(next_node)
                    paths.append(nb_path)
    for cycle in get_isolated_cycle(graph):  # need to include isolated cycles, eg 6 -> 7 -> 6
        paths.append(cycle)
    return paths


# Example inputs
# ----------------------------------------------------------------------------------------------------------------------


graph_input = [
    "0 -> 3",
    "1 -> 0",
    "2 -> 1,6",
    "3 -> 2",
    "4 -> 2",
    "5 -> 4",
    "6 -> 5,8",
    "7 -> 9",
    "8 -> 7",
    "9 -> 6"
]

nbp_input = [
    "1 -> 2",
    "2 -> 3",
    "3 -> 4,5",
    "6 -> 7",
    "7 -> 6"
]

paired_input = [
"GAGA|TTGA",
"TCGT|GATG",
"CGTG|ATGT",
"TGGT|TGAG",
"GTGA|TGTT",
"GTGG|GTGA",
"TGAG|GTTG",
"GGTC|GAGA",
"GTCG|AGAT"
]

text_input = "TAATGCCATGGGATGTT"

contig_input = [
    "ATG",
    "ATG",
    "TGT",
    "TGG",
    "CAT",
    "GGA",
    "GAT",
    "AGA"
]
