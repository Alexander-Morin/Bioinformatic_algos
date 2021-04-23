# Chapter 3 is focused on genome assembly

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
# Read1 with the node at the end of edge for Read2. If there is only one path connecting these nodes, or if all such
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


# It is noted that while the example text has TAATGCCATGGGATGTT has 3 repeated 3mers, it has no repeated (3,1)mers in
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
        read1 = text[i:i+k-d]
        read2 = text[(i+d+k):(i+2*k+d-1)]
        node_pair1.append(read1)
        node_pair2.append(read2)
    path1 = spell_edge(node_pair1)
    path2 = spell_edge(node_pair2)
    paired_path = list(map((lambda x, y: x + "|" + y), path1, path2))
    return " -> ".join(paired_path)


# The paired DBG is formed by gluing identical nodes in path graph. This results in a less tangled graph than DBG
# constructed from individual reads. The text then defines composition graph(k, d, text) as the graph consisting of the
# len(text) - (k + d + k) + 1 isolated edges labeled by the (k,d)mers in text - gluing identical nodes results in the
# same DBG as gluing identical nodes in the path graph. In practice we don't know text, but we can form the composition
# graph directly from the (k,d)mer composition of text, and from here construct a DBG and find a Eulerian path. However,
# not every Eulerian path in the paired DBG spells out a solution to the String Reconstruction from Read-Pairs problem.
# Therefore, need to be able to generate all Eulerian cycles: The idea is to transform a single labeled directed graph
# containing n >= 1 cylces to n different simple directed graphs, each with a single Eulerian cycle. This transformation
# is easily invertible: given a unique Eulerian cycle in one of the simple graphs, we can find it in the original.
# Given a node v in Graph with indegree > 1 and with incoming edge (u, v) and outgoing edge (v, w), we construct a
# simpler (u, v, w) bypass graph where edges (u, v) and (v, w) are removed, and add a new node x with edges (u, x) and
# (v, x). These new edges inherit the labels for the removed edges: given an incoming edge (u, v) into v along with k
# outgoing edges (v, w1) ... (v, wk) from v, we can construct k different bypass graphs, and each will have a unique
# Eulerian cycle. So we iteratively construct every bypass graph until we have a large family of simple graphs
# ----------------------------------------------------------------------------------------------------------------------

# TODO: not all bypass graphs get added. eg for node 2 each of 4 gets added, but terminates before getting to 6

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
        print(all_graphs)
        edges_in = count_edges(g)["Edges_in"]
        if all(edges_in == 1):
            continue
        in_gt1 = choice(edges_in[edges_in > 1].index.values)
        edges_in = [edge for edge in g if in_gt1 in g[edge]]
        edges_out = g[in_gt1]
        for edge1 in edges_in:
            for edge2 in edges_out:
                g_cp = g.copy()
                g_cp[edge1] = [edge2]
                all_graphs.append(g_cp)
            all_graphs.remove(g)
    return all_graphs


def spell_string_by_gapped_pattern(patterns, k, d):
    first_patterns = [kmer.split("|")[0] for kmer in patterns]
    second_patterns = [kmer.split("|")[1] for kmer in patterns]
    prefix_string = spell_string_by_path(first_patterns)
    suffix_string = spell_string_by_path(second_patterns)
    for i in range(k+d+1, len(prefix_string)):
        if prefix_string[i] != suffix_string[i-k-d]:
            return "There is no consensus pattern"
    return prefix_string + suffix_string[len(suffix_string) - (k+d):]



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



simple_graph_input = [
    "0 -> 1",
    "1 -> 2",
    "2 -> 3",
    "3 -> 0"
]

path_input = [
"GACC|GCGC",
"ACCG|CGCC",
"CCGA|GCCG",
"CGAG|CCGG",
"GAGC|CGGA"
]

graph = parse_graph(graph_input)
# graph = parse_graph(simple_graph_input)
# print(count_edges(graph)["Edges_in"])
# print(graph)
# print(bfs(graph))
# tt = get_bypass_graphs(graph)
# for i in tt:
#     print(i)
dd = spell_string_by_gapped_pattern(path_input, 4, 2)
print(dd)


