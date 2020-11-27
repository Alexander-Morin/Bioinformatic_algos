# Chapter 3 is focused on genome assembly

import ch1_overview as ch1
import numpy as np

# First, focus on the more simplistic case where are all reads are kmers of length k, all reads come from the same
# strand, there are no errors, and they exhibit perfect coverage (every kmer substring of the genome has a read)

# Framed as string reconstruction using a walk in an overlap graph
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




pattern = [
    "ACCGA",
    "CCGAA",
    "CGAAG",
    "GAAGC",
    "AAGCT"
    ]

# ACCGAAGCT


print(spell_string_by_path(pattern))