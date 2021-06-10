# This chapter explores translation/encoding proteins, motivated by antibiotics discovery

import pandas as pd


# Text begins by motivating how Bacillus brevis creates/encodes the antibiotic Tyrocidine B1, which is 10 aa long
# ----------------------------------------------------------------------------------------------------------------------


def get_genetic_code(string_type="DNA"):
    # https://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/
    if string_type is "RNA":
        bases = "UCAG"
    elif string_type is "DNA":
        bases = "TCAG"
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    codon_table = dict(zip(codons, amino_acids))
    return codon_table


def protein_translation(pattern, genetic_code):
    protein = []
    for i in range(0, len(pattern) - 3, 3):  # step of 3 to represent codon
        codon = pattern[i:i+3]
        acid = genetic_code[codon]
        if acid is "*":  # stop codon
            break
        protein.append(acid)
    return "".join(protein)


# Next, look at finding which DNA strings encode for a given peptide string. Need to consider all three open reading
# frames as well as both strands
# ----------------------------------------------------------------------------------------------------------------------


def get_reverse_complement(dna_string):
    complement = []
    for letter in dna_string:
        if letter == 'A':
            complement.append('T')
        elif letter == 'C':
            complement.append('G')
        elif letter == 'T':
            complement.append('A')
        elif letter == 'G':
            complement.append('C')
    return "".join(complement[::-1])



def peptide_encoding(dna_string, peptide_string="s"):
    substrings = []
    rev_comp = get_reverse_complement(dna_string)
    genetic_code = get_genetic_code(string_type="DNA")

    frame1 = protein_translation(dna_string, genetic_code)
    frame2 = protein_translation(dna_string[1:], genetic_code)
    frame3 = protein_translation(dna_string[2:], genetic_code)
    rev_frame1 = protein_translation(dna_string, genetic_code)
    return frame2







# Example inputs
# ----------------------------------------------------------------------------------------------------------------------


genetic_code = get_genetic_code(string_type="RNA")
input_text = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
print(peptide_encoding(input_text))
