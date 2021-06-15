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
    for i in range(0, len(pattern), 3):  # step of 3 to represent codon
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


def peptide_encoding(dna_string, peptide_string):
    substrings = []
    rev_comp = get_reverse_complement(dna_string)
    genetic_code = get_genetic_code(string_type="DNA")

    for i in range(len(dna_string) - len(peptide_string) * 3):  # slide window over bases size of codons in peptide
        sequence = dna_string[i:i+len(peptide_string)*3]
        peptide = protein_translation(sequence, genetic_code)
        rev_sequence = rev_comp[i:i + len(peptide_string) * 3]
        rev_peptide = protein_translation(rev_sequence, genetic_code)

        if peptide == peptide_string:
            substrings.append(sequence)
        if rev_peptide == peptide_string:
            substrings.append(get_reverse_complement(rev_sequence))  # revert back to original direction

    return substrings


# However, the text notes that applying this process to find the Tyrocidine B1 substring in the Bacillus brevis genome
# does not result in an expected 30mer. Goes on to explain that tyrocidines and gramicidins are cyclic peptides - yet
# if you were to attempt to repeat the process over the 10 linear strings of this cyclic peptide, you would still not
# find an answer. This is because these proteins are not translated by ribosomes, but instead are made by non ribosomal
# peptide (RBP) synthetase, which pieces together antibiotic peptides without any reliance on the genetic code. Notes
# the general usefulness of these NRPs. Therefore sequence these peptides using mass spectrometry. Approximate mass with
# the integer mass (notes in reality we would be using non-integer masses). Goal is to generate the experimental
# spectrum of masses after mass spec breaks down the protein. Also make simplifying assumption that the cyclic peptide
# is broken at every 2 bonds, so the resulting experimental spectrum contains all of the linear fragments of the peptide
# (called subpeptides). So NQEL -> [N,Q,E,L,NQ,QE,EL,LN,NQE,QEL,ELN,LNQ]. The theoretical spectrum of a peptide would
# be the collection of masses for all subpeptides, as well as 0, and the full peptide. First, generate linear spectrum
# ----------------------------------------------------------------------------------------------------------------------


def get_aa_mass():
    amino_acids = "GASPVTCILNDKQEMHFRYW"
    mass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    mass_dict = dict(zip(amino_acids, mass))
    return mass_dict


def linear_spectrum(peptide, aa_mass):
    prefix_mass = [0]
    for i in range(1, len(peptide) + 1):
        aa = peptide[i-1]
        mass = prefix_mass[i-1] + aa_mass[aa]
        prefix_mass.append(mass)
    linear_spectrum = [0]
    for i in range(len(peptide) - 1):
        for j in range(i + 1, len(peptide)):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    # linear_spectrum.sort
    return linear_spectrum
    # return prefix_mass




# Example inputs
# ----------------------------------------------------------------------------------------------------------------------


genetic_code = get_genetic_code(string_type="DNA")
aa_mass = get_aa_mass()
input_text = "NQEL"
print(linear_spectrum(input_text, aa_mass))