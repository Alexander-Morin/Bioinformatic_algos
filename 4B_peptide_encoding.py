# Implementing code to find which substrings of a DNA string encode a provided peptide.

# Problem 4B in the BALA textbook/Rosalind

# Input is a text file where the first line is the DNA string, the second line the peptide string
# ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
# MA

# Output
# ATGGCC
# GGCCAT
# ATGGCC

# Usage: python3 4B_peptide_encoding.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_genetic_code(string_type="DNA"):
    """
    returns a dictionary mapping the 64 codons to the 20 amino acids + stop codons
    """
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
    """
    pattern: an RNA string
    genetic_code: a dict mapping codons to amino acids/stop codons
    returns the translated protein string
    """
    protein = []
    for i in range(0, len(pattern), 3):  # step of 3 to represent codon
        codon = pattern[i:i+3]
        acid = genetic_code[codon]
        if acid is "*":  # stop codon
            break
        protein.append(acid)
    return "".join(protein)


def get_reverse_complement(dna_string):
    """
    dna_string: a string of A/C/G/Ts
    returns the reverse complement of dna string (ACTG -> CAGT)
    """
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
    """
    dna_string: a string of A/C/G/Ts
    peptide_string: a string of the single letter abbreviation of amino acids
    returns a list of the substrings of dna_string that encode for the provided peptide
    """
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


def main():
    """
    Read the input file, parse the arguments, find the substrings and print each on its own line
    """
    argv = list(sys.argv)
    input_dna_string = None
    input_peptide_string = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_dna_string is None:
                input_dna_string = str(line.replace('\n', ''))
            else:
                input_peptide_string = str(line.replace('\n', ''))

    output = peptide_encoding(input_dna_string, input_peptide_string)
    for sequence in output:
        print(sequence)


if __name__ == "__main__":
    main()
