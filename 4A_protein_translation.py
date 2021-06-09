# Implementing code to translate an RNA string into the corresponding protein

# Problem 4A in the BALA textbook/Rosalind

# Input is a text file where the first line is the RNA string to be translated
# AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

# Output is protein string
# MAMAPRTEINSTRING

# Usage: python3 4A_protein_translation.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_genetic_code(bases="UCAG", acids="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    """
    returns a dictionary mapping the 64 codons to the 20 amino acids + stop codons
    """
    # https://www.petercollingridge.co.uk/tutorials/bioinformatics/codon-table/
    codons = [a + b + c for a in bases for b in bases for c in bases]
    codon_table = dict(zip(codons, acids))
    return codon_table


def protein_translation(pattern, genetic_code):
    """
    pattern: an RNA string
    genetic_code: a dict mapping RNA codons to amino acids/stop codons
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


def main():
    """
    Read the input file, parse the arguments, and return the translated protein string
    """
    argv = list(sys.argv)
    input_text = ""

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text = line.replace('\n', '')

    code = get_genetic_code()
    output = protein_translation(input_text, code)
    print(output)


if __name__ == "__main__":
    main()
