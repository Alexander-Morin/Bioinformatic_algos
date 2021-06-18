# Implementing code to find the theoretical linear spectrum of a peptide: the integer masses of all linear subpeptides
# of a peptide, as well as 0 and the full length mass. Relies on the fact that the mass of linear subpeptides can be
# found by subtracting the masses of two prefixes. Eg,
# NQEL subpeptides -> [-, N, Q, E, L, NQ, QE, EL, NQE, QEL, NQEL]
# NQEL prefixes -> [-, N, NQ, NQE, NQEL]

# Problem 4J in the BALA textbook/Rosalind

# Input is a text file where the first line corresponds to the peptide
# NQEL

# Output is the integer masses of the linear spectrum
# 0 113 114 128 129 242 242 257 370 371 484

# Usage: python3 4J_linear_spectrum.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def get_aa_mass():
    """
    returns a dict mapping the single amino acid names to their integer mass
    """
    amino_acids = "GASPVTCILNDKQEMHFRYW"
    mass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    mass_dict = dict(zip(amino_acids, mass))
    return mass_dict


def linear_spectrum(peptide, aa_mass):
    """
    peptide: a peptide string
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of the linear spectrum of the peptide: the masses of the prefix peptides, as well as 0 and the full
    length peptide.

    based on the assumption that the mass of any subpeptide can be found by subtracting the masses of two prefixes.
    generate an array of all prefix masses in increasing order, then the mass of subpeptide beginning at position i + 1
    and ending at j can be found by prefix_mass[j] - prefix_mass[i]
    """
    prefix_mass = [0]  # init array of prefix masses with 0
    for i in range(1, len(peptide) + 1):  # get the prefix masses in increasing order, one aa at a time
        aa = peptide[i-1]
        mass = prefix_mass[i-1] + aa_mass[aa]
        prefix_mass.append(mass)
    linear_spectrum = [0]  # init the linear specturm array with 0
    for i in range(len(peptide)):  # get the mass of every subpeptide
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    linear_spectrum.sort()
    return linear_spectrum


def main():
    """
    Read the input file, parse the arguments, and print the integer masses of the linear spectrum of the peptide
    """
    argv = list(sys.argv)
    input_peptide_string = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_peptide_string is None:
                input_peptide_string = str(line.replace('\n', ''))

    aa_mass = get_aa_mass()
    output = linear_spectrum(input_peptide_string, aa_mass)

    for mass in output:
        print(mass)


if __name__ == "__main__":
    main()
