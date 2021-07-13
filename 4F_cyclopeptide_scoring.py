# Implementing code to find the score between of the theoretical cyclic spectrum of an input peptide and an input
# experimental spectrum. Because mass spec is noisy (missing or false masses), need to generate a scoring function
# to determine the alignment between a theoretical and experimental spectrum.

# Problem 4F in the BALA textbook/Rosalind

# Input is a text file where the first line is the input peptide, and the second line is the input spectrum
# NQEL
# 0 99 113 114 128 227 257 299 355 356 370 371 484

# Output is the integer score
# 11

# Usage: python3 4F_cyclopeptide_scoring.py input.txt > output.txt
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


def cyclic_spectrum(peptide, aa_mass):
    """
    peptide: a peptide string
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of the cyclic spectrum of the peptide: the masses of the subpeptides, as well as 0 and the full
    length peptide.

    builds off of the linear spectrum, which assumes that the the mass of any subpeptide can be found by subtracting the
    masses of two prefixes. in addition, must find the subpeptides that wrap around to the beginning - each of these
    subpeptides has mass equal to the difference between the full peptide and a subpeptide in the linear spectrum.
    so mass(LN) = mass(NQEL) - mass(QE)

    generate an array of all prefix masses in increasing order, then the mass of subpeptide beginning at position i + 1
    and ending at j can be found by prefix_mass[j] - prefix_mass[i]
    """
    prefix_mass = [0]  # init array of prefix masses with 0
    for i in range(1, len(peptide) + 1):  # get the prefix masses in increasing order, one aa at a time
        aa = peptide[i-1]
        mass = prefix_mass[i-1] + aa_mass[aa]
        prefix_mass.append(mass)
    peptide_mass = prefix_mass[len(prefix_mass)-1]  # full peptide mass is the last entry of prefix array
    cyclic_spectrum = [0]  # init the cyclic spectrum array with 0
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            cyclic_spectrum.append(mass)
            if i > 0 and j < len(peptide):  # account for subpeptide wrapping back around to beginning
                cyclic_spectrum.append(peptide_mass - mass)
    cyclic_spectrum.sort()
    return cyclic_spectrum


def cyclopeptide_scoring(peptide, spectrum, aa_mass):
    """
    peptide: a peptide string
    spectrum: a list of peptide inter masses
    aa_mass: dict mapping amino acids to their integer masses
    returns an integer score of how many masses in the cyclic spectrum of peptide are contained in the input spectrum
    """
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = cyclic_spectrum(peptide, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def main():
    """
    Read the input file, parse the arguments, and returns the cyclopeptide score
    """
    argv = list(sys.argv)
    input_peptide = None
    input_spectrum = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_peptide is None:
                input_peptide = line.replace('\n', '')
            elif input_spectrum is None:
                input_spectrum = line.replace('\n', '')

    aa_mass = get_aa_mass()
    input_spectrum = [int(x) for x in input_spectrum.split(" ")]
    print(cyclopeptide_scoring(input_peptide, input_spectrum, aa_mass))


if __name__ == "__main__":
    main()
