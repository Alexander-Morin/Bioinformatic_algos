# Implementing code to find the theoretical spectrum of a cyclic peptide: the collection of integer masses for all of
# its subpeptides, including 0 and the total mass. Builds off of finding a linear spectrum (which does not consider
# the subpeptides that can loop back to the beginning).

# NQEL subpeptides -> [-, N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, LNQ, NQEL]
# NQEL prefixes -> [-, N, NQ, NQE, NQEL]

# Problem 4C in the BALA textbook/Rosalind

# Input is a text file where the first line corresponds to the peptide
# NEQL

# Output is the integer masses of the linear spectrum
# 0 113 114 128 129 227 242 242 257 355 356 370 371 484

# Usage: python3 4C_cyclic_spectrum.py input.txt > output.txt
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


def main():
    """
    Read the input file, parse the arguments, and print the integer masses of the cyclic spectrum of the peptide
    """
    argv = list(sys.argv)
    input_peptide_string = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_peptide_string is None:
                input_peptide_string = str(line.replace('\n', ''))

    aa_mass = get_aa_mass()
    output = cyclic_spectrum(input_peptide_string, aa_mass)

    for mass in output:
        print(mass)


if __name__ == "__main__":
    main()
