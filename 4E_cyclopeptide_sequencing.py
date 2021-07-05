# Implementing XXX

# Problem 4E in the BALA textbook/Rosalind

# Input is XXX
# XX

# Output XXX
# XX

# Usage: python3 4E_cyclopeptide_sequencing.py input.txt > output.txt
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


def cyclic_spectrum_int(peptide_int, aa_mass):
    peptide = [int(x) for x in peptide_int.split("-")]
    prefix_mass = [0]  # init array of prefix masses with 0
    for i in range(1, len(peptide) + 1):  # get the prefix masses in increasing order, one aa at a time
        aa = int(peptide[i-1])
        mass = prefix_mass[i-1] + aa
        prefix_mass.append(mass)
    peptide_mass = prefix_mass[-1]  # full peptide mass is the last entry of prefix array
    cyclic_spectrum = [0]  # init the cyclic spectrum array with 0
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            cyclic_spectrum.append(mass)
            if i > 0 and j < len(peptide):  # account for subpeptide wrapping back around to beginning
                cyclic_spectrum.append(peptide_mass - mass)
    cyclic_spectrum.sort()
    return cyclic_spectrum


def linear_spectrum_int(peptide_int, aa_mass):
    peptide = [int(x) for x in peptide_int.split("-")]
    prefix_mass = [0]  # init array of prefix masses with 0
    for i in range(1, len(peptide) + 1):  # get the prefix masses in increasing order, one aa at a time
        aa = int(peptide[i - 1])
        mass = prefix_mass[i - 1] + aa
        prefix_mass.append(mass)
    linear_spectrum = [0]  # init the linear specturm array with 0
    for i in range(len(peptide)):  # get the mass of every subpeptide
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    linear_spectrum.sort()
    return linear_spectrum


def get_integer_mass(peptide_int):
    mass = [int(x) for x in peptide_int.split("-")]
    return sum(mass)


def grow_peptide(peptides, aa_mass):
    if len(peptides) == 0:
        new_peptides = set(aa_mass.values())
    else:
        new_peptides = set()
        for peptide in peptides:
            for aa in aa_mass.values():
                new_peptides.add("-".join([str(peptide), str(aa)]))
    return new_peptides


def is_peptide_consistent(peptide_int, spectrum, aa_mass):
    peptide_spectrum = linear_spectrum_int(peptide_int, aa_mass)
    spectrum_cp = spectrum.copy()
    for mass in peptide_spectrum.copy():
        if mass in spectrum_cp:
            peptide_spectrum.remove(mass)
            spectrum_cp.remove(mass)
    return len(peptide_spectrum) == 0


def cyclopeptide_sequencing(spectrum, aa_mass):
    output_peptide = set()
    peptides = set(spectrum).intersection(set(aa_mass.values()))
    parent_mass = spectrum[-1]
    while peptides:
        peptides = grow_peptide(peptides, aa_mass)
        for peptide in peptides.copy():
            if get_integer_mass(peptide) == parent_mass:
                if cyclic_spectrum_int(peptide, aa_mass) == spectrum:
                    output_peptide.add(peptide)
                    peptides.remove(peptide)
            elif is_peptide_consistent(peptide, spectrum, aa_mass) is False:
                peptides.remove(peptide)
    return output_peptide


def main():
    """
    Read the input file, parse the arguments, and XXX
    """
    argv = list(sys.argv)
    input_spectrum = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_spectrum is None:
                input_spectrum = str(line.replace('\n', ''))

    aa_mass = get_aa_mass()
    input_spectrum = [int(x) for x in input_spectrum.split(" ")]
    output = cyclopeptide_sequencing(input_spectrum, aa_mass)

    for peptide in output:
        print(peptide, end=" ")


if __name__ == "__main__":
    main()
