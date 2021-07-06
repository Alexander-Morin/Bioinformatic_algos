# Implementing code to find cyclical peptides (represented by the integer mass of the amino acids) that match an input
# experimental spectrum. Cyclopeptide sequencing is a branch and bound algorithm. The branching step grows candidate
# linear peptides by one amino acid at a time. The bound step then checks if the spectrum of these linear peptides
# is consistent with the experimental spectrum. Consistency means that every mass within the linear theoretical
# spectrum is contained within the experimental spectrum. If true (and the mass of the theoretical spectrum matches
# the total mass of the experimental spectrum), the cyclical spectrum is checked to see if it is a match.

# Problem 4E in the BALA textbook/Rosalind

# Input is a text file where the first line is the experimental spectrum
# 0 113 128 186 241 299 314 427

# Output is the candidate peptides
# 186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186

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


def linear_spectrum_int(peptide_int, aa_mass):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of the linear spectrum of the peptide: the masses of the prefix peptides, as well as 0 and the full
    length peptide.

    based on the assumption that the mass of any subpeptide can be found by subtracting the masses of two prefixes.
    generate an array of all prefix masses in increasing order, then the mass of subpeptide beginning at position i + 1
    and ending at j can be found by prefix_mass[j] - prefix_mass[i]
    """
    peptide = [int(x) for x in peptide_int.split("-")]
    prefix_mass = [0]  # init array of prefix masses with 0
    for i in range(1, len(peptide) + 1):  # get the prefix masses in increasing order, one aa at a time
        aa = int(peptide[i - 1])
        mass = prefix_mass[i - 1] + aa
        prefix_mass.append(mass)
    linear_spectrum = [0]  # init the linear spectrum array with 0
    for i in range(len(peptide)):  # get the mass of every subpeptide
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    linear_spectrum.sort()
    return linear_spectrum


def cyclic_spectrum_int(peptide_int, aa_mass):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
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


def get_integer_mass(peptide_int):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
    returns the sum of the integer masses
    """
    mass = [int(x) for x in peptide_int.split("-")]
    return sum(mass)


def grow_peptide(peptides, aa_mass):
    """
    peptides: a list of peptides represented by their integer masses: [129-227-257, 113-114-128]
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of peptides, where each peptide from the original list has been expanded by 1 of each of all
    of the aa integer masses
    """
    if len(peptides) == 0:
        new_peptides = set(aa_mass.values())
    else:
        new_peptides = set()
        for peptide in peptides:
            for aa in aa_mass.values():
                new_peptides.add("-".join([str(peptide), str(aa)]))
    return new_peptides


def is_peptide_consistent(peptide_int, spectrum, aa_mass):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
    spectrum: a list of the integer masses of a peptide spectrum
    aa_mass: dict mapping amino acids to their integer masses
    returns True/False whether or not all of the masses found in the linear spectrum of the peptide are found in the
    input spectrum
    """
    peptide_spectrum = linear_spectrum_int(peptide_int, aa_mass)
    spectrum_cp = spectrum.copy()
    for mass in peptide_spectrum.copy():
        if mass in spectrum_cp:
            peptide_spectrum.remove(mass)
            spectrum_cp.remove(mass)
    return len(peptide_spectrum) == 0


def cyclopeptide_sequencing(spectrum, aa_mass):
    """
    spectrum: a list of the integer masses of a peptide spectrum
    aa_mass: dict mapping amino acids to their integer masses
    returns a set of the candidate peptides whose cyclic spectrum matches the input spectrum
    """
    output_peptide = set()
    peptides = set(spectrum).intersection(set(aa_mass.values()))  # init with aa in spectrum
    parent_mass = max(spectrum)
    while peptides:
        peptides = grow_peptide(peptides, aa_mass)  # branch step
        for peptide in peptides.copy():
            if get_integer_mass(peptide) == parent_mass:
                if cyclic_spectrum_int(peptide, aa_mass) == spectrum:
                    output_peptide.add(peptide)
                    peptides.remove(peptide)
            elif not is_peptide_consistent(peptide, spectrum, aa_mass):  # bound step
                peptides.remove(peptide)
    return output_peptide


def main():
    """
    Read the input file, parse the arguments, and print the candidate peptides
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
