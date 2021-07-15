# Implementing code to XX

# Problem 4G in the BALA textbook/Rosalind

# Input is XX
# XX

# Output is XX
# XX

# Usage: python3 4G_leaderboard_cyclopeptide_sequencing.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import pandas as pd


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


def linear_peptide_score_int(peptide_int, spectrum, aa_mass):
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = linear_spectrum_int(peptide_int, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def trim_int(leaderboard, spectrum, n, aa_mass):
    scores = [linear_peptide_score_int(peptide, spectrum, aa_mass) for peptide in leaderboard]
    score_df = pd.DataFrame({"Score": scores}, index=leaderboard)
    score_df = score_df.sort_values("Score", ascending=False)
    # n_score = score_df["Score"].iloc[n]  # the score of the nth element to keep
    # leading_peptides = score_df.index[score_df["Score"] >= n_score]
    return score_df
    # return leading_peptides.tolist()


def leaderboard_cyclopeptide_seq(spectrum, n, aa_mass):
    output_peptide = set()
    leaderboard = set(spectrum).intersection(set(aa_mass.values()))  # init with aa in spectrum
    hi_score = 0
    parent_mass = max(spectrum)
    while leaderboard:
        leaderboard = grow_peptide(leaderboard, aa_mass)  # branch step
        for peptide in leaderboard.copy():
            if get_integer_mass(peptide) == parent_mass:
                peptide_score = linear_peptide_score_int(peptide, spectrum, aa_mass)
                if peptide_score > hi_score:
                    output_peptide.add(peptide)
                    hi_score = peptide_score
            elif get_integer_mass(peptide) > parent_mass:
                leaderboard.remove(peptide)
        df = trim_int(leaderboard, spectrum, n, aa_mass)
        print(df)
        # if len(leaderboard) > 0:
        #     trim_int(leaderboard, spectrum, n, aa_mass)
    return output_peptide


def main():
    """
    Read the input file, parse the arguments, and print the candidate peptides
    """
    argv = list(sys.argv)
    input_n = None
    input_spectrum = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_n is None:
                input_n = int(line.replace('\n', ''))
            elif input_spectrum is None:
                input_spectrum = line.replace('\n', '')

    aa_mass = get_aa_mass()
    input_spectrum = [int(x) for x in input_spectrum.split(" ")]
    output = leaderboard_cyclopeptide_seq(input_spectrum, input_n, aa_mass)

    for peptide in output:
        print(peptide, end=" ")


if __name__ == "__main__":
    main()
