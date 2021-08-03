# Implementing code to select the highest scoring peptides whose linear spectrums most closely align with a given
# experimental spectrum. This is used for a branch and bound approach for cyclopeptide sequencing - it is expected that
# mass spec will result in false or missing masses. Therefore, use a scoring function to determine how close candidate
# theoretical peptides are to the given spectrum, only taking n candidates.

# Problem 4L in the BALA textbook/Rosalind

# Input is a text file where the first line is the leaderboard of peptides, the second line is the experimental
# spectrum, and the third line is n, the number of top candidates to keep
# LAST ALST TLLT TQAS
# 0 71 87 101 113 158 184 188 259 271 372
# 2

# Output is the top n candidates to keep
# LAST ALST

# Usage: python3 4L_trim_leaderboard.py input.txt > output.txt
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
    linear_spectrum = [0]  # init the linear spectrum array with 0
    for i in range(len(peptide)):  # get the mass of every subpeptide
        for j in range(i+1, len(peptide)+1):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    linear_spectrum.sort()
    return linear_spectrum


def linear_peptide_score(peptide, spectrum, aa_mass):
    """
    peptide: a peptide string
    spectrum: a list of peptide inter masses
    aa_mass: dict mapping amino acids to their integer masses
    returns an integer score of how many masses in the linear spectrum of peptide are contained in the input spectrum
    """
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = linear_spectrum(peptide, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def trim(leaderboard, spectrum, n, aa_mass):
    """
    leaderboard: a list of candidate peptides
    spectrum: a list of peptide inter masses
    n: an int corresponding to how many top elements to keep
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of the top n peptides
    """
    scores = [linear_peptide_score(peptide, spectrum, aa_mass) for peptide in leaderboard]
    score_df = pd.DataFrame({"Score": scores}, index=leaderboard)
    score_df = score_df.sort_values("Score", ascending=False)
    n_score = score_df["Score"].iloc[n]  # the score of the nth element to keep
    leading_peptides = score_df.index[score_df["Score"] >= n_score]
    return leading_peptides.tolist()


def main():
    """
    Read the input file, parse the arguments, and returns top n peptides
    """
    argv = list(sys.argv)
    input_leaderboard = None
    input_spectrum = None
    input_topn = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_leaderboard is None:
                input_leaderboard = line.replace('\n', '').split(" ")
            elif input_spectrum is None:
                input_spectrum = line.replace('\n', '')
            elif input_topn is None:
                input_topn = int(line.replace('\n', ''))

    aa_mass = get_aa_mass()
    input_spectrum = [int(x) for x in input_spectrum.split(" ")]
    top_peptides = trim(input_leaderboard, input_spectrum, input_topn, aa_mass)

    for peptide in top_peptides:
        print(peptide, end=" ")


if __name__ == "__main__":
    main()
