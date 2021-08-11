# Implementing code to find the lead seequenced peptide from an input spectrum. Leaderboard sequencing keep 'n' best
# scoring candidates, but struggles when there are more frequent mass spec errors. Additionally, need to consider amino
# acids beyond the standard 20 for the topic problem (sequencing a bacterial antibiotic non ribosomal peptide), so
# generate the spectral convolution (positive differences between subpeptides) of the input spectrum to get the 'm'
# most likely amino aicds. These are then fed into leaderboard seq.

# Problem 4I in the BALA textbook/Rosalind

# Input is a text file where the first line is the leaderboard threshold, second line is convolution threshold, and
# the third line is the input spectrum
# 20
# 60
# 57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493

# Output is lead peptiden in integer mass form
# 99-71-137-57-72-57

# Usage: python3 4I_convolution_cyclopeptide_sequencing.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import numpy as np
from collections import Counter
from itertools import repeat, chain


def aa_intmass_dict(intmass_list):
    """
    dict mapping a list of amino acid integer masses
    """
    aa_dict = {}
    for i in range(len(intmass_list)):
        aa_dict[i] = intmass_list[i]
    return aa_dict


def linear_spectrum_int(peptide_int):
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


def linear_peptide_score_int(peptide_int, spectrum):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
    spectrum: a list of peptide inter masses
    aa_mass: dict mapping amino acids to their integer masses
    returns an integer score of how many masses in the linear spectrum of peptide are contained in the input spectrum
    """
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = linear_spectrum_int(peptide_int)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def cyclic_spectrum_int(peptide_int):
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


def cyclopeptide_scoring(peptide_int, spectrum):
    """
    peptide_int: a peptide string represented by its integer masses: 129-227-257
    spectrum: a list of peptide inter masses
    aa_mass: dict mapping amino acids to their integer masses
    returns an integer score of how many masses in the cyclic spectrum of peptide are contained in the input spectrum
    """
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = cyclic_spectrum_int(peptide_int)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def trim_int(leaderboard, spectrum, leaderboard_topn):
    """
    leaderboard: a set of candidate peptides
    spectrum: a list of peptide integer masses
    leaderboard_topn: an int corresponding to how many top elements to keep
    aa_mass: dict mapping amino acids to their integer masses
    returns a list of the top n peptides
    """
    n = min(leaderboard_topn, len(leaderboard))
    leaderboard = list(leaderboard)
    scores = [linear_peptide_score_int(peptide, spectrum) for peptide in leaderboard]
    sorted_ix = np.argsort(scores)[::-1]
    sorted_scores = [scores[ix] for ix in sorted_ix]
    sorted_peptides = [leaderboard[ix] for ix in sorted_ix]
    score_threshold = sorted_scores[n-1]
    for i in range(n, len(leaderboard)):
        if sorted_scores[i] < score_threshold:
            return sorted_peptides[:i]
    return sorted_peptides


def leaderboard_cyclopeptide_seq(spectrum, leaderboard_topn, aa_mass):
    """
    spectrum: a list of the integer masses of a peptide spectrum
    leaderboard_topn: an int corresponding to how many top elements to keep
    aa_mass: dict mapping amino acids to their integer masses
    returns a string corresponding to the integer peptide with the best score to the spectrum
    """
    leaderboard = set(spectrum).intersection(set(aa_mass.values()))  # init with aa in spectrum
    lead_peptide, hi_score = "", 0
    parent_mass = max(spectrum)
    while leaderboard:
        leaderboard = grow_peptide(leaderboard, aa_mass)  # branch step
        for peptide in leaderboard.copy():
            if get_integer_mass(peptide) == parent_mass:
                peptide_score = cyclopeptide_scoring(peptide, spectrum)
                if peptide_score > hi_score:
                    lead_peptide = peptide
                    hi_score = peptide_score
            elif get_integer_mass(peptide) > parent_mass:
                leaderboard.remove(peptide)
        if len(leaderboard) > 0:   # bound step
            leaderboard = trim_int(leaderboard, spectrum, leaderboard_topn)
    return lead_peptide


def spectral_convolution(spectrum):
    """
    spectrum: list of integers corresponding to the input spectrum
    returns a list of the positive differences in subpeptides of the spectrum, sorted by element occurrence
    """
    conv = []
    for i in range(0, len(spectrum)):
        for j in range(0, len(spectrum)):
            if i == j:
                next
            elif spectrum[j] - spectrum[i] > 0:
                conv.append(spectrum[j] - spectrum[i])
    conv = list(chain.from_iterable(repeat(i, c) for i, c in Counter(conv).most_common()))
    return conv


def keep_topn_convolution(conv_topn, spectrum):
    """
    conv_topn: integer of how many of the most frequently occurring subpeptides in spectrum to keep
    spectrum: list of integers corresponding to the input spectrum
    """
    conv = np.array(spectral_convolution(spectrum))
    conv = conv[(conv >= 57) & (conv <= 200)]
    peptides, counts = np.unique(conv, return_counts=True)
    sort_ix = np.argsort(counts)[::-1]
    sorted_peptides = [peptides[ix] for ix in sort_ix]
    sorted_counts = [counts[ix] for ix in sort_ix]
    for i in range(conv_topn, len(sorted_counts)):
        if sorted_counts[i] < sorted_counts[conv_topn - 1]:
            return sorted_peptides[:i]
    return sorted_peptides


def convolution_cyclopeptide_seq(conv_topn, leaderboard_topn, spectrum):
    """
    conv_topn: integer of how many of the most frequently occurring subpeptides in spectrum to keep
    leaderboard_topn: an int corresponding to how many top scoring elements in leaderboard to keep
    spectrum: list of integers corresponding to the input spectrum
    returns the peptide (represented in integer mass) with the highest score
    """
    leaderboard = keep_topn_convolution(conv_topn, spectrum)
    aa_mass = aa_intmass_dict(leaderboard)
    lead_peptide = leaderboard_cyclopeptide_seq(spectrum, leaderboard_topn, aa_mass)
    return lead_peptide


def main():
    """
    Read the input file, parse the arguments, and find the lead peptide
    """
    argv = list(sys.argv)
    input_conv_topn = None
    input_leader_topn = None
    input_spectrum = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_conv_topn is None:
                input_conv_topn = int(line.replace('\n', ''))
            elif input_leader_topn is None:
                input_leader_topn = int(line.replace('\n', ''))
            elif input_spectrum is None:
                input_spectrum = line.replace('\n', '')

    input_spectrum = [int(x) for x in input_spectrum.split(" ")]
    output = convolution_cyclopeptide_seq(input_conv_topn, input_leader_topn, input_spectrum)
    print(output)


if __name__ == "__main__":
    main()
