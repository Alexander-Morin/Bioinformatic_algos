# This chapter explores translation/encoding proteins, motivated by antibiotics discovery

import pandas as pd
import numpy as np
from collections import Counter
from itertools import repeat, chain

# Text begins by motivating how Bacillus brevis creates/encodes the antibiotic Tyrocidine B1, which is 10 aa long
# ----------------------------------------------------------------------------------------------------------------------


def get_genetic_code(string_type="DNA"):
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
    protein = []
    for i in range(0, len(pattern), 3):  # step of 3 to represent codon
        codon = pattern[i:i+3]
        acid = genetic_code[codon]
        if acid is "*":  # stop codon
            break
        protein.append(acid)
    return "".join(protein)


# Next, look at finding which DNA strings encode for a given peptide string. Need to consider all three open reading
# frames as well as both strands
# ----------------------------------------------------------------------------------------------------------------------


def get_reverse_complement(dna_string):
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


# However, the text notes that applying this process to find the Tyrocidine B1 substring in the Bacillus brevis genome
# does not result in an expected 30mer. Goes on to explain that tyrocidines and gramicidins are cyclic peptides - yet
# if you were to attempt to repeat the process over the 10 linear strings of this cyclic peptide, you would still not
# find an answer. This is because these proteins are not translated by ribosomes, but instead are made by non ribosomal
# peptide (RBP) synthetase, which pieces together antibiotic peptides without any reliance on the genetic code. Notes
# the general usefulness of these NRPs. Therefore sequence these peptides using mass spectrometry. Approximate mass with
# the integer mass (notes in reality we would be using non-integer masses). Goal is to generate the experimental
# spectrum of masses after mass spec breaks down the protein. Also make simplifying assumption that the cyclic peptide
# is broken at every 2 bonds, so the resulting experimental spectrum contains all of the linear fragments of the peptide
# (called subpeptides).  The theoretical spectrum of a peptide would be the collection of masses for all subpeptides, as
# well as 0, and the full peptide. First, generate linear spectrum, which can be found based on the assumption that the
# mass of any subpeptide can be found by subtracting the masses of two prefixes. generate an array of all prefix masses
# in increasing order, then the mass of subpeptide beginning at position i + 1 and ending at j can be found by
# prefix_mass[j] - prefix_mass[i]
# NQEL subpeptides -> [-, N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, LNQ, NQEL]
# NQEL prefixes -> [-, N, NQ, NQE, NQEL]
# ----------------------------------------------------------------------------------------------------------------------


def get_aa_mass():
    amino_acids = "GASPVTCILNDKQEMHFRYW"
    mass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]
    mass_dict = dict(zip(amino_acids, mass))
    return mass_dict


def linear_spectrum(peptide, aa_mass):
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


# The theoretical spectrum of a circular peptide can then be found by finding the linear spectrum, as well as
# the subpeptides wrapping around the end of the linearized version of peptide. Each such subpeptide
# has mass equal to the difference between mass(peptide) and the subpeptide mass identified by linear spectrum. Eg,
# for NQEL: mass(LN) = mass(NQEL) - mass(QE)
# ----------------------------------------------------------------------------------------------------------------------


def cyclic_spectrum(peptide, aa_mass):
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


# However we are trying to solve the inverse problem - given a spectrum, reconstruct the protein. Start by assuming we
# got lucky and the mass spec generated an ideal spectrum - one that coincides with the peptide's theoretical spectrum.
# Also note that from now on deal with the masses themselves - redundancies because I/L and K/Q have the same integer
# mass. Also assume that the largest mass in a spectrum corresponds to the parent mass (full peptide of interest).
# Can employ a brute force algo for cyclopeptide sequencing that tests every peptide, but this is highly impractical:
# ----------------------------------------------------------------------------------------------------------------------


def pep_to_mass_str(peptide_str, aa_mass):
    mass_string = [str(aa_mass[aa]) for aa in peptide_str]
    return "-".join(mass_string)


def get_peptide_mass(peptide, aa_mass):
    mass = 0
    for aa in peptide:
        mass += aa_mass[aa]
    return mass


def bf_cyclo_seq(spectrum, candidate_peptides, aa_mass):
    peptides = []
    for peptide in candidate_peptides:
        if get_peptide_mass(peptide) == sum(spectrum) and spectrum == cyclic_spectrum(peptide, aa_mass):
            peptides.append(peptide)
    return peptides


# So turn to branch and bound algorthims: This set of algorithms grows/branches list of candidates then uses a bounding
# step to remove hopeless candidates. Grow linear peptides that remain consistent with the theoretical spectrum,
# and if so circularize to see if the cyclical spectrum matches. A linear peptide is consistent with spectrum if every
# mass in its theoretical spectrum is contained in spectrum. If a mass appears more than once in the theoretical
# spectrum of a linear peptide, it must appear at least that many times in spectrum.
# ----------------------------------------------------------------------------------------------------------------------


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


# The brute force implementation is exponential, and while the branch and bound version is practically much faster,
# it can generate many incorrect intermediate kmers and no one has been able to prove that it is polynomial, so from a
# theoretical perspective it is just as slow. Further, cyclopeptide_sequencing assumes that the experimental spectrum
# perfectly matches the theoretical spectrum, which is not realistic - mass spec produces false and missing masses. ANY
# false/missing masses will cause this algo to throw away the correct peptide because the theoretical spectrum differs
# from the experimental. To account for this noisiness, employ a scoring function to count matching masses.
# ----------------------------------------------------------------------------------------------------------------------


def cyclopeptide_scoring(peptide, spectrum, aa_mass):
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = cyclic_spectrum(peptide, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def cyclopeptide_scoring_int(peptide_int, spectrum, aa_mass):
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = cyclic_spectrum_int(peptide_int, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


# The goal is to then adapt cyclopeptide sequencing to find the peptide with the highest score - must revise the bound
# step to include more linear candidates (account for potential misses/false masses), while eliminating those with
# insufficiently high score. Note that must calculate the score of the linear spectrum, not the cyclic. Then employ
# a trim function to retain the top N scoring peptides (while allowing ties)
# ----------------------------------------------------------------------------------------------------------------------


def linear_peptide_score(peptide, spectrum, aa_mass):
    score = 0
    spectrum_cp = spectrum.copy()
    peptide_spectrum = linear_spectrum(peptide, aa_mass)
    for mass in peptide_spectrum:
        if mass in spectrum_cp:
            score += 1
            spectrum_cp.remove(mass)
    return score


def trim(leaderboard, spectrum, n, aa_mass):
    scores = [linear_peptide_score(peptide, spectrum, aa_mass) for peptide in leaderboard]
    score_df = pd.DataFrame({"Score": scores}, index=leaderboard)
    score_df = score_df.sort_values("Score", ascending=False)
    n_score = score_df["Score"].iloc[n]  # the score of the nth element to keep
    leading_peptides = score_df.index[score_df["Score"] >= n_score]
    return leading_peptides.tolist()


# Text uses string text for peptides when talking about scoring, then switches to int for the actual seq. Annoying.
# ----------------------------------------------------------------------------------------------------------------------


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
    n = min(n, len(leaderboard))
    leaderboard = list(leaderboard)
    scores = [linear_peptide_score_int(peptide, spectrum, aa_mass) for peptide in leaderboard]
    sorted_ix = np.argsort(scores)[::-1]
    sorted_scores = [scores[ix] for ix in sorted_ix]
    sorted_peptides = [leaderboard[ix] for ix in sorted_ix]
    score_threshold = sorted_scores[n-1]
    for i in range(n, len(leaderboard)):
        if sorted_scores[i] < score_threshold:
            return sorted_peptides[:i]
    return sorted_peptides


def leaderboard_cyclopeptide_seq(spectrum, n, aa_mass):
    leaderboard = set(spectrum).intersection(set(aa_mass.values()))  # init with aa in spectrum
    lead_peptide, hi_score = "", 0
    parent_mass = max(spectrum)
    while leaderboard:
        leaderboard = grow_peptide(leaderboard, aa_mass)  # branch step
        for peptide in leaderboard.copy():
            if get_integer_mass(peptide) == parent_mass:
                peptide_score = cyclopeptide_scoring_int(peptide, spectrum, aa_mass)
                if peptide_score > hi_score:
                    lead_peptide = peptide
                    hi_score = peptide_score
            elif get_integer_mass(peptide) > parent_mass:
                leaderboard.remove(peptide)
        if len(leaderboard) > 0:   # bound step
            leaderboard = trim_int(leaderboard, spectrum, n, aa_mass)
    return lead_peptide


# Notes that the leaderboard sequencing is a heuristic that will degrade in performance when the number of errors 
# (missing or false masses) increases. A further wrinkle is that in reality, there are more than 20 of the proteinogenic
# amino acids (Selenocystein and Pyrrolysine bring it to 22), and NRPs can incorporate non-proteinogenic acids which 
# further increases chance that leaderboard incorporates an incorrect peptide with a similar weight. To get around this,
# only want to consider the relevant amino acids of a spectrum. Find the convolution of a spectrum: the positive 
# differences in masses of all subpeptides. If experimental contains NQE and NQ, you get mass of E even if it wasn't
# in the experimental spectrum. Sort the counts of occurrences of these masses to get the likely amino acids. Only
# keep aas between 57 and 200 (heuristic to keep likely aa)
# ----------------------------------------------------------------------------------------------------------------------


def spectral_convolution(spectrum):
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
    conv = np.array(spectral_convolution(spectrum))
    conv = conv[(conv > 57) & (conv < 200)]
    peptides, counts = np.unique(conv, return_counts=True)
    sort_ix = np.argsort(counts)[::-1]
    sorted_peptides = [peptides[ix] for ix in sort_ix]
    sorted_counts = [counts[ix] for ix in sort_ix]
    for i in range(conv_topn, len(sorted_counts)):
        if sorted_counts[i] < sorted_counts[conv_topn-1]:
            return sorted_peptides[:i]
    return sorted_peptides


# hacky: before was working with aa_mass as a dict that had aa names as keys and int mass as values. new aas identified
# from topn conv may be outside of these 20 aas - init a generic dict with the int masses as values again
def aa_intmass_dict(intmass_list):
    aa_dict = {}
    for i in range(len(intmass_list)):
        aa_dict[i] = intmass_list[i]
    return aa_dict


def convolution_cyclopeptide_seq(conv_topn, leaderboard_topn, spectrum):
    leaderboard = keep_topn_convolution(conv_topn, spectrum)
    aa_mass = aa_intmass_dict(leaderboard)
    lead_peptide = leaderboard_cyclopeptide_seq(spectrum, leaderboard_topn, aa_mass)
    return lead_peptide



# Example inputs
# ----------------------------------------------------------------------------------------------------------------------


# genetic_code = get_genetic_code(string_type="DNA")
aa_mass = get_aa_mass()

# spectrum = [0, 97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
# spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
# print(cyclopeptide_sequencing(spectrum, aa_mass))

# peptide = "NQEL"
# spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
# print(cyclopeptide_scoring(peptide, spectrum, aa_mass))
# print(linear_peptide_score(peptide, spectrum, aa_mass))

# leaderboard = ["LAST", "ALST", "TLLT", "TQAS"]
# spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
# top_n = 2

# top_n = 10
# spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
# peptide_int = "113-147-71-129"
# leaderboard_int = [pep_to_mass_str(peptide, aa_mass) for peptide in leaderboard]
# print(trim_int(leaderboard_int, spectrum, top_n, aa_mass))
# print(trim(leaderboard, spectrum, top_n, aa_mass))
# print(linear_peptide_score_int(peptide_int, spectrum, aa_mass))
# print(leaderboard_cyclopeptide_seq(spectrum, top_n, aa_mass))

# spectrum = [0, 137, 186, 323]
# output = spectral_convolution(spectrum)
# print(output)


# print(keep_topn_convolution(20, [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]))
tt = convolution_cyclopeptide_seq(20, 60, [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493])
print(tt)


