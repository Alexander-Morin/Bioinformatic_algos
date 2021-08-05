# Implementing code to find the spectral convolution of an input peptide: all of the positive differences of masses in 
# the spectrum. These subpeptides are then ranked in order of occurrence to indicate the most likely amino acids of
# an experimental spectrum.

# Problem 4H in the BALA textbook/Rosalind

# Input is text file containing the integer masses of the input spectrum

# Output is the spectral convolution sorted in order of occurrence

# Usage: python3 4H_spectral_convolution.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
from collections import Counter


def spectral_convolution(spectrum):
    """
    spectrum: list of integers corresponding to the input spectrum
    returns a list of integers sorted by element occurrence
    """
    conv = []
    for i in range(0, len(spectrum)):
        for j in range(0, len(spectrum)):
            if i == j:
                next
            elif spectrum[j] - spectrum[i] > 0:
                conv.append(spectrum[j] - spectrum[i])
    conv = sorted(conv, key=Counter(conv).get, reverse=True)
    return conv


def main():
    """
    Read the input file, parse the arguments, and finds the spectral convolution
    """
    argv = list(sys.argv)
    input_spectrum = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_spectrum is None:
                input_spectrum = line.replace('\n', '')
                
    input_spectrum = input_spectrum.split(" ")
    input_spectrum = [int(x) for x in input_spectrum]
    output = spectral_convolution(input_spectrum)
    
    for i in output:
      print(i, end=" ")


if __name__ == "__main__":
    main()
