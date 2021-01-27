# Implementing code to spell a string by a genomic path: output a string (text) of length k+n-1 such that that ith kmer
# in text is equal to pattern(i) (for 1 <= i <= n)
# Problem 3B in the BALA textbook/Rosalind

# Input is a text file with each line as a pattern in text
# ACCGA
# CCGAA
# CGAAG
# GAAGC
# AAGCT

# Output is the string, text
# ACCGAAGCT

# Usage: python3 3B_spell_string_by_path.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput


def spell_string_by_path(text):
    """
    text: list of DNA strings assumed to be overlapping by len(string)-1
    returns the complete glued string
    """
    string = [text[0]]
    for pattern in text[1:]:
        string.append(pattern[-1])
    return "".join(string)


def main():
    """
    Read the input file and parse the arguments.
    """
    argv = list(sys.argv)
    input_text = []

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            input_text.append(line.replace('\n', ''))

    output = spell_string_by_path(input_text)

    print(output)


if __name__ == "__main__":
    main()
