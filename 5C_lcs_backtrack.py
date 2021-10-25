# Implementing a dynamic programming approach to find the longest common subsequence between two strings. A longest
# common subsequence allows insertions/deletions of characters between the two strings. This is represented by an
# alignment graph, which can be further broken down into a scoring and backtrack matrix of dim
# len(str1)+1 x len(string2)+1. For each element of each string, if there is a match, it is a 1 in the score matrix, and
# a 0 otherwise. The backtrack matrix keeps track of these matches, where "Diagonal" represents a match, "Down"
# represents a deletion (string2 has a gap for the current character of string 1), and "Right" represents an insertion
# (a gap in string1 where string2 has a character). After filling these matrices, a recursive output function then
# starts at the sink, working backwards to find and record the "Diagonal" elements of the backtrack matrix which
# represent a match.

# Problem 5C in the BALA textbook/Rosalind

# Input is a text file where the first line is string1, and the second line is string2.
# AACCTTGG
# ACACTGTGA

# Output is the LCS between string1 and string2.
# AACTTG

# Usage: python3 5C_lcs_backtrack.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import numpy as np

sys.setrecursionlimit(1500)  # default setting was running into max recursion error


def lcs_backtrack(string1, string2):
    """
    string1: str that represents the rows of the backtrack/scoring matrix
    string2: str that represents the cols of the backtrack/scoring matrix
    returns a len(str1)+1 x len(str2)+1 backtrack matrix
    """
    scores = np.zeros([len(string1)+1, len(string2)+1])
    backtrack = np.empty([len(string1)+1, len(string2)+1], dtype="object")

    for i in range(1, len(string1)+1):
        for j in range(1, len(string2)+1):
            # fill score matrix, considering diagonal score + 1 if the strings match
            if string1[i-1] == string2[j-1]:
                scores[i, j] = max(scores[i-1, j-1] + 1, scores[i-1, j], scores[i, j-1])
            else:
                scores[i, j] = max(scores[i-1, j], scores[i, j-1])
            # fill backtrack matrix: down=deletion, right=insertion, diagonal=match
            if scores[i, j] == scores[i-1, j]:
                backtrack[i, j] = "Down"
            elif scores[i, j] == scores[i, j-1]:
                backtrack[i, j] = "Right"
            elif (scores[i, j] == scores[i-1, j-1] + 1) & (string1[i-1] == string2[j-1]):
                backtrack[i, j] = "Diagonal"

    return backtrack


def output_lcs(backtrack, string1, i, j):
    """
    backtrack: the backtrack matrix
    string1: the source string in which to find the LCS
    i: int representing the length of string1
    j: int representing the length of string2
    returns the LCS between string1 and string2 as determined by the backtrack matrix

    uses a recursive call to start at the "sink" of the backtrack matrix, working its way to the source and adding
    matching characters to output if there is a match
    """
    output = []

    def recursive_output(backtrack, string1, i, j, output):
        if i == 0 or j == 0:
            return None
        if backtrack[i, j] == "Down":
            recursive_output(backtrack, string1, i-1, j, output)
        elif backtrack[i, j] == "Right":
            recursive_output(backtrack, string1, i, j-1, output)
        else:
            recursive_output(backtrack, string1, i-1, j-1, output)
            output.append(string1[i-1])

    recursive_output(backtrack, string1, i, j, output)
    return "".join(output)


def main():
    """
    Read the input file, parse the arguments, and finds the LCS between string1 and string2
    """
    argv = list(sys.argv)
    input_string1 = None
    input_string2 = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_string1 is None:
                input_string1 = line.replace('\n', '')
            elif input_string2 is None:
                input_string2 = line.replace('\n', '')

    backtrack = lcs_backtrack(input_string1, input_string2)
    print(output_lcs(backtrack, input_string1, len(input_string1), len(input_string2)))


if __name__ == "__main__":
    main()
