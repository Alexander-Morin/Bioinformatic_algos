# Implementing dynamic programming solution to the Manhattan tourist problem: "find the longest path in a rectangular
# city." Can only move right and down, starting from the top left of the table (0,0 - the seed) and moving to the bottom
# right (i, j - the sink). Each path has a weight so choose the path with the highest max weights. DP approach finds
# all sub-solutions.

# Problem 5B in the BALA textbook/Rosalind

# Input is a text file where the first line is the dimensions of the grid, and following lines corresponds to the
# matrices of vertical and then horizontal edge weights, separated by a dash.
# 4 4
# 1 0 2 4 3
# 4 6 5 2 1
# 4 4 5 2 1
# 5 6 8 5 3
# -
# 3 2 4 0
# 3 2 4 2
# 0 7 3 3
# 3 3 0 2
# 1 3 2 2

# Output is the longest path
# 34

# Usage: python3 5B_manhattan_tourist.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import numpy as np


def manhattan_tourist(i, j, v_edge_weight, h_edge_weight):
    """
    i and j are integers corresponding to the row and col (respectively) dimensions of the DP table
    v_edge_weight and h_edge_weight are numpy array objects corresponding to the vertical and horizontal weight matrices
    returns the score at the i,jth position, which is the DP solution for the longest path/max weight
    """
    scores = np.zeros([i+1, j+1])
    for x in range(1, i+1):
        scores[x][0] = scores[x-1][0] + v_edge_weight[x-1][0]
    for y in range(1, j+1):
        scores[0][y] = scores[0][y-1] + h_edge_weight[0][y-1]
    for x in range(1, i+1):
        for y in range(1, j+1):
            scores[x][y] = max(scores[x-1][y] + v_edge_weight[x-1][y], scores[x][y-1] + h_edge_weight[x][y-1])
    return int(scores[i][j])


def main():
    """
    Read the input file, parse the arguments, and find the longest path
    """
    argv = list(sys.argv)
    input_i = None
    input_j = None
    input_vweights_list = []
    input_vweights_mat = None
    input_hweights_list = []
    input_hweights_mat = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_i is None:
                input_dim = line.replace('\n', '').split(" ")
                input_i = int(input_dim[0])
                input_j = int(input_dim[1])
            elif input_vweights_mat is None:
                input = line.replace('\n', '').split(" ")
                if input[0] == "-":
                    input_vweights_mat = np.array(input_vweights_list, dtype=int)
                else:
                    input_vweights_list.append(input)
            elif input_hweights_mat is None:
                input = line.replace('\n', '').split(" ")
                input_hweights_list.append(input)

    input_hweights_mat = np.array(input_hweights_list, dtype=int)
    output = manhattan_tourist(input_i, input_j, input_vweights_mat, input_hweights_mat)
    print(output)


if __name__ == "__main__":
    main()
