# Implementing code XX

# Problem 4M in the BALA textbook/Rosalind

# Input is XX

# Output is XX

# Usage: python3 4M_turnpike_problem.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import numpy as np


def point_differences(int_list):
    diff_list = []
    for i in range(0, len(int_list)):
        for j in range(0, len(int_list)):
            diff_list.append(int_list[j] - int_list[i])
    return sorted(diff_list)


def test_points_are_consistent(ref_list, test_list):
    ref_cp = ref_list.copy()
    test_cp = test_list.copy()
    for i in test_list:
        if i in ref_cp:
            test_cp.remove(i)
            ref_cp.remove(i)
    return len(test_cp) == 0


def add_points(candidate_set, test_points, test):
    candidate_cp = candidate_set.copy()
    test_points.append(test)
    candidate_cp.remove(test)
    return candidate_cp, test_points



def turnpike_problem(diff_list):
    diff = np.array(diff_list)
    diff = list(diff[diff >= 0])
    n_term = diff.count(0)  # 0s from integers subtracted from themselves
    points = [0, diff[-1]]
    test_points = points.copy()
    candidates = set(abs(np.array(diff_list))).difference(points)
    while candidates:
        test_candidates = candidates.copy()
        test = max(test_candidates)
        test_candidates, new_points = add_points(test_candidates, test_points, test)
        print("after expansion", test_candidates, new_points)
        if test_points_are_consistent(diff_list, new_points):
            test = max(test_candidates)
            test_candidates, new_points = add_points(test_candidates, test_points, test)
            print("if statement", test_candidates, new_points)
        else:
            candidates.remove(test)
            print("else statement", test_candidates, new_points)
    return test_points


def main():
    """
    Read the input file, parse the arguments, and find the lead peptide
    """
    argv = list(sys.argv)
    input_diff = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_diff is None:
                input_diff = line.replace('\n', '').split(" ")

    input_diff = [int(x) for x in input_diff]
    # test_diff = point_differences([0, 6, 8, 10])
    # print("ref diff", input_diff)
    # print("test_diff", test_diff)
    # print(test_points_are_consistent(input_diff, test_diff))
    print(input_diff)
    output = turnpike_problem(input_diff)
    print(output)


if __name__ == "__main__":
    main()
