# Implementing dynamic programming solution to the greedy change problem - finding the minimum amount of coins to
# equal an amount of money, given a list of coins. This problem was used to motivate thinking about dynamic programming
# as a method for comparing DNA sequences. It was contrasted with the greedy implementation of simply choosing the
# largest possible coin at each step (which isn't guaranteed to find the right answer), and the top-down/recursive
# implementation, which results in-recalculating the same values over and over again.

# Problem 5A in the BALA textbook/Rosalind

# Input is a text file where the first line is the amount of money to make change, and the second line are the coins
# 40
# 1,5,10,20,25,50

# Output is the minumum number of coins
# 2

# Usage: python3 5A_dp_change.py input.txt > output.txt
# ----------------------------------------------------------------------------------------------------------------------

import sys
import fileinput
import numpy as np


def dp_change(money, coins):
    """
    money is an integer of the amount to make change
    coins is a list of integers of the coin denominations
    """
    min_num_coins = [0] * (money + 1)
    for m in range(1, money + 1):
        min_num_coins[m] = np.inf
        for i in range(len(coins)):
            if m >= coins[i]:
                if min_num_coins[m - coins[i]] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coins[i]] + 1
    return min_num_coins[money]


def main():
    """
    Read the input file, parse the arguments, and find the lead peptide
    """
    argv = list(sys.argv)
    input_money = None
    input_coins = None

    for line in fileinput.input(argv[1]):
        if len(line) > 0:
            if input_money is None:
                input_money = int(line.replace('\n', ''))
            elif input_coins is None:
                input_coins = line.replace('\n', '').split(",")

    input_coins = [int(x) for x in input_coins]
    output = dp_change(input_money, input_coins)
    print(output)


if __name__ == "__main__":
    main()
