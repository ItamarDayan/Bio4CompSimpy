from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc



def partial_sums(lst):
    # The method returns a list of the partial sums of lst
    n = len(lst)

    sum = 0
    partial_sums = []

    for i in range(0,n-1):
        sum += lst[i]
        partial_sums.append(sum)

    return partial_sums


def calc_first_n_primes(n):
    primes = [2]
    n = int(n)
    n -= 1
    i = 3
    is_prime = True

    while n > 0:
        is_prime = True
        for k in range(2, math.ceil(i/2)+1):
            if i % k == 0:
                is_prime = False
                break

        if(is_prime):
            n -= 1
            primes.append(i)

        i += 1

    return primes


def bits_override(g1, g2):
    if g1|g2 == g1^g2:
        return False
    else:
        return True

def filter_empty_lists(the_list):
    return [elem for elem in the_list if len(elem)!=0]
