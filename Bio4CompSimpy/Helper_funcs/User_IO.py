from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc
from Bio4CompSimpy.Helper_funcs import Math_funcs
from Bio4CompSimpy.Models.Exact_cover import Exact_cover
from Bio4CompSimpy.Models.SSP import SSP



def manual_SSP_input():
    # number of elements as input
    lst_length = int(input("Enter the length of the set and a set of numbers for the SSP: "))
    lst = []

    # iterating till the range
    for i in range(0, lst_length):
        ele = int(input())
        lst.append(ele)  # adding the element

    lst.sort()

    num_of_agents = input("Enter amount of agents")

    return lst, num_of_agents


def primes_input():
    lst_length = int(input("Enter the amount of first prime numbers to use as a set for the SSP"))
    lst = Math_funcs.calc_first_n_primes(lst_length)
    num_of_agents = input("Enter amount of agents")

    return lst, num_of_agents


def Exact_cover_input():
    lst_length = int(input("Enter number of groups\n"))
    size = int(input("Enter total number of items\n"))
    lst = []

    for i in range(0, lst_length):
        print("For group number " + str(i) + ", Enter 1 if item is in group and 0 o.w\n")
        lst.append(Exact_cover.encode_group(size))

    lst.sort()

    num_of_agents = input("Enter number of agents")

    return lst, num_of_agents, size


def menu():
    # creating an empty list
    lst = []
    mode = int(
        input("Choose simulation mode;\n 1: Manual SSP \n 2: Run SSP on first n primes \n 3: Exact Cover Problem"))
    #   TODO: Add option to perform simulation and results read on a manual written simulation file
    if mode == 1:
        lst, num_of_agents = manual_SSP_input()
        scrptG = SSP(lst, num_of_agents, "manual_SSP_generated_bcs_code.bc")  # Create a SSP object

    if mode == 2:
        lst, num_of_agents = primes_input()
        scrptG = SSP(lst, num_of_agents, "primes_SSP_generated_bcs_code.bc")  # Create a SSP object

    if mode == 3:
        lst, num_of_agents, size = Exact_cover_input()
        scrptG = Exact_cover(lst, num_of_agents, "Exact_cover_generated_bcs_code.bc", size)  # Create Exact_Cover obj

    tagged = input("tag agents? y: Yes, n: No")
    if tagged == "y":
        scrptG.set_tagged(True)
    elif tagged == "n":
        scrptG.set_tagged(False)

    num_of_threads = input("Enter number of threads to use")
    scrptG.set_num_of_threads(num_of_threads)

    num_of_simulations = input("Enter number of simulations to run")
    scrptG.set_num_of_simulations(num_of_simulations)

    return scrptG
