from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc


from Bio4CompSimpy.Helper_funcs import Math_funcs
from Bio4CompSimpy.Models.Script_generator import Script_generator


class Exact_cover(Script_generator):

    # ------------------------------------------Create the bcs code (Exact Cover)-------------------------------------#
    #   Model explanation:
    #   At first, the setup process is started. It produces beacons at y levels that
    #   correspond to the partial sums of the sets encoding. For each set there is a calculation
    #   of what x level bits don't override the set encoding. A set transmits over his beacon
    #   an x value if and only if this x value bits do not override the set encoding's bits.
    #
    #   After the setup is done, the agents are created, each with a unique serial number, and then
    #   begin to travel the network. At each step an agent checks his y level:
    #   If there is no active beacon on this y channel that transmits x, then the agent proceeds to travel in the last direction that he took.
    #
    #   If there is an active beacon on this y channel that transmits x, then it means the agent reached a split junction,
    #   and then there is a 50% chance that he will split down, and 50% chance that he will split diagonally.
    #   If there is an active beacon on this y channel but it does not transmit x, then it means the agent reached
    #   a blocked junction, then it will proceed in the last direction that it took.
    #
    #   When an agent reaches a y level that is equal to the exact cover groups encoding sum, it creates a beacon
    #   on a channel called "done" and transmits through it his x level.

    def __init__(self, list, num_of_agents, bcs_file_name, group_size):
        super().__init__(list, num_of_agents, bcs_file_name)
        self.group_size = group_size

    def print_setup(self):

        max_junctions = int(math.pow(2, self.group_size)) - 1
        partial_sums = self.get_partial_sums()

        setup = "{s0![0],fast} \n \t \t \t \t"

        for i in range(0, self.list_length - 1):
            junctions_on_level = min(max_junctions, partial_sums[i])

            for x in range(0, junctions_on_level + 1):  # until junctions_on_level including
                if not Math_funcs.bits_override(self.list[i + 1], x):
                    setup += ".{s" + str(i + 1) + "![" + str(x) + "],fast}"

            setup += "\n \t \t \t \t"

        return setup

    @staticmethod
    def encode_group(n):
        group = 0
        size = int(n)
        for i in range(0, size):
            bit = int(input())
            while bit != 0 and bit != 1:
                print("invalid input. Enter 1 if item is in group or 0 else\n")
                bit = int(input())

            group += int(math.pow(2, i)) * bit

        return group

