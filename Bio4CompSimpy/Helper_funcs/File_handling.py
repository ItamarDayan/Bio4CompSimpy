from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc


def run_shell_command(cmd):

    cmd = str(cmd)

    start = time.time()  # get time before computation

    subprocess.run(cmd, shell=True)  # Run BCS simulation on the generated code

    end = time.time()  # get time after computation

    elapsed_time = end - start  # estimate the computation time

    return elapsed_time


def write_new_file(file_name, data_to_write):

    file_name = str(file_name)
    data_to_write = str(data_to_write)

    # ---------------------------Create a new file for the generated bcs code ---------------------#
    file = open(file_name, "w")  # Open a new text file with write premissions
    file.write(data_to_write)  # Write the bcs code to the file
    file.close()  # close the file
    # ----------------------------------------------------------------------------------------------#


