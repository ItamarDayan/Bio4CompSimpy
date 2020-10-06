from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc
from Helper_funcs import *
import Models.Script_generator


class SSP(Models.Script_generator.Script_generator):

    # ------------------------------------------bcs model notes (SSP)-------------------------------------#
    #   Model explanation:
    #   At first, the setup process is started. It produces beacons at y levels that
    #   correspond to the partial sums of the numbers given, these are the split junctions.
    #   These beacons transmit the value 1.
    #
    #   After the setup is done, the agents are created, each with a unique serial number, and then
    #   begin to travel in the system. At each step an agent checks his y level:
    #   If there is no active beacon on this y channel, then the agent proceeds to travel in the last direction that he took.
    #   If there is an active beacon on this y channel, then it means the agent reached a split junction,
    #   and then there is a 50% chance that he will split down, and 50% chance that he will split diagonally
    #
    #   When an agent reaches a y level that is equal to the SSP numbers sum, it creates a beacon
    #   on a channel called "done" and transmits through it his x level.

    def __init__(self, list, num_of_agents, bcs_file_name):
        super().__init__(list, num_of_agents, bcs_file_name)

    def print_setup(self):
        setup = "{s0![0],fast} \n \t \t \t \t"

        for i in range(1, self.list_length):
            setup += ".{s" + str(i) + "![0],fast}"
        setup += '\n \t \t \t \t'
        return setup

    def create_bcs_code(self):
        # Optimized simulation code for SSP (reduced beacon transmits)

        bcs_code = (self.print_constants() + '''\n\n        
                Setup[''' + self.print_setup_params() + '''] = ''' + self.print_setup() + '''.{start![1] ,fast};


                P[x,y,sum,last,serial] = {start?[1], r}.
                (
                [y==sum]->{done![x],fast} ||

                [y!=sum]->{y?[0..sum], fast}.({splitDown,1}.P[x,y+1,sum,0,serial] + {splitDiag,1}.P[x+1,y+1,sum,1,serial]) ||

                [y!=sum]->{~y?[0..sum],fast}.( [last == 0] -> {continueStraight, r}.P[x,y+1,sum,0,serial] + [last == 1]-> {continueStraight, r}.P[x+1,y+1,sum,1,serial])
                );


                Setup[''' + self.print_partial_sums() + ''']''' + "\n" + self.print_processes_start()) + ';'

        # ---------------------------------------------------------------------------------------------------#
        self.bcs_code = bcs_code

    # ***************************************************************************************************#
