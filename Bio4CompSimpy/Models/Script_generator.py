from pathlib import Path
import subprocess
import csv
import numpy
import matplotlib.pyplot as plt
import math
import time
import abc

from Bio4CompSimpy.Helper_funcs import Math_funcs, File_handling
from Bio4CompSimpy.Helper_funcs.Math_funcs import filter_empty_lists


class Script_generator:

    def __init__(self, list, num_of_agents, bcs_file_name):
        self.list = list
        self.list_sum = sum(list)
        self.list_length = len(list)
        self.num_of_agents = int(num_of_agents)
        self.bcs_file_name = bcs_file_name
        self.bcs_code = ""
        self.DIR_path = Path().absolute()
        self.elapsed_time = 0
        self.results = [[]]
        self.tagged = True  # set default
        self.agents_info = {}
        self.num_of_threads = 1
        self.num_of_simulations = 1
        self.sim_results = []

    # **************************** Print functions ******************************************************#

    def print_partial_sums(self):  # TODO: put this func in a differnet module(maybe math module) and change its name
        # The function returns a string that describes the partial sums of self.list
        partial_sums_list = self.get_partial_sums()

        partial_sums = "0"

        for item in partial_sums_list:
            partial_sums += ", " + str(item)

        return partial_sums

    def print_setup_params(self):
        setup_params = "s0"
        for i in range(1, self.list_length):
            setup_params += ",s" + str(i)

        return setup_params

    @abc.abstractmethod
    def print_setup(self):
        '''This method prints the bcs setup. Implementation in inheriting classes'''
        return

    def print_sim_results(self, i):
        results = self.results[i]
        elapsed_time = self.elapsed_time

        for counter, value in enumerate(results):
            print("Number of agents on output slot " + str(counter) + ": " + str(value) + '\n')

        print("Elapsed computation time: " + str(elapsed_time) + " seconds\n")

    def print_all_sims_results(self):

        print("\n--------------------------------------RESULTS--------------------------------------------\n\n")

        for i in range(0, len(self.results)):
            print("__________________________________Simulation #%d results_______________________________\n" % (i + 1))
            self.print_sim_results(i)
            print("\n________________________________________________________________________________________\n")

    def print_start_P_untagged(self):

        string = "|| num_of_agents*P[0,0," + str(self.list_sum) + ",0,0]\n"  # all agents have serial = 0
        return string

    def print_start_P_tagged(self):

        serial = 0
        string = ""

        for i in range(0, self.num_of_agents):
            string += "||P[0,0," + str(self.list_sum) + ",0," + str(serial) + "]\n"
            serial += 1

        return string

    def print_processes_start(self):

        if self.tagged:
            return self.print_start_P_tagged()
        else:
            return self.print_start_P_untagged()

    def print_constants(self):
        string = '''
        fast = 1000;
        r = 1;
        num_of_agents = ''' + str(self.num_of_agents) + ";"
        return string

    def create_bcs_code(self):

        bcs_code = self.print_constants() + '''


                Setup[''' + self.print_setup_params() + '''] = \n \t \t \t \t''' + self.print_setup() + '''.{start![1] ,fast};


                P[x,y,sum,last,serial] = {start?[1], r}.
                (
                [y==sum]->{done![x],fast} ||

                [y!=sum]->{y?[x], fast}.({splitDown,fast}.P[x,y+1,sum,0,serial] + {splitDiag,fast}.P[x+1,y+1,sum,1,serial]) ||

                [y!=sum]->{~y?[0],fast}.( [last == 0] -> {continueStraight, fast}.P[x,y+1,sum,0,serial] + [last == 1]-> {continueStraight, fast}.P[x+1,y+1,sum,1,serial]) ||

                [y!=sum]->{y?[0],fast}.{~y?[x],fast}.{blockedSplitDown,fast}.P[x,y+1,sum,0,serial]
                );


                Setup[''' + self.print_partial_sums() + ''']''' + "\n" + self.print_processes_start() + ";"

        # ---------------------------------------------------------------------------------------------------#
        self.bcs_code = bcs_code

    # ***************************************************************************************************#

    # ***************************** Interpret results functions *****************************************#

    def interpret_results_tagged(self, sim_results):
        # This function gets results of *one* simulation (csv obj) and returns slots count and agents info dict
        agents_info = {}
        DIR_path = self.DIR_path

        results = numpy.zeros(self.list_sum + 1, dtype=numpy.uint64)

        for line in sim_results:

            if "splitDown" in line or "splitDiag" in line or "done" in line:  # If this line is in our interest

                x = line[line.index("x") + 1]  # get x
                x = int(x)

                y = line[line.index("y") + 1]  # get y
                y = int(y)

                serial = line[line.index("serial") + 1]  # get the proccess serial number

                if not serial in agents_info:  # If the proccess is not in our dictionary, add it
                    agents_info[serial] = ""

                if ("done" in line):
                    results[x] += 1

                if "splitDown" in line:
                    agents_info[serial] += "[x:" + str(x) + ", y:" + str(y) + ", Down]"

                if "splitDiag" in line:
                    agents_info[serial] += "[x:" + str(x) + ", y:" + str(y) + ", Diag]"

                # TODO: create a class "agent_info" so we can easily get info on agents instead of this dictionary
        return results, agents_info

    def interpret_results_untagged(self, sim_results):
        # ----------------------------Read and interpret the results file----------------------------------------------#
        list_sum = self.list_sum
        results = numpy.zeros(list_sum + 1, dtype=numpy.uint64)

        for line in sim_results:
            if "P" in line:  # if this line describes a proccess (agent)
                if "done" in line:  # if the agent is done
                    index = int(line[line.index("x") + 1])
                    results[index] += 1  # add 1 to the slot counting
        # -------------------------------------------------------------------------------------------------------------#
        return results

    def interpret_results(self):

        sims_res_file = open(str(self.DIR_path) + '/simulationOutput')
        sims_res_string = sims_res_file.read()
        sims_lst = self.split_simulations(sims_res_string)

        results = [[]]
        agent_info = [{}]

        for i, sim in enumerate(sims_lst):

            sim_csv = [[]]
            sim_lines = sim.split('\n')
            for k, line in enumerate(sim_lines):
                sim_csv[k] += line.split('\t')
                sim_csv.append([])

            if self.tagged:
                results[i], agent_info[i] = self.interpret_results_tagged(sim_csv)
                results.append([])
                agent_info.append({})

            else:  # not tagged
                results[i] = self.interpret_results_untagged(sim_csv)
                results.append([])

        self.results = filter_empty_lists(results)  # remove redundent empty lists
        self.agents_info = filter_empty_lists(agent_info)

    def split_simulations(self, simulationOutput):
        simulationOutput = str(simulationOutput)
        sim_lst = simulationOutput.split(">=======")
        sim_lst = sim_lst[1:]  # the file starts with a ">=====" so we get rid of the first empty string
        return sim_lst

    # ***************************************************************************************************#

    # ***************************** Set functions *******************************************************#

    def set_num_of_threads(self, num):
        self.num_of_threads = str(num)

    def set_num_of_simulations(self, num):
        self.num_of_simulations = str(num)

    def set_tagged(self, bool):
        self.tagged = bool

    def set_results(self, results):
        self.results = results

    # **************************************************************************************************#

    # **************************** Files and simulation handling ***************************************#

    def run_simulation(self):
        # ------------------------ Run bcs simulation -----------------------------------#
        cmd = (str(self.DIR_path) + '/../bin/bcs -t ' + str(self.num_of_threads) + ' -s '  # Create cmd
               + str(self.num_of_simulations) + " " + self.bcs_file_name)

        elapsed_time = File_handling.run_shell_command(cmd)

        self.elapsed_time = elapsed_time
        # -------------------------------------------------------------------------------#

    def generate_and_run(self):
        self.create_bc_file()
        self.run_simulation()

    def create_bc_file(self):
        File_handling.write_new_file(self.bcs_file_name, self.bcs_code)

    def save_plot(self, results, i):

        plt.figure()
        plt.bar(range(0, len(results)), results)
        plt.xlabel('Output slots')
        plt.ylabel('Number of agents')
        fig_name = 'resultsFig#' + str(i)
        plt.savefig(fig_name)

    def save_plots(self):
        for i, res in enumerate(self.results):
            self.save_plot(res, i)

    # **************************************************************************************************#

    # ***************************** Browsing functions *************************************************#

    def browse_result_slots(self, indx):
        # this function lets the user browse the number of agents in each slot
        results = self.results[indx]

        while True:
            try:
                answer = int(input("To get a specific result, enter the slot number; to go back enter -1 \n"))

            except:
                print("invalid input\n")
                continue

            if answer == -1:
                break

            try:
                print("Number of agents on output slot " + str(answer) + ": " + str(results[int(answer)]) + '\n')

            except:
                if answer >= len(results) or answer < -1:
                    print("Index out of results range\n")
                    continue
                else:
                    print("unknown error")

    def browse_results_tagged(self, indx):
        # this function lets the user browse the number of agents in each slot and the path of every agent
        agents_info = self.agents_info[indx]
        while True:
            answer = input(
                "To see a specific slot, enter 1;\nTo see an agent's route, enter 2 \nto go back enter -1 \n\n")
            answer = int(answer)
            if answer == -1:
                return

            if answer == 1:
                self.browse_result_slots(indx)

            if answer == 2:
                answer = input("Enter agent's serial number ")
                answer = int(answer)
                if answer > self.num_of_agents - 1 or answer < 0:
                    print("Invalid serial number")
                    continue

                print("Route taken by agent: ")
                print(agents_info.get(str(answer)) + '\n')

    def browse_sim_results(self, indx):
        if self.tagged:
            self.browse_results_tagged(indx)
        else:
            self.browse_result_slots(indx)

    def browse_main_menu(self):

        while True:
            try:
                sim_to_browse_index = int(input("Enter index of simulation to browse. to exit enter -1\n")) - 1

                if sim_to_browse_index == -2:
                    return
                else:
                    self.browse_sim_results(sim_to_browse_index)

            except:
                print("invalid input. try again\n")
                continue

    # **************************************************************************************************#

    def get_partial_sums(self):
        # The method returns a list of the partial sums of self.list
        lst = self.list
        return Math_funcs.partial_sums(lst)

