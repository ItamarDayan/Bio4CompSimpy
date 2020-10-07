import os
import pathlib
import sys

#------------------ Setup enviroment variable--------------------#
path = str(pathlib.Path().absolute().parent)
try:
    sys.path.index(path)
except ValueError:
    sys.path.append(path)
#----------------------------------------------------------------#


from Bio4CompSimpy.Helper_funcs import User_IO
#   ______________________________Main______________________________________________

if __name__ == "__main__":

    scrptG = User_IO.menu()

    scrptG.create_bcs_code()

    scrptG.generate_and_run()

    scrptG.interpret_results()

    scrptG.save_plots()

    scrptG.print_all_sims_results()

    scrptG.browse_main_menu()
