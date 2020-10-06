
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