##################################################
# Spock 0: Initialisation, sort the files
##################################################

from global_utils.MetaDataClass import BorgCollective
from Spock0_save_files.sort_images import sort_the_images
import os

def run_Spock0_save_files():
    meta = BorgCollective(spock=0)
    meta.read_tcf("Spock0_save_files.tcf")

    os.makedirs("Spock0_calib_output", exist_ok=True)
    
    meta.write("Spock0_calib_output/reduction_log_Spock0.h5")

    # Function that sorts the files into lists
    sort_the_images(meta)


if __name__ == "__main__":

    run_Spock0_save_files()