##################################################
# Spock 0: Initialisation, sort the files
##################################################

from MetaDataClass import BorgCollective
from sort_images import sort_the_images
import os

if __name__ == "__main__":

    meta = BorgCollective(spock=0)
    meta.read_tcf("Spock0_save_files.tcf")

    os.makedirs("Spock0_calib_output", exist_ok=True)
    
    meta.write("Spock0_calib_output/reduction_log_Spock0.h5")

    # Function that sorts the files into lists
    sort_the_images(meta)
