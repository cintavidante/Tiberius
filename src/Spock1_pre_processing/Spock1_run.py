####################################################################
# Spock 1: Pre-processing the data
####################################################################

import xarray as xr
import os
from global_utils.utils import get_bool_attr
from global_utils.MetaDataClass import BorgCollective
from global_utils.InstrumentDataClass import load_instrument_config
from Spock1_pre_processing.master_bias import create_master_bias
from Spock1_pre_processing.master_flat import create_master_flat_pixel_mask
#from ... import locate_cosmics
#from ... import bad_pixel_mask

def run_Spock1_pre_processing():
    # Read in the metadata from previous spock and combine it with current spock
    meta = BorgCollective(spock=1)
    meta.read_tcf("Spock1_pre_processing.tcf")

    previous_meta = xr.open_dataset("Spock0_calib_output/reduction_log_Spock0.h5", engine='h5netcdf')
    meta.combine(previous_meta)

    # Read in the Instrument static parameters
    instr_cfg = load_instrument_config(meta.ds.attrs["instrument"].strip("'"))

    # Build the 
    os.makedirs("Spock1_pre_processing", exist_ok = True)

    if get_bool_attr(meta.ds.attrs["skip_bias"]) == False:
        print('Starting Master-Bias.')
        create_master_bias(meta, instr_cfg)

    if get_bool_attr(meta.ds.attrs["skip_flats_bad_pixel_mask"]) == False:
        print('\n')
        print('Starting Flats and Bad-Pixel-Mask.')
        create_master_flat_pixel_mask(meta, instr_cfg)

    #if get_bool_attr(meta.ds.attrs["skip_cosmics"]) == False:
    #


    meta.write("Spock1_pre_processing/reduction_log_Spock1.h5")


if __name__ == "__main__":

    run_Spock1_pre_processing()
    