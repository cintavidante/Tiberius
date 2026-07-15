###############################################################
# Flats-BadPixelMask: Master Flat Computation
###############################################################
"""
Computation of the Master Flat.
"""
from __future__ import annotations

import logging
import sys, os
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
from astropy.io import fits

from core.product import BorgCollective, find_existing_run, next_run_id
from core.registry import update_registry, find_latest_product
from core.instruments import load_instrument_config
from core.html_report import write_registry_html
from core.html_report import _write_table_preview

logger = logging.getLogger("Spock-1")

# ---------------------------------------------------------------------------
# FLAT AND BAD MASK HELPERS
# ---------------------------------------------------------------------------

def bias_subtraction(master_bias, data):

    return data - master_bias

# -----------------------------------

def save_fits(data, filename, overwrite=False):

    hdu = fits.PrimaryHDU(data)
    hdu.writeto(filename, overwrite=overwrite)

    return

def _frame_stats(arr: np.ndarray) -> Dict[str, float]:
    """
    Compute the Mean, Shape, Max_count
    
    Parameters
    ----------
    arr: numpy.ndarray
        Array with the bias frame

    Returns
    -------
    Dict[str, float]
        Dictionary with the statistics of the frame
    """
    return {
        "mean": float(np.mean(arr)),
        "shape": np.shape(arr),
        "max_count": float(np.max(arr[:,20:-20])), # Need to clip extreme edges which can have very high counts
    }

# ---------------------------------------------------------------------------
# FLAT COMBINATION FUNCTIONS
# ---------------------------------------------------------------------------

def combine_flats_1window(
    flats_list: List[str],
    flat_list_name: str, 
    master_bias: pd.DataFrame, 
    instr_cfg, 
    idx: List[int], 
    sat_limit: int
) -> Tuple[np.ndarray, dict[str, Any]]:
    """
    Master Flat for one window data.

    Parameters
    ----------
    flats_list: List[str]
        List of flat field frames
    flat_list_name: str
        Name of the flat field list
    master_bias: pandas.DataFrame
        Master Bias 
    instr_cfg
        Instrument parameters
    idx: List[int]
        List of the data extension
    sat_limit: int
        Saturation Limit
    
    Returns
    -------
    numpy.ndarray

    List[Dict[str, Any]]
    """

    id_slice = slice(*instr_cfg.id_slice)
    per_frame_rows: List[Dict[str, Any]] = []
    stack: List[np.ndarray] = []

    for n_frame, fp in enumerate(flats_list, start=1):
        with fits.open(fp) as hdul:
            arr = hdul[idx[0]].data

        frame_id = str(fp)[id_slice]
        row = {"path": str(fp), "frame_id": frame_id, **_frame_stats(arr), "list_name": flat_list_name}
        if np.max(arr[:,20:-20]) <= sat_limit:
            stack.append(arr - master_bias)
            row["selected"] = True
        else:
            logger.info(f'--- ingoring the potentially saturated frame {frame_id}', extra={"data_name": "flats"})
            row["selected"] = False
        hdul.close()
        per_frame_rows.append(row)

    mean_combine = np.mean(np.stack(stack, axis=0), axis=0)

    master_flat = mean_combine/np.median(mean_combine)

    return master_flat, per_frame_rows

# -----------------------------------

def combine_flats_2windows(flats_list, master_bias, project, instrument, instr_cfg, idx_list, sat_limit, 
                           showfig=False, savefig=False, flat_output=os.getcwd()):

    id_slice = slice(*instr_cfg.id_slice)
    flat_data = [[],[]]

    plt.figure(figsize=(10,6))

    for n_frame, f in enumerate(flats_list):

        hdul = fits.open(f)
        frame_id = f[id_slice]

        for w, idx in enumerate(idx_list):

            data_frame = hdul[idx].data

            print(f'File {n_frame+1}/{len(flats_list)} ; {frame_id} ; window {w} ; mean={np.mean(data_frame)} ; shape={np.shape(data_frame)} ; max_count={np.max(data_frame[:,20:-20])}')
            
            if np.max(data_frame[:,20:-20]) <= sat_limit:
                flat_data[w].append(data_frame-master_bias[w])
            else:
                print('--- ingoring potentially saturated frame')
            
            plt.subplot(1,2, w+1)

            vmin,vmax = np.nanpercentile(data_frame,[10,90])
            plt.imshow(data_frame, vmin=vmin, vmax=vmax, cmap='hot')
            if w == 0:
                plt.xlabel("X pixel")
                plt.ylabel("Y pixel")
            else:
                plt.yticks(visible=False)

        plt.colorbar()
        plt.suptitle(f'Frame{n_frame+1:03d} ; {frame_id} ; {project} ; {instrument}')
        if savefig:
            plt.savefig(flat_output + f'/flat_frame_{n_frame+1:03d}.png')
        
        if showfig:
            plt.show(block=False)
            plt.pause(1)
        
        plt.clf()

        hdul.close()
    
    plt.close()

    flat_data = np.array(flat_data)

    mean_combine = [np.mean(flat_data[0],axis=0),np.mean(flat_data[1],axis=0)]

    return mean_combine


def create_master_flats(
    borg: BorgCollective,
    flat_files: list[str],
    master_bias: pd.DataFrame,
    extensions: dict[str, Any],
    instr_cfg,) -> None:
    """
    Create the master flat(s)
    
    Parameters
    ----------
    flat_files: list[str],
        List with the flat field files
    master_bias: pd.DataFrame
        DataFrame containing the Master Bias
    extensions: dict[str, Any]
        Extension dictionary containing extension and nwin
    instr_cfg
        Instrument config parameters
    
    Returns
    -------
    None
    """

    satlimit = borg.meta.parameters["spock1_flats"]["saturation_limit"]
    nwin = extensions['nwin']
    flat_list_name = borg.meta.parameters["spock1_flats"]["list"]
    
    # Create master flat for different windows
    if nwin == 1:
        master_flat, per_frame_rows = combine_flats_1window(flat_files, flat_list_name, master_bias, instr_cfg, extensions['idx_list'], satlimit)
        borg.ds = xr.Dataset(
            {"master_flat": xr.DataArray(master_flat.astype(np.float32), dims=("y", "x"))}
        )
    else:
        master_flat, per_frame_rows = combine_flats_2windows(flat_files, master_bias, instr_cfg, extensions['idx_list'], satlimit)
        borg.ds = xr.Dataset(
            {"master_flat": xr.DataArray(master_flat.astype(np.float32), dims=("window", "y", "x"))}
        )

    qa = pd.DataFrame(per_frame_rows)
    nfiles = qa["selected"].sum()

    summary = {
        "path": "__Master Flat__",
        "frame_id": "",
        "n_files": int(nfiles),
        **_frame_stats(master_flat),
        "list_name": borg.meta.parameters["spock1_flats"]["list"],
    }

    if instr_cfg.instrument_id == "WHT_ACAM":
        summary["n_windows"] = int(nwin)

    qa = pd.concat([qa, pd.DataFrame([summary])], ignore_index=True)
    borg.set_qa(qa)

    return

# ---------------------------------------------------------------------------
# Optional: Testing the script
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    meta = xr.Dataset()
    meta.attrs['instrument'] = 'EFOSC2'
    meta.attrs['projet_name'] = 'Planet'
    meta.attrs['inputdir_Spock1'] = 'Spock0_calib_output'
    meta.attrs['outputdir_Spock1'] = 'Spock1_pre_processing'
    
    meta.attrs['flat_list'] = '/Users/cividante/Tiberius/hatp65/Spock0_calib_output/SKY_FLAT_Gr11_27arcsec_list'
    meta.attrs['flat_savefits'] = True
    meta.attrs['flat_showfig'] = True
    meta.attrs['flat_savefig'] = True
    meta.attrs['flat_saturation_limit'] = 55000
    meta.attrs['flat_box_width'] = 9

    meta.attrs['bad_pixel_mask'] = True
    meta.attrs['bad_pixel_mask_savefits'] = True
    meta.attrs['bad_pixel_mask_showfig'] = True
    meta.attrs['bad_pixel_mask_savefig'] = True
    meta.attrs['bad_pixel_mask_medianfilter_cutoff'] = 5


    create_master_flat_pixel_mask(meta, None)
