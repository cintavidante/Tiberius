######################################################
# Flats-BadPixelMask: Master Flat Smoothing Models
######################################################
"""
This script contains the running functions for the smoothing models
"""
from __future__ import annotations

import logging
import sys, os
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
from astropy.io import fits
from scipy.ndimage import median_filter, gaussian_filter
import matplotlib.pyplot as plt

from core.product import BorgCollective, find_existing_run, next_run_id
from core.registry import update_registry
from core.html_report import write_registry_html
from core.html_report import _write_table_preview

# Initialise the Logger. Will not be saved for Spock1
logger=logging.getLogger("Spock-1")

# =======================================
# Helper functions
# =======================================
def test_smoothing_widths(flat_data, nwindows):

    """Run a series of median filters with differing box widths to find optimal width."""

    if nwindows > 1:
        nrows,ncols1 = np.shape(flat_data[0])
        _,ncols2 = np.shape(flat_data[1])
        std1 = []
        std2 = []
    else:
        nrows,ncols = np.shape(flat_data)
        std = []

    for i in range(1,nrows//10,2):

        box_width = i

        if nwindows == 1:
            MF = median_filter(flat_data[:,ncols//2],box_width)
            residuals = flat_data[:,ncols//2]/MF
            # rms.append(np.sqrt(np.mean(residuals**2)))
            std.append(np.std(residuals))

        else:
            MF1 = median_filter(flat_data[0][:,ncols1//2],box_width) # just use single column for plotting
            MF2 = median_filter(flat_data[1][:,ncols2//2],box_width)

            residuals1 = flat_data[0][:,ncols1//2]/MF1
            residuals2 = flat_data[1][:,ncols2//2]/MF2

            std1.append(np.std(residuals1))
            std2.append(np.std(residuals2))

    plt.figure(figsize=(10,8))
    plt.xlabel('Bin width')
    plt.ylabel('Standard deviation')
    plt.title("Plotting bin width vs. standard deviation - Please enter the DESIRED BIN WIDTH!")
    if nwindows == 1:
        plt.plot(np.arange(1,nrows//10,2),std,'ko')
    else:
        plt.plot(np.arange(1,nrows//10,2),std1,'bo',label='Window 1')
        plt.plot(np.arange(1,nrows//10,2),std2,'ro',label='Window 2')
        plt.legend(loc='upper left')
    
    plt.show(block=False)
    logger.info("Please enter the desired bin width", extra={"data_name": "flats"})
    best_box_width = input("Input of best box_width (integer): ").strip().lower()

    logger.info(f"Input of the best box_width: {best_box_width}", extra={"data_name": "flats"})
    plt.close("all")

    return best_box_width

# -----------------------------------

def median_smooth(self, flat_data, nwindows, box_width):
    """Smooth the flat using a running median, and evaluated at each column individually"""

    if nwindows == 1:
        _, ncols = np.shape(flat_data)
        # Now find running median for each column individually
        MF_reshaped =  np.array([median_filter(flat_data[:,i],box_width) for i in range(ncols)]).transpose()

        MF_reshaped[MF_reshaped == 0] = 1 # have to replace 0s which occur near the edges, the bias over subtracts to give negative values which mess up this line

        divided_sky_flat = flat_data/MF_reshaped
        normalised_sky_flat = divided_sky_flat/divided_sky_flat.mean()

        self.ds= xr.Dataset(
            {"master_flat_median_smooth": xr.DataArray(normalised_sky_flat.astype(np.float32), dims=("y", "x"))}
        )

    else:
        _, ncols1 = np.shape(flat_data[0])
        _, ncols2 = np.shape(flat_data[1])
        MF1_reshaped = np.array([median_filter(flat_data[0][:,i],box_width) for i in range(ncols1)]).transpose()
        divided_sky_flat1 = flat_data[0]/MF1_reshaped
        normalised_sky_flat1 = divided_sky_flat1/divided_sky_flat1.mean()

        MF2_reshaped = np.array([median_filter(flat_data[1][:,i],box_width) for i in range(ncols2)]).transpose()
        divided_sky_flat2 = flat_data[1]/MF2_reshaped
        normalised_sky_flat2 = divided_sky_flat2/divided_sky_flat2.mean()

        normalised_sky_flat = np.array([normalised_sky_flat1,normalised_sky_flat2])

        self.ds= xr.Dataset(
            {"master_flat_median_smooth": xr.DataArray(normalised_sky_flat.astype(np.float32), dims=("window", "y", "x"))}
        )

        
    return normalised_sky_flat

# -----------------------------------

def gaussian_smooth(self, flat_data, nwindows, sigma):

    """Smooth the flat using a Gaussian filter"""

    if nwindows == 1:
        GF = gaussian_filter(flat_data,sigma)
        divided_sky_flat = flat_data/GF
        normalised_sky_flat = divided_sky_flat/divided_sky_flat.mean()

        self.ds= xr.Dataset(
            {"master_flat_gaussian_smooth": xr.DataArray(normalised_sky_flat.astype(np.float32), dims=("y", "x"))}
        )

    else:
        GF1 = gaussian_filter(flat_data[0],sigma)
        divided_sky_flat1 = flat_data[0]/GF1
        normalised_sky_flat1 = divided_sky_flat1/divided_sky_flat1.mean()

        GF2 = gaussian_filter(flat_data[1],sigma)
        divided_sky_flat2 = flat_data[1]/GF2
        normalised_sky_flat2 = divided_sky_flat2/divided_sky_flat2.mean()

        normalised_sky_flat = np.array([normalised_sky_flat1,normalised_sky_flat2])

        self.ds= xr.Dataset(
            {"master_flat_gaussian_smooth": xr.DataArray(normalised_sky_flat.astype(np.float32), dims=("window", "y", "x"))}
        )
    
    return normalised_sky_flat


# ==============================================================
# Main Runner functions
# ==============================================================

def run_median_smooth(
    flat_run_dir,
    nwin,
    sp1_params,
    context,
    storage_format,
    overwrite,
):    
    # Get Master Flat
    df = BorgCollective.open_product(flat_run_dir)
    master_flat = df["master_flat"]

    # Test box width
    box_width = int(test_smoothing_widths(master_flat, nwin))

    borg = BorgCollective(
        spock_name="Median Smooth",
        spock=1,
        storage_format=storage_format,
    )

    borg.set_parent("spock0", str(flat_run_dir))

    borg.add_parameters({"spock1_flat": {
        "list": sp1_params.get("flat_list"),
        "model": "median",
        "box_width": box_width,
    }})

    borg.set_config_hash_from_parameters()
    borg.set_input_hash_from_parents([flat_run_dir])

    # Output location
    registry_dir = ( flat_run_dir / "models" / "median_smooth" )

    # Find the run with same hashes if existing
    existing = find_existing_run(registry_dir, borg.meta.product_id.config_hash, borg.meta.product_id.input_hash)
    if existing and not overwrite:
        logger.info(f"Skipped: Median_Smooth already exists -> {existing}", extra={"data_name", "flats"})
        return registry_dir / existing
    elif existing and overwrite:
        run_id = existing
    else: 
        run_id = next_run_id(registry_dir) 
    borg.meta.parameters["run_id"] = run_id

    # Running directory
    run_dir = registry_dir / run_id

    master_flat_median_smooth = median_smooth(borg, master_flat, nwin, box_width)

    borg.publish(run_dir)

    update_registry(registry_dir)

    hdu = fits.PrimaryHDU(master_flat_median_smooth)
    hdu.writeto(run_dir / 'master_flat_median_smooth.fits', overwrite=True)
    
    # Update the registry html file
    subtitle=f"Instrument: {context.get('instrument_id')} | Planet: {context.get('project_name')} | Date: {context.get('project_date')}"
    write_registry_html(registry_dir, spock=1, registry_parquet="_registry.parquet", out_html="_registry.html", title="Tiberius-Spock1-Flats-MedianSmooth", subtitle=subtitle)
   
    logger.info(f"New Median Smooth: {run_id}", extra={"data_name": "flats"})
    
    return run_dir


def run_gaussian_smooth(
    flat_run_dir,
    nwin,
    sp1_params,
    context,
    storage_format,
    overwrite,
):
    # Get Master Flat
    df = BorgCollective.open_product(flat_run_dir)
    master_flat = df["master_flat"]

    # Test box width
    sigma = sp1_params.get("flat_gaussian_sigma")

    borg = BorgCollective(
        spock_name="Median Smooth",
        spock=1,
        storage_format=storage_format,
    )

    borg.set_parent("spock0", str(flat_run_dir))

    borg.add_parameters({"spock1_flat": {
        "list": sp1_params.get("flat_list"),
        "model": "gaussian",
        "sigma": sigma,
    }})

    borg.set_config_hash_from_parameters()
    borg.set_input_hash_from_parents([flat_run_dir])

    # Output location
    registry_dir = ( flat_run_dir / "models" / "gaussian_smooth" )

    # Find the run with same hashes if existing
    existing = find_existing_run(registry_dir, borg.meta.product_id.config_hash, borg.meta.product_id.input_hash)
    if existing and not overwrite:
        logger.info(f"Skipped: Gaussian_Smooth already exists -> {existing}", extra={"data_name", "flats"})
        return registry_dir / existing
    elif existing and overwrite:
        run_id = existing
    else: 
        run_id = next_run_id(registry_dir) 
    borg.meta.parameters["run_id"] = run_id

    # Running directory
    run_dir = registry_dir / run_id

    master_flat_gaussian_smooth = gaussian_smooth(borg, master_flat, nwin, sigma)

    borg.publish(run_dir)

    update_registry(registry_dir)

    hdu = fits.PrimaryHDU(master_flat_gaussian_smooth)
    hdu.writeto(run_dir / 'master_flat_gaussian_smooth.fits', overwrite=True)
    
    # Update the registry html file
    subtitle=f"Instrument: {context.get('instrument_id')} | Planet: {context.get('project_name')} | Date: {context.get('project_date')}"
    write_registry_html(registry_dir, spock=1, registry_parquet="_registry.parquet", out_html="_registry.html", title="Tiberius-Spock1-Flats-GaussianSmooth", subtitle=subtitle)
   
    logger.info(f"New Gaussian Smooth: {run_id}", extra={"data_name": "flats"})
    
    return run_dir

