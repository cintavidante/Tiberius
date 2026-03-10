#####################################################
# BIAS EYEBALL FUNCTIONS
#####################################################
"""
This script contains the eyeball function. It contains the selection mode, for which
each single fits file is showed in consolse and selected for good or bad.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union
import hashlib, json, sys
import logging

import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
from astropy.io import fits

from core.product import BorgCollective
from core.instrument_utils import get_fits_image_extensions_from_config

plt.rcParams["image.origin"] =  "lower"

# Initialise the Logger. Will not be saved for Spock1
logger=logging.getLogger("Spock-1")

def _selection_hash_from_df(sel_df: pd.DataFrame) -> str:
    """
    Calculate the SHA-256 input hash from the selected files.

    Parameters
    ----------
    sel_df: pandas.DataFrame
        Selection of the input files in form of DataFrame
    
    Returns
    -------
    str
        10-character hex digit
    """
    rows = (
        sel_df[["path", "selected"]]
        .assign(
            path=lambda d: d["path"].astype(str),
            selected=lambda d: d["selected"].astype(bool),
        )
        .sort_values("path")
        .to_dict(orient="records")
    )
    s = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:10]

def _load_stage0_raw_index(spock0_product_dir:Path)-> pd.DataFrame:
    """
    Load the Spock0 raw index from its path.
    
    Parameters
    ----------
    spock0_product_dir: Path
        Path to the Spock0 raw_index.parquet file

    Returns
    -------
    pandas.DataFrame
        Raw index Table of the Spock0 folder
    """
    p = spock0_product_dir / "raw_index.parquet"
    if not p.exists():
        logger.error(f"Stage0 raw_index.parquet not found: {p}")
        sys.exit()
    df = pd.read_parquet(p)
    if "list_name" not in df.columns or "path" not in df.columns:
        logger.error("Stage0 raw_index.parquet must contain columns: list_name, path")
        sys.exit()
    return df

def _show_frame(arrs: List[np.ndarray], title: str) -> None:
    """
    Plot the single frames.
    
    Parameters
    ----------
    arrs: list[np.ndarray]
        List of Array(s), correspond to the seperate windows
    title: str
        Title of the Image Imshow Frames
    
    Returns
    -------
    None
    """
    plt.figure(figsize=(10,5))
    if len(arrs) == 1:
        a = arrs[0]
        plt.imshow(a, vmin=np.median(a)*0.99, vmax=np.median(a)*1.01, cmap="viridis")
        plt.title(title)
        plt.colorbar(label="ADU")

    else:
        for i, a in enumerate(arrs):
            plt.subplot(1, len(arrs), i+1)
            plt.imshow(a, vmin=np.median(a)*0.99, vmax=np.median(a)*1.01, cmap="viridis")
            plt.title(title)
            plt.colorbar(label="ADU")
    
    plt.tight_layout()
    plt.show(block=False)

def run_bias_eyeball(
        *,
        spock0_run_dir: Path,
        bias_list_name: str,
        instr_cfg,
        sel_dir: Path,
        cand_hash: str,
) -> Path:
    """
    Main bias eyeball runner. 
    
    Parameters
    ----------
    spock0_run_dir: Path
        Path to the Spock0 run directory.
    bias_list_name: str
        Name of the bias list wanted to compute the master bias 
    instr_cfg:
        Instrument configs from yaml file
    sel_dir: Path
        Path to the output folder for the selection parquet files
    cand_hash: str
        Input hash of the complete Bias List
    
    Returns
    -------
    Path
        Path to the file, containing the list with the selected bias frames 
    """
    df0 = _load_stage0_raw_index(spock0_run_dir)
    files = df0.loc[df0["list_name"] == bias_list_name, "path"].astype(str).tolist()
    if not files:
        logger.error(f"No files found for list_name: {bias_list_name!r}")
        sys.exit()
    nfiles = len(files)
    logger.info("Check the files and enter good or bad.")

    id_slice = slice(*instr_cfg.id_slice)
    idx_list = get_fits_image_extensions_from_config(instr_cfg, files[0])

    rows: List[Dict[str, Any]] = []

    for i, fp in enumerate(files, start=1):
        with fits.open(fp) as hdul:
            arrs = [hdul[idx].data for idx in idx_list]

        _show_frame(arrs, title=f"{i}/{len(files)}   {Path(fp).name}")
        logger.info(f"Frame {i}/{nfiles}: {Path(fp).name[id_slice]}")
        ans = input ("good/bad? [g,b]: ").strip().lower()
        logger.info(f"good/bad?: {ans}")
        plt.close("all")

        selected = (ans in {"", "g"})
        rows.append({"path": str(fp), "selected": bool(selected)})

    sel_df = pd.DataFrame(rows)
    sel_hash = _selection_hash_from_df(sel_df)

    sel_dir.mkdir(parents=True, exist_ok=True)
    sel_path = sel_dir / f"selection__cand={cand_hash}__sel={sel_hash}.parquet"

    sel_df.to_parquet(sel_path, index=False)
    logger.info(f"Stored selection: {sel_path}")
    return sel_path


