#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk

#############################################################################
# Master Bias
#############################################################################
"""
Computation of the Master Bias.
For 1 and 2 window datasets.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import logging
import numpy as np 
import pandas as pd
import xarray as xr
from astropy.io import fits

from core.product import BorgCollective
from core.instrument_utils import get_fits_image_extensions_from_config

# Initialise logger
logger = logging.getLogger("Spock-1")

def _load_stage0_raw_index(stage0_product_dir:Path)-> pd.DataFrame:
    p = stage0_product_dir / "_info/raw_index.parquet"
    if not p.exists():
        raise FileNotFoundError(f"Stage0 raw_index.parquet not found: {p}")
    df = pd.read_parquet(p)
    if "list_name" not in df.columns or "path" not in df.columns:
        raise ValueError("Stage0 raw_index.parquet must contain columns: list_name, path")
    return df

def _resolve_bias_files(
        stage0_product_dir: Path,
        bias_list_name: str,
        selection_parquet: Optional[Path] = None,
) -> Tuple[List[str], Optional[pd.DataFrame], Optional[Dict[str, bool]]]:
    """
    Returns selected bias paths and selection dataframe if provided
    """
    df0 = _load_stage0_raw_index(stage0_product_dir)

    all_bias = df0.loc[df0["list_name"] == bias_list_name, "path"].astype(str).tolist()
    if not all_bias:
        raise ValueError(f"No files found in Stage0 raw_index for list_name={bias_list_name!r}")
    
    if selection_parquet is None:
        return all_bias, None, None
    
    sel = pd.read_parquet(selection_parquet)
    if "path" not in sel.columns or "selected" not in sel.columns:
        raise ValueError("selection.parquet must contain columns: path, selected")
    
    selected_map = dict(zip(sel["path"].astype(str), sel["selected"].astype(bool)))
    chosen = [p for p in all_bias if bool(selected_map.get(p, False))]

    if not chosen:
        raise ValueError("Selection file resulted in 0 selected frames. Check selection.parquet")
    
    return chosen, sel, selected_map


def _frame_stats(arr: np.ndarray) -> Dict[str, float]:
    """
    Compute the Mean, Median, Var, Std or the frame.
    
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
        "median": float(np.median(arr)),
        "var": float(np.var(arr)),
        "std": float(np.std(arr)),
    }


def combine_biases_1window(
        bias_files: List[str],
        list_name: str,
        idx_list: List[int],
        id_slice: slice,
        sel_mode: bool,
        sel_map: Dict[str,bool],
        sel_file: str,
) -> Tuple[List[np.ndarray], List[Dict[str,Any]]]:
    """
    Compute the master bias of 1 windowed datasets.
    
    Parameters
    ----------
    bias_files: list[str]
        List of the used full bias files
    list_name: str
        Name of the used bias list
    idx_list: list[int]
        Extension of the fits file with the science data
    id_slice: slice
        Slice for shortening the file names.
    sel_mode: bool
        True if selection mode is used (eyeball function)
    sel_map: dict[str, bool] or {}
        Selection map of the bias files.
    sel_file: str
        Path to the selection file
    
    Returns
    -------
    List[numpy.ndarray]
        Master Bias
    List[Dict[str,Any]]
        List of Dictionaries containing single frame parameters
    """
    per_frame_rows: List[Dict[str, Any]] = []
    stack: List[np.ndarray] = []
    idx = idx_list[0]

    for i, fp in enumerate(bias_files, start=1):
        with fits.open(fp) as hdul:
            arr = hdul[idx].data
        stack.append(arr)

        frame_id = str(fp)[id_slice]
        row = {"path": str(fp), "frame_id": frame_id, **_frame_stats(arr), "selection_mode": sel_mode, "list_name": list_name}
        if sel_mode:
            row["selection_file"] = sel_file
            row["selected"] = bool(sel_map.get(str(fp), False))
        else:
            row["selected"] = True
        per_frame_rows.append(row)
        
    master = np.median(np.stack(stack, axis=0), axis=0)
    return master, per_frame_rows

def combine_biases_2window(
        bias_files: List[str],
        list_name: str,
        idx_list: List[int],
        id_slice: slice,
        sel_mode: bool,
        sel_map: Dict[str,bool],
        sel_file: str,
) -> Tuple[List[List[np.ndarray]], Dict[str,Any]]:
    """
    Compute 2 master bias of 2 windowed datasets. (For ACAM 2 windowed mode)
    
    Parameters
    ----------
    bias_files: list[str]
        List of the used full bias files
    list_name: str
        Name of the used bias list
    idx_list: list[int]
        Extension of the fits file with the science data
    id_slice: slice
        Slice for shortening the file names.
    sel_mode: bool
        True if selection mode is used (eyeball function)
    sel_map: dict[str, bool] or {}
        Selection map of the bias files.
    sel_file: str
        Path to the selection file
    
    Returns
    -------
    List[List[numpy.ndarray]]
        Master Bias in Form (nwin, N, N)
    List[Dict[str,Any]]
        List of Dictionaries containing single frame parameters for 2 windows
    """
    per_frame_rows: List[Dict[str, Any]] = []
    stacks: List[List[np.ndarray]] = [[] for _ in idx_list]

    for i, fp in enumerate(bias_files, start=1):
        with fits.open(fp) as hdul:
            frame_id = str(fp)[id_slice]
            for w, idx in enumerate(idx_list, start=1):
                arr = hdul[idx].data
                stacks[w-1].append(arr)
                row = {"path": str(fp), "frame_id": frame_id, "n_window":w, **_frame_stats(arr), "selection_mode": sel_mode, "selection_file": sel_file, "list_name": list_name}
                if sel_mode:
                    row["selection_file"] = sel_file
                    row["selected"] = bool(sel_map.get(str(fp), False))
                else:
                    row["selected"] = True

                per_frame_rows.append(row)

    master = np.stack([np.median(np.stack(s, axis=0), axis=0) for s in stacks], axis=0)

    return master, per_frame_rows

def create_master_bias(
        borg: BorgCollective,
        instr_cfg,
        context: Dict[str, str],
        selection_parquet: Optional[Union[str, Path]] = None, 
) -> None:
    """
    Main runner script for Master Bias.

    Parameters
    ----------
    borg: BorgCollective
        BorgCollective Instance
    instr_cfg: InstrumentConfig
        Instrument Config containing telescope specific parameters
    context: Dict[str, str]
        Context Parameters of the observation
    selection_parquet: Optional[Union[str,Path]], default=None
        If selection_mode, Path to the selection parquet file

    Returns
    -------
    None
    """
    if "spock0" not in borg.meta.parents:
        raise ValueError("borg.meta.parents['spock0'] must point to Stage0 product_dir")
    stage0_dir = Path(borg.meta.parents["spock0"])

    selection_path = Path(selection_parquet) if selection_parquet is not None else None
    bias_list_name = borg.meta.parameters['spock1_bias']['list']
    bias_files, sel_df, selected_map = _resolve_bias_files(stage0_dir, bias_list_name, selection_path)

    test_file = bias_files[0]
    idx_list = get_fits_image_extensions_from_config(instr_cfg, test_file)
    nwin = len(idx_list)

    id_slice = slice(*instr_cfg.id_slice)
    selection_mode = bool(selection_path is not None)
    selection_file = str(selection_path)
    selected_map = selected_map or {}

    if nwin == 1:
        master_bias, per_frame_rows = combine_biases_1window(bias_files, bias_list_name, idx_list, id_slice, selection_mode, selected_map, selection_file)
        borg.ds = xr.Dataset(
            {"master_bias": xr.DataArray(master_bias.astype(np.float32), dims=("y", "x"))}
        )

    else:
        master_bias, per_frame_rows = combine_biases_2window(bias_files, bias_list_name, selected_map, idx_list, id_slice, selection_mode, selection_file)
        borg.ds = xr.Dataset(
            {"master_bias": xr.DataArray(master_bias.astype(np.float32), dims=("window", "y", "x"))}
        )

    
    qa = pd.DataFrame(per_frame_rows)
    summary = {
        "path": "__Master Bias__",
        "frame_id": "",
        "n_files": int(len(bias_files)),
        "var": float(np.var(master_bias)),
        "mean": float(np.mean(master_bias)),
        "median": float(np.median(master_bias)),
        "std": float(np.std(master_bias)),
        "list_name": bias_list_name,
        "selection_mode": selection_mode,
    }
    if selection_mode:
        summary["selection_file"] = selection_file
    if context.get("instrument_id") == "WHT_ACAM":
        summary["n_windows"] = int(nwin)

    qa = pd.concat([qa, pd.DataFrame([summary])], ignore_index=True)
    borg.set_qa(qa)


#### MISSING main()
    