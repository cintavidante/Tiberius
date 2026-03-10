############################################
# Spock 0: Initialisation - Compute File
############################################
"""
Sort raw FITS/FIT/FITS.Z files into logical groups based on headers.
Saves the result inside the BorgCollective for later publishement.
"""

from __future__ import annotations

import glob
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, Callable, Tuple

import pandas as pd
import xarray as xr
from astropy.io import fits

from core.product import BorgCollective
from core.html_report import _write_table_preview

# Get the correct logger
logger=logging.getLogger("Spock-0")

# -----------------------------------------------------------------------------
# Header parsing
# -----------------------------------------------------------------------------

def efosc_headers(hdr) -> Optional[Dict[str, Any]]:
    """
    Build the Dictionary List from the EFOSC2 headers

    Parameters
    ----------
    hdr
        Fits file
    
    Returns
    -------
    dict[str, Any]
        Dictionary with the header informations
    """
    try:
        grism = hdr["HIERARCH ESO INS GRIS1 NAME"]
        filt = hdr["HIERARCH ESO INS FILT1 NAME"]
        slit = hdr["HIERARCH ESO INS SLIT1 NAME"]
        obj = hdr["OBJECT"]
    except Exception:
        return None

    if obj == "SKY,FLAT":
        obj = "SKY_FLAT"
    if obj == "WAVE":
        obj = "ARC"
    if obj in {"FOCUS", "OTHER"}:
        return None

    filt_tag = "" if filt == "Free" else "_" + str(filt)

    slit_tag = str(slit)
    if slit_tag == "slit#1.0":
        slit_tag = "1arcsec"
    if slit_tag == "slit#15.0":
        slit_tag = "15arcsec"
    if slit_tag == "Special#1":
        slit_tag = "27arcsec"

    grism_tag = str(grism)
    if grism_tag == "Gr#11":
        grism_tag = "Gr11"
    if grism_tag == "Gr#13":
        grism_tag = "Gr13"

    list_name = f"{obj}_{grism_tag}_{slit_tag}{filt_tag}_list"

    return {
        "object": obj,
        "grism": grism_tag,
        "slit": slit_tag,
        "filter": "" if filt == "Free" else str(filt),
        "list_name": list_name,
    }


def acam_headers(
        hdr,
        position1: Optional[str] = None,
        window1: Optional[str] = None,
        position2: Optional[str] = None,
        window2: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Parse the headers for ACAM to receive the ACAM List dictionary.

    Parameters
    ----------
    hdr
        Fits file
    position1: str or None
        Position of the first window. Might be "1","2","3","4" or None
    window1: str or None
        Window corresponding to position1. Might be in this format "[a:b,c:d]" or None
    position2: str or None
        Position of the second window
    window2: str or None
        Window corresponding to position2

    Returns    position1: str or None
        Position of the first window. Might be "1","2","3","4" or None
    window1: str or None
        Window corresponding to position1. Might be in this format "[a:b,c:d]" or None
    position2: str or None
        Position of the second window
    window2: str or None
        Window corresponding to position2
    -------
    dict[str, Any]
        Dictionary with the header informations
    """
    try:
        # Object name
        obj = hdr["Object"]
        obj = obj.replace(" ", "_").replace("/", "_") # replace with spaces, replace forward slashes from FOCRUN

        filename = obj

        # Slit
        slit = hdr.get("ACAMSLI", "CLR")
        if slit != "CLR":
            if slit == "SLIT":
                slit_tag = "40arcsec"
            else:
                slit_tag = slit.replace(".", "p") + "arcsec"
            filename += "_" + slit_tag

        # Grism
        grism = hdr.get("ACAMDISP")
        #if grism != 'NONE':
            #list_name += '_'+grism

        # Filter 1
        wheel1 = hdr.get("ACAMWH1", "CLEAR")
        if wheel1 != "CLEAR":
            filename += "_" + wheel1

        wheel2 = hdr.get("ACAMWH2", "CLEAR")
        if wheel2 != "CLEAR":
            filename += "_" + wheel2

        # Blocking filters etc.
        mask = hdr.get("ACAMMASK")
        #if mask != 'CLR':
            #list_name += '_'+mask

        # Readout speed
        readout_speed = hdr.get("CCDSPEED", "UNKNOWN")
        filename += "_" + str(readout_speed)


    except Exception:
        return None

    winsecs: Dict[str, str] = {}
    for i in range(1, 5): 
        key = f"WINSEC{i}"
        val = hdr.get(key)
        if val and "enabled" in str(val).lower():
            winsecs[key] = str(val.split(", ")[0].strip())
        else:
            winsecs[key] = 'disabled'
    if window1 is not None and position1 is not None:
        key = f"WINSEC{position1}"
        if winsecs.get(key, "") != window1:
            filename += "_wrong_window1"

    if window2 is not None and position2 is not None:
        key = f"WINSEC{position2}"
        if winsecs.get(key, "") != window2:
            filename += "_wrong_window2"

    
    return {
        "object": obj,
        "grism": grism,
        "list_name": filename + "_list",
        "winsec1": winsecs.get("WINSEC1"),
        "winsec2": winsecs.get("WINSEC2"),
        "winsec3": winsecs.get("WINSEC3"),
        "winsec4": winsecs.get("WINSEC4"),
    }


# -----------------------------------------------------------------------------
# Index/Table building
# -----------------------------------------------------------------------------

def classify_fits_file(
    file_path: Union[str, Path],
    instrument: str,
    position1: Optional[str] = None,
    window1: Optional[str] = None,
    position2: Optional[str] = None,
    window2: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Classify the fits files. Read the headers inside parameters.

    Parameters
    ----------
    file_path: str or Path
        Specific file path
    instrument: str
        Instrument used
    window1: str or None
        Window corresponding to position1. Might be in this format "[a:b,c:d]" or None
    position2: str or None
        Position of the second window
    window2: str or None
        Window corresponding to position2

    Returns
    -------
    dict[str, Any] or None
        If existant, the dictionary containing all necessary header parameters 
    """
    file_path = Path(file_path)
    try:
        with fits.open(file_path) as f:
            hdr = f[0].header

            inst = instrument.upper()
            if inst == "NTT_EFOSC2":
                parsed = efosc_headers(hdr)
            elif inst == "WHT_ACAM":
                parsed = acam_headers(
                    hdr,
                    window1= window1,
                    position1 = position1,
                    window2 = window2,
                    position2 = position2,
                )
            else:
                raise ValueError(f"Unknown instrument: {instrument}")
            
            if parsed is None:
                return None

            rec: Dict[str, Any] = {
                "path": str(file_path),
                "instrument": inst,
                "list_name": parsed["list_name"],
            }

            # Add parsed fields
            rec.update(parsed)

            # Useful generic header fields
            rec["exptime"] = hdr.get("EXPTIME", hdr.get("EXPOSURE", None))
            rec["mjd"] = hdr.get("MJD-OBS", None)
            rec["date_obs"] = hdr.get("DATE-OBS", None)

            return rec

    except Exception:
        return None


def build_raw_index(
    input_folder: Union[str, Path],
    instrument: str,
    *,
    pattern: str,
    position1: Optional[str] = None,
    window1: Optional[str] = None,
    position2: Optional[str] = None,
    window2: Optional[str] = None,
) -> pd.DataFrame:
    """
    Creates the Table containing all Spock 0 parameters and inputs.

    Parameters
    ----------
    input_folder: str or Path
        Folder Path of the raw data input folder
    instrument: str
        Used instrument
    pattern: str
        Pattern used to browse files inside the folder. Can be <fits, fit, fits.Z>
    window1: str or None
        Window corresponding to position1. Might be in this format "[a:b,c:d]" or None
    position2: str or None
        Position of the second window
    window2: str or None
        Window corresponding to position2
    
    Returns
    -------
    pandas.DataFrame
        Table containing all header informations of the full observation
    """
    input_folder = Path(input_folder)
    all_files = sorted(glob.glob(str(input_folder / pattern)))
    if all_files == []:
        logger.warning("NO files found! Check your inputs!")
    rows: List[Dict[str, Any]] = []
    for fp in all_files:
        rec = classify_fits_file(
            fp,
            instrument,
            position1 = position1,
            window1 = window1,
            position2 = position2,
            window2 = window2,
        )
        if rec is not None:
            rows.append(rec)

    df = pd.DataFrame(rows)

    if df.empty:
        df = pd.DataFrame(
            columns=["path", "instrument", "list_name", "object", "grism", "slit", "filter", "exptime", "mjd", "date_obs"]
        )

    return df

# -----------------------------------------------------------------------------
# BorgCollective-powered runner
# -----------------------------------------------------------------------------


def run_spock0_sort(borg: BorgCollective) -> pd.DataFrame:
    """
    Sorting algorithm for Input Files inside the *input_folder*.
    Builds the DataFrame containing all important observation parameters, lists and input files.
    Saves them inside borg and returns the pandas.DataFrame
    
    Parameters
    ----------
    borg: BorgCollective
        The instance that contains all parameters of the current spock stage.
    
    Returns
    -------
    pandas.DataFrame
        Table containing all informations of Spock 0.
    """
    sp0 = borg.meta.parameters.get("spock0", {})
    spc = borg.meta.parameters.get("context", {})
    if not sp0 or not spc:
        raise ValueError("BorgCollective.meta.parameters must include a 'spock0' and 'context' dict.")

    instrument_id = str(spc["instrument_id"])
    input_dir = Path(sp0["input_dir"])
    pattern = str(sp0.get("pattern", "*.fits"))

    # For ACAM only
    position1 = sp0.get("position_1", None)
    window1 = sp0.get("window_1", None)
    position2 = sp0.get("position_2", None)
    window2 = sp0.get("window_2", None)

    # Build Stage0 table
    df = build_raw_index(
        input_dir,
        instrument_id,
        pattern=pattern,
        position1 = position1,
        window1 = window1,
        position2 = position2,
        window2 = window2,
    )

    # Query Analytics: Count the number of files inside a list and save them inside borg.
    counts = (
        df.groupby("list_name").size().reset_index(name="n_files")
        if not df.empty else pd.DataFrame(columns=["list_name", "n_files"])
    )
    borg.set_qa(counts.assign(stage="spock0_index"))

    return df


# -----------------------------------------------------------------------------
# Example standalone run
# -----------------------------------------------------------------------------

def main():
    borg = BorgCollective(stage="spock0_index", spock=0, storage_format="zarr")
    borg.add_parameters(
        {
            "spock0": {
                "telescope": "NTT",
                "instrument": "EFOSC2",
                "input_dir": "date_folder/raw",
                "pattern": "*.fits",
                "output_root": "Spock0_index",
                "write_legacy_lists": True,
                "legacy_lists_dirname": "compat_lists",
                "clobber_legacy": True,
            }
        }
    )

    out = run_spock0_sort(borg, overwrite=True)
    print(f"Stage0 finished: {out}")


if __name__ == "__main__":
    main()