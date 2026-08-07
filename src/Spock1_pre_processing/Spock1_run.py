####################################################################
# Spock 1: Pre-Processing - Run File
####################################################################
"""
Documentation:

Currently Working for NTT_EFOSC2 and WHT_ACAM instruments.
BIAS, FLATS, BAD-PIXEL-MASK, COSMIC RAYS

Primary outputs:
- raw_index.parquet   (one row per FITS file)
- meta.json           
- qa.parquet          (Metrics)
- _SUCCESS            (publish marker)
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

import pandas as pd
from astropy.io import fits

from core.product import BorgCollective, find_existing_run, next_run_id
from core.registry import update_registry, find_latest_product
from core.instruments import load_instrument_config
from core.html_report import write_registry_html
from core.html_report import _write_table_preview
from core.instrument_utils import get_fits_image_extensions_from_config

from Spock1_pre_processing.bias.bias_compute import create_master_bias
from Spock1_pre_processing.bias.bias_eyeball import run_bias_eyeball
from Spock1_pre_processing.flats_badpixelmask.flats_master_compute import create_master_flats
from Spock1_pre_processing.flats_badpixelmask.flats_models_compute import run_median_smooth, run_gaussian_smooth
#from Spock1_pre_processing.flats_badpixelmask.bad_pixel_mask_compute import run_bad_pixel_masks

# Initialise the Logger. Will not be saved for Spock1
logger=logging.getLogger("Spock-1")

# ============================================================
# Sub-Spock: Bias
# ============================================================

def write_master_bias_fits(tmp_dir: Path, borg: BorgCollective) -> None :
    """
    Write the master bias as .Fits file to run_id folder"
    
    Parameters
    ----------
    tmp_dir: Path
        Path of the temporary run_id folder
    borg: BorgCollective
        BorgCollective instance
    
    Returns
    -------
    None
    """
    data = borg.ds["master_bias"].values.astype("float32")
    hdu = fits.PrimaryHDU(data)

    hdu.header["PIPELINE"] = "TIBERIUS"
    hdu.header["SPOCK"] = borg.meta.spock
    hdu.header["NAME"] = borg.meta.spock_name
    hdu.header["CFG_HASH"] = borg.meta.product_id.config_hash
    hdu.header["INP_HASH"] = borg.meta.product_id.input_hash

    hdu.writeto(tmp_dir / "master_bias.fits", overwrite=True)

def write_master_flat_fits(tmp_dir: Path, borg: BorgCollective) -> None :
    """
    Write the master flat as .Fits file to run_id folder"
    
    Parameters
    ----------
    tmp_dir: Path
        Path of the temporary run_id folder
    borg: BorgCollective
        BorgCollective instance
    
    Returns
    -------
    None
    """
    data = borg.ds["master_flat"].values.astype("float32")
    hdu = fits.PrimaryHDU(data)

    hdu.header["PIPELINE"] = "TIBERIUS"
    hdu.header["SPOCK"] = borg.meta.spock
    hdu.header["NAME"] = borg.meta.spock_name
    hdu.header["CFG_HASH"] = borg.meta.product_id.config_hash
    hdu.header["INP_HASH"] = borg.meta.product_id.input_hash

    hdu.writeto(tmp_dir / "master_flat.fits", overwrite=True)

def spock1_bias_artifacts_writer(self, df: pd.DataFrame) -> None:
    """
    Helper function that returns the writer function which is accepted by BorgCollective.publish()
    
    Parameters
    ----------
    df: pandas.DataFrame
        Data Frame that stores the raw index
    
    Returns
    -------
    None
    """

    def _writer(tmp_dir: Path) -> function:
        """
        Write supplementary Data to extra_writers for BorgCollective.publish().

        Parameters
        ----------
        tmp_dir: Path
            Path to the run_id folder

        Returns
        -------
        function
            Returns function itself
        """
        write_master_bias_fits(tmp_dir, self)
        tmp_dir=tmp_dir / "_info"
        df.to_parquet(tmp_dir / "raw_index.parquet", index=False)
        _write_table_preview(df, tmp_dir / "raw_index_preview.html", title="Bias raw_index")
        if self.qa is not None and not self.qa.empty:
            _write_table_preview(self.qa, tmp_dir / "qa_preview.html", title="Bias QA")
    
    return _writer


def spock1_flats_artifacts_writer(self, df: pd.DataFrame) -> None:
    """ box width, MADS and 
    Helper function that returns the writer function which is accepted by BorgCollective.publish()
    
    Parameters
    ----------
    df: pandas.DataFrame
        Data Frame that stores the raw index
    
    Returns
    -------
    None
    """

    def _writer(tmp_dir: Path) -> function:
        """
        Write supplementary Data to extra_writers for BorgCollective.publish().

        Parameters
        ----------
        tmp_dir: Path
            Path to the run_id folder

        Returns
        -------
        function
            Returns function itself
        """
        write_master_flat_fits(tmp_dir, self)
        tmp_dir=tmp_dir / "_info"
        df.to_parquet(tmp_dir / "raw_index.parquet", index=False)
        _write_table_preview(df, tmp_dir / "raw_index_preview.html", title="Bias raw_index")
        if self.qa is not None and not self.qa.empty:
            _write_table_preview(self.qa, tmp_dir / "qa_preview.html", title="Bias QA")
    
    return _writer




def run_bias(
    parent_run_dir: Path,
    sp1_params: dict[str, Any],
    context: dict[str, str],
    storage_format: str,
    overwrite: bool,
) -> Optional[Path]:
    """
    Bias runner. Receives the Master Bias, stores all in BorgCollective and publishes the files.
    2 possible variants, Normal and Eyeball.
    
    Parameters
    ----------
    parent_run_dir: Path
        Path to the Spock0 running folder
    sp1_params: dict[str, Any]
        Spock1 Parameters from config file
    context: dict[str, str]
        Context Parameters saved in BorgCollective
    storage_format: str
        Storage Format to chose from .zarr or .h5
    overwrite: bool
        Overwrite argument. Only use for dev.
        
    Returns
    -------
    Path
        Path of the bias running folder
    """
    logger.info("Start: Running bias", extra={"data_name": "bias"})

    bias_list_name = Path(sp1_params.get("bias_list")).name

    # Load Spock0 raw index
    spock0_index = pd.read_parquet(parent_run_dir / "_info/raw_index.parquet")
    bias_files = spock0_index[spock0_index["list_name"] == bias_list_name]["path"].tolist()
    if not bias_files:
        logger.error(f"No bias files found for list {bias_list_name}", extra={"data_name":"bias"})
        sys.exit(1)

    selection_mode = sp1_params.get("bias_selection_mode", False)

    # Eyeball Mode
    selection_path: Optional[Path] = None
    selected_map = None
    if selection_mode:
        logger.info("Start: Selection mode", extra={"data_name": "bias"})
        instr_cfg = load_instrument_config(context.get("instrument_id"))

        tmp = BorgCollective(spock_name="__tmp__", spock=1, storage_format=storage_format)
        cand_hash = tmp.set_input_hash_from_manifest(bias_files, strong=False)

        selection_dir = (
            Path(sp1_params.get("outputdir_Spock1", "Spock1_pre_processing"))
            / "bias"
            / "_selections"
            / bias_list_name
        )
        selection_dir.mkdir(parents=True, exist_ok=True)

        selection_path = run_bias_eyeball(
            spock0_product_dir = parent_run_dir,
            bias_list_name = bias_list_name,
            instr_cfg=instr_cfg,
            selection_dir = selection_dir,
            cand_hash = cand_hash,
        )
        selection_df = pd.read_parquet(selection_path)
        selected_map = dict(zip(selection_df["path"].astype(str), selection_df["selected"].astype(bool)))

    # Initialise Bias Borgcollective
    borg = BorgCollective(
        spock_name="Bias",
        spock=1,
        storage_format=storage_format,
    )

    # Attach parent Spock
    borg.set_parent("spock0", str(parent_run_dir))

    # Add parameters relevant to bias only
    if selection_mode:
        borg.add_parameters({"spock1_bias": {
            "list": sp1_params.get("bias_list"),
            "selection_mode": selection_mode,
            "selection_file": str(selection_path),
        }})
    else:
        borg.add_parameters({"spock1_bias": {
            "list": sp1_params.get("bias_list"),
            "selection_mode": selection_mode,
        }})

    input_files = list(bias_files)

    if selection_path is not None and selection_path.exists():
        input_files.append(selection_path)

    # Hash computation
    borg.set_config_hash_from_parameters()
    borg.set_input_hash_from_manifest(input_files, strong=False)

    # Output location
    registry_dir = (
        Path(sp1_params.get("outputdir_Spock1", "Spock1_pre_processing"))
        / "bias"
    )
    
    # Find the run with same hashes if existing
    existing = find_existing_run(registry_dir, borg.meta.product_id.config_hash, borg.meta.product_id.input_hash)
    if existing and not overwrite:
        logger.info(f"Skipped: Master Bias already exists -> {existing}", extra={"data_name", "bias"})
        return registry_dir / existing
    elif existing and overwrite:
        run_id = existing
    else: 
        run_id = next_run_id(registry_dir) 
    borg.meta.parameters["run_id"] = run_id

    # Running directory
    run_dir = registry_dir / run_id

    # Instrument Config
    instr_cfg = load_instrument_config(context.get("instrument_id"))

    # Compute bias (fills borg.ds + borg.qa)
    create_master_bias(borg, instr_cfg, context, selection_parquet=selection_path)

    # Data frame for raw_list
    bias_index = pd.DataFrame({
        "path": bias_files,
        "list_name": bias_list_name,
    })
    if selected_map is not None:
        bias_index["selected"] = bias_index["path"].astype(str).map(lambda p: bool(selected_map.get(p, False)))

    # Publish the results to running directory
    borg.publish(run_dir, overwrite=overwrite, extra_writers=[spock1_bias_artifacts_writer(borg, bias_index)])

    # Update registry table
    update_registry(registry_dir)

    # Update the registry html file
    subtitle=f"Instrument: {context.get('instrument_id')} | Planet: {context.get('project_name')} | Date: {context.get('project_date')}"
    write_registry_html(registry_dir, spock=1, registry_parquet="_registry.parquet", out_html="_registry.html", title="Tiberius-Spock1-Bias", subtitle=subtitle)
   
    logger.info(f"New Bias: {run_id}", extra={"data_name": "bias"})
    return run_dir

# ============================================================
# Sub-spock: Flats + Bad Pixel Mask
# ============================================================

def run_flats(
    parent_run_dir: Path,
    bias_dir: Path,
    sp1_params: dict[str, Any],
    context: dict[str, str],
    storage_format: str,
    overwrite: bool,
) -> Tuple[Optional[Path], dict[str, any]]:
    """
    Flats and Bad-Pixel Mask Runner.
    
    Parameters
    ----------
    parent_run_dir: Path
        Path to the parent run directory (Spock0)
    bias_dir: Path
        Path to the bias run folder
    sp1_params: dict[str, Any]
        Dictionary containing Spock1 parameters
    context: dict[str,str]  
        Dictrionary of the context parameters of the instrument
    storage_format: str
        Format to save the output science data
    overwrite: bool
        Overwrite possibility, should only be used during dev

    Returns
    -------
    Path
        Path to the Flat outputs and 
    dict[str, any]
        Dictionary of the used extensions
    """  
    #if bool(sp1_params.get("skip_flats_bad_pixel_mask", False)):  #Maybe not needed? Depends on KECK needs a testing
    #    print("Skipping flats.")
    #    return None

    logger.info("Start: Running Flats", extra={"data_name": "flats"})

    flat_list_name = Path(sp1_params.get("flat_list")).name

    spock0_index = pd.read_parquet(parent_run_dir / "_info/raw_index.parquet")
    flat_files = spock0_index[spock0_index["list_name"] == flat_list_name]["path"].tolist()

    # Instrument Config
    instr_cfg = load_instrument_config(context.get("instrument_id"))

    test_file = flat_files[0]
    idx_list = get_fits_image_extensions_from_config(instr_cfg, test_file)
    nwin = len(idx_list)
    extensions = {'idx_list': idx_list, 'nwin': nwin}

    df = BorgCollective.open_product(bias_dir)
    master_bias = df["master_bias"]

    borg = BorgCollective(
        spock_name="Flats",
        spock=1,
        storage_format=storage_format,
    )

    borg.set_parent("spock0", str(parent_run_dir))

    borg.add_parameters({"spock1_flats": {
        "list": sp1_params.get("flat_list"),
        "saturation_limit": sp1_params.get("flat_saturation_limit"),
    }})

    borg.set_config_hash_from_parameters()
    borg.set_input_hash_from_manifest(flat_files, strong=False)

    # Output location
    registry_dir = (
        Path(sp1_params.get("outputdir_Spock1", "Spock1_pre_processing"))
        / "flats"
    )
    
    # Find the run with same hashes if existing
    existing = find_existing_run(registry_dir, borg.meta.product_id.config_hash, borg.meta.product_id.input_hash)
    if existing and not overwrite:
        logger.info(f"Skipped: Master Flat already exists -> {existing}", extra={"data_name", "flats"})
        return registry_dir / existing
    elif existing and overwrite:
        run_id = existing
    else: 
        run_id = next_run_id(registry_dir) 
    borg.meta.parameters["run_id"] = run_id

    # Running directory
    run_dir = registry_dir / run_id

    create_master_flats(borg, flat_files, master_bias, extensions, instr_cfg)

    # Data frame for raw_list
    flats_index = pd.DataFrame({
        "path": flat_files,
        "list_name": flat_list_name,
    })


    flats_index["selected"] = borg.qa["selected"]

    borg.publish(run_dir, overwrite=overwrite, extra_writers=[spock1_flats_artifacts_writer(borg, flats_index)])

    # Update registry table
    update_registry(registry_dir)

    # Update the registry html file
    subtitle=f"Instrument: {context.get('instrument_id')} | Planet: {context.get('project_name')} | Date: {context.get('project_date')}"
    write_registry_html(registry_dir, spock=1, registry_parquet="_registry.parquet", out_html="_registry.html", title="Tiberius-Spock1-Flats", subtitle=subtitle)
   
    logger.info(f"New Flats: {run_id}", extra={"data_name": "flats"})

    return run_dir, extensions
    
def run_flat_models(
    flat_run_dir: Path,
    sp1_params: dict[str, Any],
    extensions: dict[str, Any],
    context: dict[str, str],
    storage_format: str,
    overwrite: bool,
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Master Flat post processing models.


    Parameters
    ----------
    flat_run_dir: Path
        Path to the parent flat run directory
    sp1_params: dict[str, Any]
        Dictionary containing Spock1 parameters
    extensions: dict[str, Any]
        Dictionary containing window extensions (nwin, list(extensions))
    context: dict[str,str]  
        Dictrionary of the context parameters of the instrument
    storage_format: str
        Format to save the output science data
    overwrite: bool
        Overwrite possibility, should only be used during dev

    Returns
    -------
    Path
        Path to the Model output
    Path
        Path to the Model output
    """
    
    logger.info("Start: Running Flat Smoothing Models", extra={"data_name": "flats"})

    # Median smooth
    median_run_dir = run_median_smooth(flat_run_dir, extensions["nwin"], sp1_params, context, storage_format, overwrite)
    
    # Gaussing smooth
    # ## Note: Gaussian smooth not currently used, despite good perfomance, as it removes features in x as well as y
    gaussian_run_dir = run_gaussian_smooth(flat_run_dir, extensions["nwin"], sp1_params, context, storage_format, overwrite)
                     
    return median_run_dir, gaussian_run_dir

# Under construction
#def run_bad_pixel_mask(
#        parent_run_dir,
#        sp1_params,
#        context,
#        storage_format,
#        overwrite,):
#
#    return run_dir

# ============================================================
# Main Spock1 Runner (Orchestrator)
# ============================================================

def run_Spock1_pre_processing(
    config_path: str = "Spock1_pre_processing.tcf",
    parent_run_dir: Union[str, Path] = None,
    overwrite: bool = False,
) -> None:
    """
    Main Spock-1 Running function. Whole Pre-Processing Stage.

    Parameters
    ----------
    config_path: str
        Config file Path
    parent_run_dir: str or Path, default=None
        Path to the parent running folder
    overwrite: bool, default=False

    Returns
    -------
    None
    """
    logger.info("[Start] Pre-Processing.")

    if parent_run_dir is None:
        logger.info("No parent folder specified. Use latest Spock0 product.")

        # Here would be maybe better to check for Folder starting with Spock0
        parent_run_dir = find_latest_product(Path("Spock0_initialisation"))

    parent_run_dir = Path(parent_run_dir)

    # --------------------------------------------------
    # Create orchestration borg (No Storage, Only to read Spock1 Parameters)
    # --------------------------------------------------
    orchestrator = BorgCollective(
        spock_name="spock1_pre_processing",
        spock=1,
        storage_format="zarr",
    )

    orchestrator.read_tcf(
        config_path,
        artifact_name=Path(config_path).name,
        prefix="spock1",
    )

    sp1 = orchestrator.meta.parameters.get("spock1", {})

    # --------------------------------------------------
    # Inherit context from Spock0
    # --------------------------------------------------
    parent_meta = BorgCollective.load_meta(parent_run_dir)
    context = parent_meta.parameters.get("context", {})
    orchestrator.add_parameters({"context": context})
    instrument = context.get("instrument_id")
    storage_format = context.get("storage_format")

    # --------------------------------------------------
    # Run Sub-Spocks, Bias, Flats, Cosmics
    # --------------------------------------------------

    bias_dir = run_bias(
        parent_run_dir,
        sp1,
        context,
        storage_format,
        overwrite,
    )

    flat_dir, extensions = run_flats(
        parent_run_dir,
        bias_dir,
        sp1,
        context,
        storage_format,
        overwrite,
    )

    flat_model_dirs = run_flat_models(
        flat_dir,
        sp1,
        extensions,
        context,
        storage_format,
        overwrite,
    )

#    bad_pixel_dirs = run_bad_pixel_masks(
#        flat_dir,
#        sp1,
#        context,
#        storage_format,
#        overwrite,
#    )

    logger.info("[Completed] Pre-Processing.")
    logger.info(f"[Completed] Bias Product: {bias_dir}")
    #print(f"Bias product: {bias_dir}")
    #print(f"Flats product: {flats_dir}")

    return {
        "bias": bias_dir,
        "flats": flat_dir,
    }


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":

    run_Spock1_pre_processing(
        config_path="Spock1_pre_processing.tcf",
        parent_run_dir="Spock0_calib_output/config_hash=XXX/input_hash=YYY",
        overwrite=True,
    )