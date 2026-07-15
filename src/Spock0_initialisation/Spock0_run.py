##################################################
# Spock 0: Initialisation - Run File
##################################################
"""
Documentation:

Currently Working for NTT_EFOSC2 and WHT_ACAM instruments

Primary outputs:
- raw_index.parquet   (one row per FITS file)
- meta.json           
- qa.parquet          (Metrics)
- _SUCCESS            (publish marker)
"""

from __future__ import annotations

import argparse
import logging
import sys, glob
from pathlib import Path
from typing import Dict

import pandas as pd

from core.product import BorgCollective, find_existing_run, next_run_id
from core.registry import update_registry
from core.html_report import _write_table_preview, write_registry_html
from core.logger import setup_logging

from Spock0_initialisation.sort_images import run_spock0_sort  

# Initialise the Logger. Will not be saved for Spock0
logger=logging.getLogger("Spock-0")


def _pick_storage_format(sp0_params: Dict[str, str]) -> str:
    """
    Extract storage format from the Spock0 config dict.
    Will be added to the Borg context parameters so that it will be available in further Spock stages.

    Parameters
    ----------
    sp0_params: dict[str, str]
        Dictionary of the sp0 parameters from the .tcf file
    
    Returns
    -------
    str
        Storage format as defined in config file
    """
    candidate = sp0_params.get("storage_format")

    if candidate is None:
        raise ValueError(
            "Missing required configuration parameter 'storage_format'."
            "Supported values: zarr, hdf5, h5, netcdf."
        )
    
    fmt = str(candidate).strip().lower()

    valid_formats = {"zarr", "hdf5", "h5", "netcdf"}

    if fmt not in valid_formats:
        raise ValueError(
            f"Invalid storage_format '{candidate}'. "
            f"Supported values: {', '.join(sorted(valid_formats))}."
        )

    if fmt in {"h5", "netcdf"}:
        return "hdf5"

    return fmt

def spock0_artifacts_writer(self, df: pd.DataFrame) -> None:
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
        #self.INFO_DIRNAME
        tmp_dir = Path(str(tmp_dir) + '/' +self.INFO_DIRNAME)
        df.to_parquet(tmp_dir / "raw_index.parquet", index=False)

            # html preview (human) – limit rows so it stays fast
        _write_table_preview(df.head(2000), tmp_dir / "raw_index_preview.html", "Stage0 raw_index preview (first 2000 rows)")

            # qa parquet (machine)
        counts = (
            df.groupby("list_name").size().reset_index(name="n_files")
            if not df.empty else pd.DataFrame(columns=["list_name", "n_files"])
        )
        counts.to_parquet(tmp_dir / "qa.parquet", index=False)

            # qa html preview (human)
        _write_table_preview(counts, tmp_dir / "qa_preview.html", "Stage0 QA preview (group counts)")
    
    return _writer


def run_Spock0_initialisation(config_path: str = "Spock0_initialisation.tcf", *, overwrite: bool = False) -> Path:
    """
    Main Spock0 Runner.
    Sorts the input files, saved under the input_folder defined in the config file.
    Returns the meta, parquet and _SUCCESS files.

    Parameters
    ----------
    config_path: str, default="Spock0_initialisation.tcf"
        Relative config path
    overwrite: bool, default="False"
        Overwrite argument, should only be used for dev

    Returns
    -------
    Path   
        Path of the current run_id folder
    """
    config_path = str(config_path)

    # Parse config file first once to decide storage_format before creating the real BorgCollective
    tmp = BorgCollective(spock_name="__tmp__", spock=0, storage_format="zarr")
    sp0_try = tmp.read_tcf(config_path, artifact_name=None, prefix=None)  
    storage_format = _pick_storage_format(sp0_try)

    # Create real BorgCollective using configured storage format
    borg = BorgCollective(spock_name="Initialisation", spock=0, storage_format=storage_format)

    # Read config file and store original config file as artifact
    borg.read_tcf(config_path, artifact_name=Path(config_path).name, prefix="spock0")

    # Get keys that are necessary for stage0
    sp0 = borg.meta.parameters.get("spock0", {})

    telescope = sp0.get("telescope")
    instrument = sp0.get("instrument")
    input_dir = sp0.get("inputdir_Spock0")
    output_root = sp0.get("outputdir_Spock0", "Spock0_initialisation")
    project_name = sp0.get("project_name")
    project_date = sp0.get("project_date")
    pattern = sp0.get("input_format")

    window_1 = sp0.get("window_1")
    window_2 = sp0.get("window_2")
    position_1 = sp0.get("position_1")
    position_2 = sp0.get("position_2")

    # Make sure important parameters are defined
    if telescope is None or instrument is None or input_dir is None or pattern is None:
        raise ValueError(
            "Spock0 config must define at least:\n"
            "  telescope <NTT|WHT|KECK>\n"
            "  instrument <EFOSC2|ACAM>\n"
            "  inputdir_Spock0 <relative Path>\n"
            "  input_format <fits|fit|fits.Z>\n"
        )

    # Instrument ID
    instrument_id = str(telescope) + '_' + str(instrument)
    # Pattern of the input files
    pattern = "*." + str(pattern)


    # Core parameters that will be available at each Spock Stage
    context = {
        "instrument_id" : str(instrument_id),
        "project_name" : str(project_name),
        "project_date" : str(project_date),
        "storage_format" : str(storage_format),
    }

    sp0_params = {
        "input_dir": str(input_dir),
        "output_root": str(output_root),
        "pattern": str(pattern),
    }

    # These parameters are only available for ACAM
    if instrument_id == "WHT_ACAM":
        sp0_params["position_1"] = int(position_1)
        sp0_params["window_1"] = str(window_1)
        sp0_params["position_2"] = int(position_2)
        sp0_params["window_2"] = str(window_2)
    
    # Initialise Logger, will not be stored as folder will be overwritten
    setup_logging(f"{output_root}/reduction.log")
    logger.info("[Start] Initialisation.")
    logger.info(f"Instrument used: {instrument_id}")

    # Write sp0_params and context keys back into borg meta so downstream always sees uniform config
    borg.add_parameters({"context": context})
    borg.add_parameters({"spock0": sp0_params})

    # Calculate the hashes for version control
    input_files = sorted(glob.glob(str(Path(input_dir) / pattern)))
    borg.set_config_hash_from_parameters()
    borg.set_input_hash_from_manifest(input_files, strong=False)

    # Check if run was already done
    output_root = Path(output_root)
    existing = find_existing_run(output_root, borg.meta.product_id.config_hash, borg.meta.product_id.input_hash)

    if existing and not overwrite:
        logger.info("[Skipped] Run already exists.")
        return output_root / existing
    elif existing and overwrite:
        output_root = output_root / existing
    elif not existing and Path(output_root).is_dir():
        logger.warning("ATTENTION! Only use full dataset and don't delete/change/add files.")
        logger.warning("If you want to delete/add/change files, then make a new project folder!")
        sys.exit()
    else: 
        run_id = next_run_id(output_root)
        borg.meta.parameters["run_id"] = run_id
        
    # Run Stage-0
    df = run_spock0_sort(borg)

    # Publish the results: meta.json + qa.parquet + product + raw_index.parquet + _SUCCESS
    from Spock0_initialisation.Spock0_run import spock0_artifacts_writer
    borg.publish(
        output_root,
        overwrite=overwrite,
        write_product=False,
        extra_writers=[spock0_artifacts_writer(borg, df)],
    )

    # Update Registry Table. In Spock0 will only build it once. Only should contain 1 line.
    update_registry(output_root)

    subtitle = f"Instrument: {instrument_id} | Planet: {project_name} | Date: {project_date}"
    write_registry_html(output_root, spock=0, registry_parquet="_registry.parquet", out_html="_registry.html", title="Tiberius-Spock0", subtitle=subtitle)
    logger.info("[Completed] Initialisation.")
    return output_root


def main():
    parser = argparse.ArgumentParser(description="Run Spock0 (file sorting / indexing)")
    parser.add_argument("--config", default="Spock0_save_files.tcf", help="Path to Spock0 config file")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing Spock0 product (dev mode)")
    args = parser.parse_args()

    run_Spock0_save_files(args.config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()