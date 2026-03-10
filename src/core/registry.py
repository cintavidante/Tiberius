#######################################################################
# Registry
# This is a parquet file containing the basic informations of all runs
#######################################################################

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union
import json
import pandas as pd

def update_registry(spock_root: Union[str, Path]) -> None:
    """
    Update the registry file

    Parameters
    ----------
    spock_root: str or Path
        Path to the stage root

    Returns
    -------
    None
    """
    stage_root = Path(stage_root)

    rows = []

    for meta_path in spock_root.glob("**/_info/meta.json"):
        run_dir = meta_path.parent.parent

        # only successful products will be included
        if not (run_dir / "_info/_SUCCESS").exists():
            continue

        meta = json.loads(meta_path.read_text())
        
        rows.append({
            "run_id": meta.get("run_id"),
            "product_dir": str(run_dir),
            "spock": meta.get("spock"),
            "spock_name": meta.get("spock_name"),
            "created_at": meta.get("created_at"),
            "config_hash": meta.get("product_id", {}).get("config_hash"),
            "input_hash": meta.get("product_id", {}).get("input_hash"),
            **flatten_parameters(meta.get("parameters", {})),
        })

    df = pd.DataFrame(rows)

    registry_path = spock_root / "_registry.parquet"
    df.to_parquet(registry_path, index=False)


def flatten_parameters(params: Any):
    """
    Flatten the parameters, to be able to be easyly browsed (for example with SQL)
    
    Parameters
    ----------
    params: Any

    Returns
    -------
    None
    """
    out = {}

    def recurse(prefix, d):
        for k, v in d.items():
            key = f"{prefix}__{k}" if prefix else k
            if isinstance(v, dict):
                recurse(key, v)
            else:
                out[key] = v

    recurse("", params)
    return out


def find_latest_product(spock_root: Union[str,Path]) -> Path:
    """
    Find the latest product that was runned.
    
    Parameters
    ----------
    spock_root: Path or str
        Path to the stage root folder

    Returns
    -------
    None
    """
    registry = Path(spock_root) / "_registry.parquet"
    df = pd.read_parquet(registry)
    df = df.sort_values("created_at", ascending=False)
    return Path(df.iloc[0]["product_dir"])