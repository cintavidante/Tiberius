##################################################
# BorgCollective Class
##################################################

"""
BorgCollective (Data-Product Manager of the Tiberius Pipeline)


Simple Overview
---------------
This module implements a "data product" abstraction used by the Tiberius pipeline.
Each processing step (Spock stage or sub-stage) produces a *product folder* that contains:
- product.h5/ or product.zarr:   xarray Dataset (scientific outputs)
- meta.json:                     provenance + parameters + hashes + lineage
- qa.parquet:                    tabular Query Analytics
- configs:                       optional archived config file(s)
- _SUCCESS:                      marker to signal the completion of the current run


Core TASKS
----------
1) Seperate Outputs
    - ds (xarray.Dataset) stores numeric arrays
    - meta.json stores provenance, paramers, hashes and lineage
    - qa.parquet stores Analytics from the run in a query-friendly parquet file (To use for example SQL)

2) Reproducibility
    - config_hash: what parameters were used in the .tcf files
    - input_hash: what input folder was used
    Enables full version control of the results

3) Lineage
    Each *product folder* points to the parent folder from previous stage (inside meta.parents)
    Enables full version control over multiple stages

4) Storage Format
    Both formats are fully API compatible: .zarr and .h5

5) Publishing
    Results are first saved inside a temporary folder, then renamed after successful compilation.
    At the end the _SUCCESS marker is published inside the run folder.

    
Dependencies
------------
- xarray
- pandas
- pydantic
- h5netcdf (if using HDF5 with xarray)
- zarr (if using Zarr with xarray)
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union
import hashlib
import json
import os

import xarray as xr
import pandas as pd

try:
    from importlib.metadata import version, PackageNotFoundError
except Exception:  # pragma: no cover
    version = None
    PackageNotFoundError = Exception

try:
    from pydantic import BaseModel, Field
except Exception:  # pragma: no cover
    BaseModel = object  # type: ignore
    Field = lambda default=None, **kwargs: default  # type: ignore


# ---------------------------------------
# Utilities
# ---------------------------------------

def utc_now_iso() -> str:
    """
    Return the current UTC time in ISO-8601 format. 
    Enables easy way to navigate through different runs.

    Parameters
    ----------
    None

    Returns
    -------
    str
        Timestamp in this format "2026-01-01T12:00:00Z"
    """
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _canonical_json(obj: Any) -> str:
    """
    Serialise an object to a canonical JSON string for hashing.
    This is used to ensure that logically same dictionaries will produce the same hash.
    - stable key order (sort_keys=True)
    - compact seperators to avoid whitespace differences
    
    Parameters
    ----------
    obj: Any
        JSON-serialisable object (dict/list)
    
    Returns
    -------
    str
        Canonical JSON string
    """
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_hex(s: Union[str, bytes]) -> str:
    """
    Compute SHA-256 hash of a string and return the hex.

    Parameters
    ----------
    s: str or bytes
        Input that will be hashed
    
    Returns
    -------
    str
        64-character hex digit
    """
    if isinstance(s, str):
        s = s.encode("utf-8")
    return hashlib.sha256(s).hexdigest()


def short_hash(hex_str: str, n: int = 10) -> str:
    """
    Shorten a hex digit.

    Parameters
    ----------
    hex_str: str
        Full hex digit
    n: int, default=10
        Number of characters to keep

    Returns
    -------
    str
        Shortened hex digit
    """
    return hex_str[:n]


def file_fingerprint(path: Union[str, Path], algo: str = "sha256", chunk_size: int = 1024 * 1024) -> str:
    """
    Compute a strong content hash of a file.
    This is slow but powerful. Reads the entire file and hashes its content. 

    Parameters
    ----------
    path: str or Path
        File Path
    algo: str, default=sha256
        Hash algorithm
    chunk_size: int, default=1024*1024
        Block size for hashing
    
    Returns
    -------
    str
        Hex digits of the content
    """
    p = Path(path)
    h = hashlib.new(algo)
    with p.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def manifest_hash_fast(paths: List[Union[str, Path]]) -> str:
    """
    Compute a fast content hash of a file metadata. Good for dev/iteration.
    This is fast but not ideal. Use it only for testing or if you are sure you will not open or modify files!
    It uses:
        - file path
        - file size
        - modification time

    Parameters
    ----------
    paths: list[str or Path]
        Input files
    
    Returns
    -------
    str
        SHA-256 hex digit of the canonical json file
    """
    rows = []
    for x in paths:
        p = Path(x)
        st = p.stat()
        rows.append({"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns})
    rows.sort(key=lambda r: r["path"])
    return sha256_hex(_canonical_json(rows))


def ensure_dir(p: Union[str, Path]) -> Path:
    """
    Ensure a directory exists.
    
    Parameters
    ----------
    p: str or Path
        Directory Path

    Returns
    -------
    Path
        Path object of the existing/built directory
    """
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def atomic_write_text(path: Union[str, Path], text: str) -> None:
    """
    Atomically write a text file.
    Writes the product to a temporary file in the same directory and then replaces the target file.
    Minimises the problem that files might be only partially written.
    
    Parameters
    ----------
    path: str or Path
        Destination directory
    text: str
        String to write
    
    Returns
    -------
    None
    """
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)

def _parse_scalar(v: str) -> Any:
    """
    Parse a scalar config value from a TCF file into a Python Type.

    Parameters
    ----------
    v: str
        String from the TCF file
    
    Returns
    -------
    Any
        Whatever Python value, Boolean, Integer, Float, ...
    """
    s = v.strip()
    low = s.lower()

    # booleans
    if low in {"true", "t", "yes", "y", "1", "on"}:
        return True
    if low in {"false", "f", "no", "n", "0", "off"}:
        return False

    # none/null
    if low in {"none", "null"}:
        return None

    # int
    try:
        if low.startswith("0") and len(low) > 1 and low[1].isdigit():
            # keep leading-zero strings as strings (e.g. IDs)
            raise ValueError
        return int(s)
    except Exception:
        pass

    # float
    try:
        return float(s)
    except Exception:
        pass

    # string
    return s


def read_tcf_like(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Reads a simple whitespace key/value config file (.tcf style).
    - ignores empty lines
    - ignores full-line comments starting with '#'
    - removes inline comments after '#'

    Parameters
    ----------
    path: str or Path
        Path to the config file
    
    Returns
    -------
    dict[str, Any]
        Parsed Parameters of the config file
    """
    d: Dict[str, Any] = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "#" in line:
            line = line.split("#", 1)[0].strip()
        parts = line.split(None, 1)
        if len(parts) == 2:
            k, v = parts
            d[k] = _parse_scalar(v)
    return d


def find_existing_run(stage_root: Union[str, Path], config_hash: str, input_hash: str) -> Union[str, None]:
    """
    Find the run which corresponds to config and input hash.
    
    Parameters
    ----------
    stage_root: str or Path
        Directory where parquet file is stored
    config_hash: str
        Config hash of the current product
    input_hash: str
        Input hash of the current product
    
    Returns
    -------
    str or None
        If not none, then the existing run_id to the hashes.
    """
    registry = stage_root / "_registry.parquet"

    if not registry.exists():
        return None

    df = pd.read_parquet(registry)

    match = df[
        (df["config_hash"] == config_hash)
        & (df["input_hash"] == input_hash)
    ]

    if match.empty:
        return None
    return match.iloc[0]["run_id"]


def next_run_id(stage_root: Path) -> str:
    """
    Find the next run_id.
    
    Parameters
    ----------
    stage_root: Path
        Directory where registry is store.

    Returns
    -------
    str
        Next run_id
    """

    registry = stage_root / "_registry.parquet"

    if not registry.exists():
        return "run_0001"
    
    df = pd.read_parquet(registry)

    if df.empty:
        return "run_0001"

    nums = (
        df["run_id"]
        .str.replace("run_", "")
        .astype(int)
    )

    n = nums.max() + 1

    return f"run_{n:04d}"


# ----------------------------
# Pydantic models
# ----------------------------

class StorageInfo(BaseModel):
    """
    Storage format for the xarray dataset.
    
    Attributes
    ----------
    format: str
        Storage format (zarr or hdf5)
    data_path: str
        Relative filename/path where the dataset is stored
    """
    format: str = Field(default="zarr", description="zarr or hdf5")
    data_path: str = Field(default="product.zarr", description="Relative path to the dataset store/file")

class ProductId(BaseModel):
    """
    Identity of a data product.
    
    Attributes
    ----------
    config_hash: str
        Hash of the stage parameters
    input_hash: str
        Hash of the stage inputs
    """
    config_hash: str
    input_hash: str

class MetaRecord(BaseModel):
    """
    Metadata stored in meta.json.
    This is the provenance file for the pipeline stages.
    Keep it stable and add fields as your pipeline matures.

    Attributes
    ----------
    spock_name: str
        Name of the current stage/sub-stage
    spock: int
        Spock stage number
    created_at: str
        UTC time in ISO format
    pipeline_version: str
        Package version of Tiberius if installed over pip
    git_commit: str or None
        Commit identifier
    product_id: ProductId
        config and input hashes
    storage: StorageInfo
        Dataset format and relative Path of the storage file
    parents: dict[str, str]
        Parent files used for this stage (only direct parent folders!)
    parameters: dict[str, Any]
        Parameters used for this specific product
    inputs_manifest_hash: str or None
        Hash of the manifest
    inputs: list[dict] or None
        Optional detailed manifest
    """
    spock_name: str
    spock: int
    created_at: str = Field(default_factory=utc_now_iso)

    pipeline_version: str = "unknown"
    git_commit: Optional[str] = None 

    product_id: ProductId
    storage: StorageInfo

    parents: Dict[str, str] = Field(default_factory=dict)

    parameters: Dict[str, Any] = Field(default_factory=dict)

    inputs_manifest_hash: Optional[str] = None
    inputs: Optional[List[Dict[str, Any]]] = None  # optional detailed manifest


# ----------------------------
# Storage Backend
# ----------------------------

@dataclass(frozen=True)
class StorageBackend:
    """
    Backend defines how to read/write xarray datasets.

    This Class resolves storage format differences, so that the pipeline can call at each time
    backend.write() or backend.read(). Currently working for data formats Zarr and HDF5.

    Parameters
    ----------
    format: str
        "zarr" or "h5"/"hdf5"/"netcdf"
    """
    format: str

    def default_data_name(self) -> str:
        """
        Return the default dataset filename for this backend.

        Returns
        -------
        str
            "product.zarr" or "product.h5"
        """
        return "product.zarr" if self.format.lower() == "zarr" else "product.h5"

    def write(self, ds: xr.Dataset, path: Union[str, Path], mode: str = "w") -> None:
        """
        Write a xarray Dataset.
        
        Parameters
        ----------
        ds: xarray.Dataset
            Dataset to write
        path: str or Path
            Output path for the Dataset
        mode: str, default="w"
            Zarr writing mode, (ignore for h5 files)
        """
        path = Path(path)
        fmt = self.format.lower()
        if fmt == "zarr":
            ds.to_zarr(path, mode=mode)
        elif fmt in {"h5", "hdf5", "netcdf"}:
            ds.to_netcdf(path, engine="h5netcdf", mode="w")
        else:
            raise ValueError(f"Unsupported storage format: {self.format}")

    def read(self, path: Union[str, Path]) -> xr.Dataset:
        """
        Read a xarray Dataset.
        
        Parameters
        ----------
        path: str or Path
            Input path
        
        Returns
        -------
        xarray.Dataset
            Opened Dataset able to be processed
        """
        path = Path(path)
        fmt = self.format.lower()
        if fmt == "zarr":
            return xr.open_zarr(path)
        elif fmt in {"h5", "hdf5", "netcdf"}:
            return xr.open_dataset(path, engine="h5netcdf")
        else:
            raise ValueError(f"Unsupported storage format: {self.format}")


# ----------------------------
# BorgCollective
# ----------------------------

class BorgCollective:
    """
    BorgCollective is the Tiberius product manager.
    Each processing step builds an BorgCollective that contains:
        - ds: xarray.Dataset
        - meta: MetaRecord with provenance/parameters/lineage/config_hashes
        - qa: pandas.DataFrame
    Then calls borg.publish(*run_dir*) to store the files
    At the end writes "_SUCCESS" to the product folder as cache marker
    """

    def __init__(
        self,
        spock_name: str,
        spock: int,
        storage_format: str = "h5",
        pipeline_pkg_name: str = "Tiberius",
        parents: Optional[Dict[str, str]] = None,
    ):
        """
        Create a new BorgCollective instance.
        
        Parameters
        ----------
        spock_name: str
            Name of the Spock Stage
        spock: int
            Spock Stage number
        storage_format: str, default="h5"
            Dataset storage backend
        pipeline_pkg_name: str
            Package name of the used pipeline
        parents: dict[str, str] or None
            Optional initial parent pointers
        """
        self.ds = xr.Dataset()

        # Determine pipeline version if installed as package
        pipeline_version = "unknown"
        if version is not None:
            try:
                pipeline_version = version(pipeline_pkg_name)
            except PackageNotFoundError:
                pipeline_version = "unknown"

        self.backend = StorageBackend(storage_format)

        # Start with placeholder hashe, they will be set afterwards after parsing the input parameter file
        product_id = ProductId(config_hash="UNSET", input_hash="UNSET")
        storage = StorageInfo(format=storage_format, data_path=self.backend.default_data_name())

        self.meta = MetaRecord(
            spock_name=spock_name,
            spock=spock,
            pipeline_version=pipeline_version,
            product_id=product_id,
            storage=storage,
            parents=parents or {},
            parameters={},
        )

        # Query Analytics table
        self.qa = pd.DataFrame()

        # Keep original config artifacts (optional)
        self._config_artifacts: Dict[str, str] = {} 

    # ----------------------------
    # Config & parameter handling
    # ----------------------------

    def add_parameters(self, params: Mapping[str, Any], *, prefix: Optional[str] = None) -> None:
        """
        Add parameter values to meta.parameters
        Attention: This will overwrite if you call this twice!

        Parameters
        ----------
        params: Mapping[str, Any]
            Parameters to add
        prefix: str or None
            If provided, will add a grouping.
            Example prefix="bias" -> meta.parameters["bias"] = Dict
        """
        if prefix:
            self.meta.parameters[prefix] = dict(params)
        else:
            # merge flat dict
            self.meta.parameters.update(dict(params))
    
    def read_tcf(self, filename: Union[str, Path], *, artifact_name: Optional[str] = None, prefix: Optional[str] = None) -> Dict[str, str]:
        """
        Read a .tcf style config file and store it into meta.parameters.

        Parameters
        ----------
        filename: str or Path
            Path to the config file
        artifact_name: str or None
            If provided, will archive the raw config file under product_dir/configs/
        prefix: str or None
            If provided, store parsed parameters under meta.parameters[prefix]
        
        Returns
        -------
        dict[str, Any]
            Parsed key/value dictionary
        """
        filename = Path(filename)
        parsed = read_tcf_like(filename)
        self.add_parameters(parsed, prefix=prefix)

        if artifact_name:
            self._config_artifacts[artifact_name] = filename.read_text(encoding="utf-8")

        return parsed

    def set_git_commit(self, git_commit: str) -> None:
        """
        Record the git commit in meta.json.
        
        Parameters
        ----------
        git_commit: str
            Git commit hash
        
        Returns
        -------
        None
        """
        self.meta.git_commit = git_commit

    # ----------------------------
    # Hashing (identity)
    # ----------------------------

    def set_config_hash_from_parameters(self, *, include_keys: Optional[List[str]] = None, exclude_keys: Optional[List[str]] = None) -> str:
        """
        Compute and store the product config_hash from meta.parameters.
        This version control identifies what exact parameters from the config file were used for the product.
        Attention! Do not use parameters like plots, or verbose.

        Parameters
        ----------
        include_keys: list[str] or None
            Optional whitelist from meta.parameters to include.
        exclude_keys: list[str] or None
            Optional blacklist from meta.parameters to exclude.

        Returns
        -------
        str
            Short hash stored into meta.product_id.config_hash
        """
        params = dict(self.meta.parameters)

        if include_keys is not None:
            params = {k: params[k] for k in include_keys if k in params}
        if exclude_keys is not None:
            for k in exclude_keys:
                params.pop(k, None)

        canonical = _canonical_json(params)
        cfg_hash = short_hash(sha256_hex(canonical), 10)
        self.meta.product_id.config_hash = cfg_hash
        return cfg_hash

    def set_input_hash_from_manifest(self, input_files: List[Union[str, Path]], *, strong: bool = False) -> str:
        """
        Compute and store the product input_hash from the list of input files.
        This version control identifies what exact files were used to start the spock stage.

        Parameters
        ----------
        input_files: list[str of Path]
            List of files used to compute the product.
        strong: bool, default=False
            Boolean to use fast manifest (from path, size, mtime) (False) or strong content hashing (True)
        
        Returns
        -------
        str
            Short hash stored into meta.product_id.input_hash
        """
        paths = [Path(p) for p in input_files]
        if strong:
            rows = []
            for p in sorted(paths, key=lambda x: str(x)):
                rows.append({"path": str(p), "sha256": file_fingerprint(p)})
            man_hash = sha256_hex(_canonical_json(rows))
            self.meta.inputs = rows
        else:
            man_hash = manifest_hash_fast(paths)
            self.meta.inputs = [{"path": str(p)} for p in sorted(paths, key=lambda x: str(x))]

        inp_hash = short_hash(man_hash, 10)
        self.meta.product_id.input_hash = inp_hash
        self.meta.inputs_manifest_hash = man_hash
        return inp_hash

    def set_input_hash_from_parents(self, parent_run_dirs: List[Path]) -> str:
        """
        Compute and set the product input hash from the parents hashes
        
        Paramters
        ---------
        parent_run_dirs
            List to the parent run directories.

        Returns
        -------
        str
            Hash stored into meta.product_id.input_hash
        """
        rows: List[Dict[str, Any]] = []
        for run_dir in sorted(parent_run_dirs, key = lambda p: str(p)):
            meta = BorgCollective.load_meta(run_dir)
            rows.append({
                "run_dir": str(run_dir),
                "spock": meta.spock,
                "spock_name": meta.spock_name,
                "config_hash": meta.product_id.config_hash,
                "input_hash": meta.product_id.input_hash,
            })
        
        man_hash = sha256_hex(_canonical_json(rows))
        inp_hash = short_hash(man_hash, 10)

        self.meta.inputs = rows
        self.meta.inputs_manifest_hash = man_hash
        self.meta.product_id.input_hash = inp_hash
        return inp_hash

    # ----------------------------
    # Lineage / Parents
    # ----------------------------

    def set_parent(self, name: str, product_path: Union[str, Path]) -> None:
        """
        Record a pointer to a parent product folder.
        
        Parameters
        ----------
        name: str
            Parent folder key e.g. "spock0", "bias"
        product_path: str of path
            Path to the parent product folder

        Returns
        -------
        None
        """
        self.meta.parents[name] = str(product_path)

    # ----------------------------
    # Dataset attributes (keep minimal!)
    # ----------------------------

    def apply_identity_attrs_to_dataset(self) -> None:
        """
        Write minimal identity into the ds.attrs.
        Full provenance stays in meta.json.

        Returns
        -------
        None
        """
        self.ds.attrs.update(
            {
                "spock_name": self.meta.spock_name,
                "spock": self.meta.spock,
                "created_at": self.meta.created_at,
                "pipeline_version": self.meta.pipeline_version,
                "config_hash": self.meta.product_id.config_hash,
                "input_hash": self.meta.product_id.input_hash,
            }
        )

    # ----------------------------
    # Query Analysis handling
    # ----------------------------

    def set_qa(self, qa: pd.DataFrame) -> None:
        """
        Set the QA dataframe for the product.

        Parameters
        ----------
        qa: pandas.DataFrame
            QA table. Stored in qa.parquet during publish()
        """
        self.qa = qa.copy()

    def add_qa_row(self, row: Mapping[str, Any]) -> None:
        """
        Append a row to the QA table.
        
        Parameters
        ----------
        row: Mapping[str, Any]
            Row that will be appended to the DataFrame
        
        Returns
        -------
        None
        """
        self.qa = pd.concat([self.qa, pd.DataFrame([dict(row)])], ignore_index=True)

    # ----------------------------
    # Writing / publishing
    # ----------------------------
    INFO_DIRNAME = "_info" # If you change this also change it at the top of this script!

    @classmethod
    def info_dir(cls, run_dir: Union[str, Path]) -> Path:
        """
        Return the path of the run folder
        
        Parameters
        ----------
        run_dir: str or Path
            run directory
        
        Returns
        -------
        Path
            Path of the run directory
        """
        return Path(run_dir) / cls.INFO_DIRNAME

    @classmethod
    def _success_path(cls, run_dir: Union[str, Path]) -> Path:
        """
        Return the path of the _SUCCESS cache marker file.
        
        Parameters
        ----------
        product_dir: str or Path
            Product folder
            
        Returns
        -------
        path
            Product_dir/_SUCCESS
        """
        return cls.info_dir(run_dir) / "_SUCCESS"
    
    def exists_and_success(self, run_dir: Union[str, Path]) -> bool:
        """
        Cache hit check. Checks if directory exists and has "_SUCCESS" marker

        Parameters
        ----------
        run_dir: str or path
            Run folder path
        
        Returns
        -------
        bool
            True if _SUCCESS exists else False
        """
        return self._success_path(run_dir).exists()

    def write_meta_json(self, run_dir: Union[str, Path]) -> Path:
        """
        Writes meta.json into the product directory.
        
        Parameters
        ----------
        run_dir: str or Path
            Run directory where meta.json should be saved
        
        Returns
        -------
        path
            Path of the written json file
        """
        info = ensure_dir(self.info_dir(run_dir))
        path = info / "meta.json"
        atomic_write_text(path, self.meta.model_dump_json(indent=2) if hasattr(self.meta, "model_dump_json") else json.dumps(self.meta.__dict__, indent=2))
        return path

    def write_qa_parquet(self, run_dir: Union[str, Path]) -> Optional[Path]:
        """
        Write qa.parquet into the product directory.
        
        Parameters
        ----------
        run_dir: str or path
            Run directory where qa.parquet should be saved
        
        Returns
        -------
        path
            Path of the written parquet file if QA is not empty
        """
        if self.qa is None or self.qa.empty:
            return None
        info = ensure_dir(self.info_dir(run_dir))
        path = info / "qa.parquet"
        self.qa.to_parquet(path, index=False)
        return path

    def write_config_artifacts(self, run_dir: Union[str, Path]) -> List[Path]:
        """
        Write archived config artifacts into run_000*/config.

        Parameters
        ----------
        run_dir: str or Path
            Run directory where the artifacts should be saved

        Returns
        -------
        list[path]
            List of paths of the written artifacts. This might be empty.
        """
        out: List[Path] = []
        if not self._config_artifacts:
            return out
        info = ensure_dir(self.info_dir(run_dir))
        cfg_dir = ensure_dir(info / "configs")
        for name, text in self._config_artifacts.items():
            p = cfg_dir / name
            atomic_write_text(p, text)
            out.append(p)
        return out

    def write_product(self, run_dir: Union[str, Path], *, mode: str = "w") -> Path:
        """
        Write the xarray to the product directory using the storage backend.

        Parameters
        ----------
        run_dir: str or Path
            Run directory
        mode: str, default="w"
            Zarr mode, ignored for HDF5
        
        Returns
        -------
        Path
            Path to the stored Dataset
        """
        run_dir = ensure_dir(run_dir)
        data_path = Path(run_dir) / self.meta.storage.data_path
        self.backend.write(self.ds, data_path, mode=mode)
        return data_path

    def publish(
        self,
        run_dir: Union[str, Path],
        *,
        overwrite: bool = False,
        write_reports: bool = False,  # placeholder: for reports
        write_product: bool = True,
        extra_writers: Optional[List] = None
    ) -> Dict[str, Optional[Path]]:
        """
        Publish a complete run directory.
        Writes to a temporary folder first and only after successful stage run, the output will be renamed and
        the "_SUCCESS" cache hit file is stored.

        Parameters
        ----------
        run_dir: str or Path
            Final run directory
        overwrite: bool, default=False
            - True: replace the existing product_dir. Should only be used for dev.
            - False: if product_dir exists, then no publish will be executed.
        write_reports: bool, default=False
            reserved for plotting, still in development
        extra_writers: list or None
            Optional list of callbacks. Use for additional files, like raw_index.parquet, previews,...

        Returns
        -------
        dict[str, Path or None]
            Path of _SUCCESS
        """
        run_dir = Path(run_dir)
        info = self.info_dir(run_dir)

        # Cache hit behavior
        if self.exists_and_success(run_dir) and not overwrite:
            return {
                "product": run_dir / self.meta.storage.data_path,
                "meta": info / "meta.json",
                "qa": info / "qa.parquet" if (info / "qa.parquet").exists() else None,
                "success": self._success_path(run_dir),
            }

        # Prepare tmp dir
        tmp_dir = run_dir.parent / (run_dir.name + "._tmp")
        if tmp_dir.exists():
            for root, dirs, files in os.walk(tmp_dir, topdown=False):
                for f in files:
                    Path(root, f).unlink(missing_ok=True)
                for d in dirs:
                    Path(root, d).rmdir()
            tmp_dir.rmdir()
        ensure_dir(tmp_dir)

        # Make identity attrs available on the dataset itself (minimal)
        self.apply_identity_attrs_to_dataset()

        # Write artifacts to tmp folders
        if write_product:
            out_product = self.write_product(tmp_dir, mode="w")
        out_meta = self.write_meta_json(tmp_dir)
        out_qa = self.write_qa_parquet(tmp_dir)
        self.write_config_artifacts(tmp_dir)

        if extra_writers:
            for fn in extra_writers:
                fn(tmp_dir)

        # Write success marker LAST in tmp
        self._success_path(tmp_dir).write_text("", encoding="utf-8")

        # Finalize: replace target
        if run_dir.exists() and overwrite:
            # remove old contents
            for root, dirs, files in os.walk(run_dir, topdown=False):
                for f in files:
                    Path(root, f).unlink(missing_ok=True)
                for d in dirs:
                    Path(root, d).rmdir()
            run_dir.rmdir()

        # atomic-ish rename (good on local filesystems)
        tmp_dir.replace(run_dir)
        info= self.info_dir(run_dir)

        return {
            "product": run_dir / self.meta.storage.data_path,
            "meta": info / "meta.json",
            "qa": info / "qa.parquet" if out_qa is not None else None,
            "success": self._success_path(run_dir),
        }

    # ----------------------------
    # Reading existing products
    # ----------------------------

    @classmethod
    def load_meta(cls, run_dir: Union[str, Path]) -> MetaRecord:
        """
        Load and parse a meta.json file.
        
        Parameters
        ----------
        product_dir: str or Path
            Product_folder that contains meta.json
        
        Returns
        -------
        MetaRecord
            Parsed metadata recored
        """
        meta_path = cls.info_dir(run_dir) / "meta.json"
        obj = json.loads(meta_path.read_text(encoding="utf-8"))
        return MetaRecord(**obj)

    @classmethod
    def open_product(cls, run_dir: Union[str, Path]) -> xr.Dataset:
        """
        Open an existing product dataset.
        It reads the meta.json file to find the dataset type and filename.

        Parameters
        ----------
        run_dir: str or Path
            Run folder
        
        Returns
        -------
        xarray.Dataset
            Opened Dataset ready for further processing
        """
        meta = cls.load_meta(run_dir)
        backend = StorageBackend(meta.storage.format)
        data_path = Path(run_dir) / meta.storage.data_path
        return backend.read(data_path)