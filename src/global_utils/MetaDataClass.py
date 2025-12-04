import xarray as xr
from datetime import datetime
from importlib.metadata import version, PackageNotFoundError


class BorgCollective:
    " A class that contains Tiberius MetaData "

    def __init__(self, spock: int):
        "Initialise an empty dataset with base attributes"
        try:
            pipeline_version = version("Tiberius")
        except PackageNotFoundError:
            pipeline_version = "unknown"

        self.ds = xr.Dataset()
        self.ds.attrs = {
            "pipeline_version": pipeline_version,
            "created": datetime.utcnow().isoformat() + "Z",
            "spock": spock
            }


    def read_tcf(self, filename: str):
        "Read .tcf files and stores them as global attributes"
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                if "#" in line:
                    line = line.split('#', 1)[0].strip()

                parts = line.split(None, 1)
                if len(parts) == 2:
                    key, value = parts
                    self.ds.attrs[key] = value
                else:
                    print(f"Warning: cannot parse line: {line}")


    def combine(self, old_ds: xr.Dataset):
        "Merge two xarrays into one xarray"
        merged = xr.merge([self.ds, old_ds])
        merged.attrs = {**self.ds.attrs, **old_ds.attrs}

        merged.attrs["pipeline_version"] = self.ds.attrs["pipeline_version"]
        merged.attrs["created"] =  self.ds.attrs["created"]
        merged.attrs["spock"] = self.ds.attrs["spock"]

        self.ds = merged

    
    def write(self, path: str):
        "Save the dataset to an .h5 file"
        self.ds.to_netcdf(path, engine='h5netcdf')
