import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------

def read_tcf(ds, filename: str):
        "Read .tcf files and stores them as global attributes"
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                # if "#" in line:
                #     line = line.split('#', 1)[0].strip()

                parts = line.split(None, 1)
                if len(parts) == 2:
                    key, value = parts
                    ds.attrs[key] = value
                else:
                    print(f"Warning: cannot parse line: {line}")

# ---------------------------------------------------------------------------

if __name__ == "__main__":

    # Create a dummy xarray dataset
    ds = xr.Dataset()

    # ds.attrs['instrument'] = 'EFOSC2'  # or 'ACAM'
    # ds.attrs['telescope'] = 'NTT'      # or 'WHT'
    # ds.attrs['project'] = 'project_folder'
    # ds.attrs['path'] = '/Users/cividante/Tiberius/project_folder/stage0/HAT-P-65_Free_27arcsec_V#641_list'

    read_tcf(ds, 'test.tcf')

    print(ds.attrs['a_path_somewhere'])

    science_list = np.loadtxt(ds.attrs['a_path_somewhere'], dtype=str)

    print(science_list)