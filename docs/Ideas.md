# Ideas to Tiberius:

## How should the pipeline look like? First Ideas.
1) Should be NOT NESTED and CLEAR FOLDER STRUCTURE !!!!
2) Instrument parameters not hardcoded inside the pipeline -> inside maybe yaml file?
3) No console parametrers. Only very essential arguments inside the console!
4) No single script runnings but a Tiberius run script that runs calls all stages
5) Meta knowledge what we could load from previous stages
6) There should be a debugging mode !!! (overwrite argument that will be given through the complete pipeline)
7) HTML table to show metrics, input_files and parameters 
8) Logging too to see what exactly was compiled.
9) Please check what file format beside .h5 and .fits it should have? -> maybe .zarr v.2 ?
10) Add a _SUCCESS clause when a stage was finished

## BorgCollective
- Class
- Self Functions:
    - calcualte hashes for product_id = {config and input hash}
    - add xr dataset
    - add parquet files (maybe call Query Analytics)
    - publish all results to the run_directory
    - check if a run exists and it was succsssfully runned
    - add parameters from .tcf files
    - add parent run dir for version control
    - from reading files
        - open a result: open the .h5 or .zarr file
        - read meta: read the meta 
- Check if a run exists
- Find the next run_id

## Registry
- Registry as summary from all parquet files, so you can browse all runs and compare
- Update registry function -> add runs
- These into HTML

## Logging
- Logging should save the logs at stage layer and project root (all)
- Logging format time | stage number | stage name | message -> also substage maybe for flats?

## Plotting
- Should best not be inside the computation but only afterwards to fasten the processes
- Maybe clause --reports in console if savefigs should be runned
- Also maybe --only reports or --run_id AND ALSO --stage... for running only one specific stage

## No X Forwarding
- Maybe do some arguments in .tcf for this case. 
- We could also keep a second version of .tcf file for the case X Forwarding does not work
- Use it as global argument, and in case of no X11, then plt.savefig instead of plt.show inside a temporary folder that will be deleted afterwards again. If reports then written to correct folder. Only case for bias, and flat (test_median_smooth)

## Structure Overall
- run Tiberius script
- Spock0_initialisation: sort files
- Spock1_pre_processing:
    - bias (master_bias + eyeball mode)
    - flats (master_flat + smoothing functions)
    - bad-pixel-mask (MADS, STD and Median)
    - cosmic-rays (+ groups)
- Spock2_reductions: combination of results of Spock1
- Spock3_extraction
- Spock4_reduction_notebooks
- Spock5_fitting

## Folder structure
- Stage folder
    - run_0001
        - _info
            - pa_files
            - html_files
            - _SUCCESS
            - config_file
                - Stage_config_file.tcf
        - product.h5/product.zarr
    - registry.html
    - registry.parquet

## Folder Structure on Flats 
Flat folder sturcture is quite not straightforward, as it needs to be unfortunately nested.
This should actually keep the version control and at the same time help with good browsing possibilities.
- Folder Structure
    - flats
        - run_0001
            - masterflat.fits
            - info
            - ...
            - models
                - gaussian_smooth
                    - run_0001
                    - ...
                - median_smooth
                    - run_0001
                    - ...
