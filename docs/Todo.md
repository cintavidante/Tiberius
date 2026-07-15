# TODO Tiberius 2.0

## Python Version
- [x] Update 3.8 -> 3.12
- [x] Modernise Github pip installation
- [ ] Add possibility if forwarding on ssh not working!

## Structure BorgCollective - Brain of Tiberius
- [ ] Build the BorgCollective
    - [x] Add Debugging mode -> overwrite
    - [x] MetaData saved inside JSON files
    - [x] Context MetaData (telescope, instrument, date, project_name, storage_format) loaded from previous Stage
    - [x] Results saved in form of parquet files, Query Analysis and Raw_index
    - [x] Build 64 character hex digits for input parameters + files, to get version control.
    - [x] Publish function to publish results, and add a SUCCESS marker file.
    - [x] Add registry that contains all major input parameters and major outputs
    - [x] Write HTML files (help of CHAT-GPT)
    - [ ] Lineage integration, all folders that were used in the production of a run folder, to trace back to previous stages.
- [x] Function to trace back if a concrete folder exists, makes use of Hashes + Next run_id
- [x] Add Logging on master and stage layer
- [ ] Add the reports (save images)

## Spock 0 - Initialisation
- [x] Structure BorgCollective, no reports here needed (plots)
- [x] Documentation
- [x] Add context parameters to .tcf file.
- [x] Add storage format to .tcf file (zarr of h5)
- [ ] Ask James about KECK data

## Spock 1 - Pre Processing
- [ ] Documentation
- [ ] Bias
    - [x] Structure BorgCollective
    - [x] Add eyeball mode for 2 window fits files (ACAM)
    - [ ] Ask to add np.shape to raw_index parquet
    - [ ] Ask James about KECK data -> Skip bias needed?
    - [ ] Reports (images) 
- [ ] Flats & Bad-Pixel-Mask
    - [ ] Flats
        - [x] Structure BorgCollective
        - [ ] Add 2 window version
        - [ ] Add slicing for .tcf file input of list
        - [ ] Ask James about KECK data, -> Skip stage needed? If No then make one input if no flats available
        - [x] Add saturation limit to .tcf file
        - [ ] Models
            - [x] Structure BorgCollective
            - [x] Add Gaussian Model to .tcf (+ sigma-dev)
            - [x] Add Median Model to .tcf (+ box-width)
            - [ ] Add possibility from .tcf to use the box_width from the .tcf.
        - [ ] Reports
    - [ ] Bad-Pixel-Mask
        - [ ] Structure BorgCollective 
        - [x] Add Medianfilter cutoff to .tcf file
        - [ ] Reports
- [ ] Cosmics
    - [ ] Structure BorgCollective 
    - [x] Add groups
    - [x] Add ACAM compatibility
    - [x] Fix Bug with Bias
    - [ ] Ask James about KECK data


## Spock 2 - Reduced Images


## Spock 3 - Extraction 
- [x] Add different left/right sky background widths


## Spock 4 - Reduction Notebooks
- [ ] Talk to group maybe way to automize it a bit more? 

## Spock 5 - Fitting
- [x] Add Spotrod model
- [x] Fix Bug for GP

