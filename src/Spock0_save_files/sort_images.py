#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk 

#################################################################################
# Initialisation: Sorting the files
#################################################################################

## Cinta's note:
## Adjust so that if I run this python script in date_folder/raw, the list will be 
#$ saved in date_folder/calib_files

# ---------------------------------------------------------------------------

import glob
import os
from astropy.io import fits

# ---------------------------------------------------------------------------

def efocs_headers(hdr):

    """
    A function to read EFOSC2 headers and create appropriate list filenames.
    """

    try:
        grism = hdr['HIERARCH ESO INS GRIS1 NAME']
        filt = hdr['HIERARCH ESO INS FILT1 NAME']
        slit = hdr['HIERARCH ESO INS SLIT1 NAME']
        obj = hdr['OBJECT']
    except:
        return None

    if obj == 'SKY,FLAT':
        obj = 'SKY_FLAT'

    if obj == 'WAVE':
        obj = 'ARC'

    if obj == 'FOCUS' or obj == 'OTHER':
        return None

    if filt == 'Free':
        filt = ''
    else:
        filt = '_' + filt

    if slit == 'slit#1.0':
        slit = '1arcsec'
    if slit == 'slit#15.0':
        slit = '15arcsec'
    if slit == 'Special#1':
        slit = '27arcsec'
        
    if grism == 'Gr#11':
        grism = 'Gr11'
    if grism == 'Gr#13':
        grism = 'Gr13'

    # Create the filename
    filename = obj + '_' + grism + '_' + slit + filt + '_list'

    return filename

# ---------------------------------------------------------------------------

def acam_headers(hdr, window1=None, position1=None, window2=None, position2=None):

    try:
        # Object name
        obj = hdr['Object']
        obj = obj.replace(" ", "_") # replace white space
        obj = obj.replace("/", "_") # replace forward slashes from FOCRUN

        filename = obj
        
        # Slit
        slit = hdr['ACAMSLI']
        if slit != 'CLR':
            if slit == 'SLIT':
                slit = '40'
            else:
                slit = slit.replace('.','p')
            
            slit = slit+'arcsec'
            filename += '_'+slit
        
        # Grism
        grism = hdr['ACAMDISP']
        #if grism != 'NONE':
            #list_name += '_'+grism
        
        # Filter 1
        wheel1 = hdr['ACAMWH1']
        if wheel1 != 'CLEAR':
            filename += '_'+wheel1
        
        wheel2 = hdr['ACAMWH2']
        if wheel2 != 'CLEAR':
            filename += '_'+wheel2
        
        # Blocking filters etc.
        mask = hdr['ACAMMASK']
        #if mask != 'CLR':
            #list_name += '_'+mask
        
        # Readout speed
        readout_speed = hdr['CCDSPEED']
        filename += '_'+readout_speed
    
    except:
        return None
        
    if window1 is not None:
        if hdr['WINSEC%d'%position1] != window1 +', enabled':
            filename += '_wrong_window1'
    
    if window2 is not None:
        if hdr['WINSEC%d'%position2] != window2 +', enabled':
            filename += '_wrong_window2'
    
    return filename + '_list'

# ---------------------------------------------------------------------------

def get_list(file_names, output_folder, inst=None):

    # Pick instrument
    if inst == 'EFOSC2':
        parse_header = efocs_headers
    elif inst == 'ACAM':
        parse_header = efocs_headers
    else:
        raise ValueError(f"Unknown instrument: {inst}")

    # Enumerate files
    for file_number, file in enumerate(file_names):

        # Open fits file
        f = fits.open(file)
        hdr = f[0].header

        # Create filename based on header
        filename = parse_header(hdr)

        # Check if filename is None
        if filename is None:
            f.close()
            continue
        
        # Output the file one directory above the Raw directory
        filepath = os.path.join(output_folder, filename)
        
        try:
            file_list = open(filepath,'a')
        except:
            file_list = open(filepath,'w')
        
        # Write or append to file
        file_list.write(file + ' \n')
        file_list.close()
        
        f.close()

# ---------------------------------------------------------------------------

def sort_the_images(meta):

    # Read from xarray
    attr = getattr(meta, "ds", meta).attrs
    inst = attr['instrument']
    input_dir = attr['inputdir_Spock0']
    output_dir = attr['outputdir_Spock0']

    # Current working directory
    pwd = os.getcwd()

    # Open input files
    input_folder = os.path.join(pwd, input_dir)
    output_folder = os.path.join(pwd, output_dir)

    # Get all FITS files in input folder
    all_files = sorted(glob.glob(os.path.join(input_folder, "*.fits")))

    # Make sure that one overwrites
    for f in glob.glob(os.path.join(output_dir, "*_list")):
        open(f, 'w').close()   

    # Get file lists
    get_list(all_files, output_folder, inst)

    print('List generation complete.')

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    
    import xarray as xr

    # Example usage
    ds = xr.Dataset()
    ds.attrs['instrument'] = 'EFOSC2'
    ds.attrs['inputdir'] = 'input'
    ds.attrs['outputdir'] = 'Spock0_calib_output'
    
    sort_the_images(ds)

# # Current working directory
# pwd = os.getcwd()

# # Parent working directory (one directory above)
# parent_dir = os.path.dirname(pwd)

# # Directory of list
# list_dir = os.path.join(parent_dir, "calib_files/file_lists")

# if args.clobber:
#      preexisting = [open(i,'w') for i in glob.glob(os.path.join(list_dir, "*_list"))]
#      [i.close() for i in preexisting]
# else:
# 	pass





