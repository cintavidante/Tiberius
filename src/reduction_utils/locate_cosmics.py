#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk

import argparse
import pickle
import os
import copy
import glob
import shutil

import numpy as np
import xarray as xr
# import imageio.v3 as iio
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.ndimage import median_filter
from datetime import datetime, timezone

# Prevent matplotlib plotting frames upside down
plt.rcParams['image.origin'] = 'lower'

# ------------------------------------------------------------------------------------------

def save_to_xarray(array, array_dir, name, pixelclip, bias_path):

    # Make a xarray
    bad_pixels_da = xr.DataArray(array, name=name)
        
    # Add attributes (metadata)
    bad_pixels_da.attrs['locate_cosmic_sigma'] = pixelclip
    bad_pixels_da.attrs['master_bias'] = bias_path
    bad_pixels_da.attrs['date'] = datetime.now(timezone.utc).isoformat()

    # Output as h5 
    bad_pixels_da.to_netcdf(array_dir + f"/{name}.h5", engine="h5netcdf")

# ---------------------------------------------

def locate_zeroes(nframes, plot_neg_cosmic_map, plot_neg_med_map, plot_dir, showfig, savefig):

    neg_plot_dir = os.path.join(plot_dir, 'negative_median')
    os.makedirs(neg_plot_dir, exist_ok=True)

    for i in range(nframes):

        if i < 9:
            n = '0{}'.format(i+1)
        else:
            n = i + 1

        plt.figure(figsize=[15, 10])
        plt.imshow(plot_neg_cosmic_map[i], cmap='viridis', origin='lower')
        plt.colorbar()
        plt.title("Negative median values for cosmic frame {}".format(i+1))
        plt.xlabel("Pixel column")
        plt.ylabel("Pixel row")
        
        # Save the plot
        filename = os.path.join(neg_plot_dir, 'cosmic_negative_median_frame_{}.png'.format(n))

        if savefig:
            plt.savefig(filename)
        if showfig:
            plt.show(block=False)
            plt.pause(1e-6)

        plt.close()


        # -------

        plt.figure(figsize=[15, 10])
        plt.imshow(plot_neg_med_map[i], cmap='viridis', origin='lower')
        plt.colorbar()
        plt.title("Negative median values in all frame {}".format(i+1))
        plt.xlabel("Pixel column")
        plt.ylabel("Pixel row")
        
        # Save the plot
        filename = os.path.join(neg_plot_dir, 'negative_median_frame_{}.png'.format(n))

        if savefig:
            plt.savefig(filename)
        if showfig:
            plt.show(block=False)
            plt.pause(1e-6)

        plt.close()

# ---------------------------------------------

def locate_bad_frames(image_data, pixel_row, pixel_col, cut_off_arr, plot_path,
                      showfig=False, savefig=False):
    """The function that locates the frames/pixels where cosmics are located.

    Inputs:
    image_data - the image data, dimensions of nframes x nrows x ncols
    pixel_row - the row number of the pixel under consideration
    pixel_col - the column number of the pixel under consideration
    cut_off - the sigma cut off / outlier rejection threshold
    showfig - True/False - plot the outlier identififcation for this pixel?

    Returns:
    bad_frames - the array of frames for which this pixel is an outlier"""

    nframes, nrows, ncols = image_data.shape
    pixel = image_data[:,pixel_row,pixel_col].astype(float) # make sure that the science data is correctly defined as floats
    median = median_filter(pixel,3) # take a running median across 3 frames for each pixel

    # deal with edge effects of the running median
    median[0] = np.median((median[0],median[1]))
    median[-1] = np.median((median[-1],median[-2]))

    residuals = pixel - median # calulate residuals
    good_frames = ((residuals <= cut_off_arr * np.nanstd(residuals)) & 
                    (residuals >= -cut_off_arr * np.nanstd(residuals))) # locate the good frames based on residuals array

    bad_frames = ~good_frames # flip the sign to find the outliers

    # ignore the nans (saturated pixels)
    bad_frames[~np.isfinite(residuals)] = False

    if showfig or savefig: # plot output

        plt.figure()
        plt.subplot(211)
        plt.plot(pixel,label="Pixel value")
        plt.plot(median,label="Running median")
        plt.title("Pixel [%d,%d]"%(pixel_row,pixel_col))
        plt.ylabel("Counts (ADU)")
        plt.plot(np.arange(nframes)[bad_frames],pixel[bad_frames],"rx",label="Flagged outliers")
        plt.legend()

        plt.subplot(212)
        plt.ylabel("Residuals")
        plt.xlabel("Frame")
        plt.plot(residuals)
        plt.axhline(cut_off_arr*np.nanstd(residuals),ls='--',color='k')
        plt.axhline(cut_off_arr*-np.nanstd(residuals),ls='--',color='k',label="cut-off")
        plt.plot(np.arange(nframes)[bad_frames],residuals[bad_frames],"rx",label="Flagged outliers")
        plt.legend()

        if savefig:
            plt.savefig(plot_path + '/locate_cosmic.png')
        if showfig:
            plt.show()


    return bad_frames, median

# ---------------------------------------------

def plot_cosmic_frames(cosmic_pixels, cut_off_arr, plotdir, showfig=False, savefig=False):

    """A function that plots all cosmics frames"""
    plt.figure(figsize=[15, 12])

    for i, c in enumerate(cosmic_pixels):

        # Calculate percentage of cosmic pixels
        ncosmics = np.count_nonzero(c)
        n_pixels = c.shape[0] * c.shape[1]
        cosmics = ncosmics / n_pixels * 100

        # Plot figures
        plt.imshow(c, cmap='viridis', aspect="auto", interpolation='none')
        plt.colorbar()
        
        plt.title("Frame {}; Percentage: {:.3f}%, Sigma: {}".format(i+1, cosmics, cut_off_arr[i]),
                  fontsize=16)
        plt.xlabel("Pixel column")
        plt.xlabel("Pixel row")
        
        if i < 9:
            n = '0{}'.format(i+1)
        else:
            n = i + 1

        filename = os.path.join(plotdir, 'cosmic_frames_{}.png'.format(n))

        if savefig:
            plt.savefig(filename)
        if showfig:
            plt.show(block=False)
            plt.pause(1e-6)
        plt.clf()
    
    return 

# ---------------------------------------------

def check_cosmic_frames(cosmic_pixels, frame_cut_off, plot_dir, showfig=False, savefig=False):
    """A function that plots and optionally resets cosmic pixels for frames where a disproportionate number of pixels have been flagged as cosmics.

    Inputs:
    cosmic_pixels - the array of all cosmic flagged pixels, dimensions of nframes x nrows x ncols

    Returns:
    cosmic_pixels - the new array of all cosmic flagged pixels, taking into the account the user-defined reset frame masks"""

    # Make plotting path
    plot_path = os.path.join(plot_dir, "too_many_cosmic")

    # Create folder
    os.makedirs(plot_path, exist_ok=True)
    
    nframes, nrows, ncols = cosmic_pixels.shape

    ncosmics = []

    for i,c in enumerate(cosmic_pixels):
        ncosmics.append(len(np.where(c)[0]))

    median_cosmics = np.nanmedian(ncosmics)
    if median_cosmics == 0:
        median_cosmics = 1

    print("Median number of cosmics per frame = %d (%.3f%%)"%(median_cosmics,100*median_cosmics/(nrows*ncols)))

    incorrectly_flagged_cosmics = []

    for i, c in enumerate(cosmic_pixels):

        ncosmics = len(np.where(c==1)[0])

        if ncosmics > frame_cut_off*median_cosmics:

            print("Integration %d has %.2fX the median number of cosmics, somethings up"%(i,ncosmics/median_cosmics))

            plt.figure()
            plt.imshow(c,cmap='Greys', interpolation='none',aspect="auto")
            incorrectly_flagged_cosmics.append(i)
            plt.title("Integration %d"%i)
            plt.ylabel("Pixel row")
            plt.xlabel("Pixel column")

            if i < 9:
                n = '0{}'.format(i+1)
            else:
                n = i + 1

            if savefig:
                plt.savefig(plot_path + f'/suspicious_cosmic_frames_{n}.png')

            if showfig:
                plt.show(block=False)

            plt.clf()

            reset_mask = input("Reset mask for integration %d? [y/n]: "%i)
            if reset_mask == "y":
                print("...resetting mask\n")
                cosmic_pixels[i] = np.zeros_like(c)

    return cosmic_pixels

# ---------------------------------------------

def replace_cosmics(cosmic_pixels, medians, science_list, nints, cut_off_name, instrument, 
                    cleaned_direc=os.getcwd(), master_bias=None):

    """
    Function to replace cosmic pixels with median values in the science fits files.

    Parameters:
    ----------
    cosmic_pixels: numpy array (nframes x nrows x ncols)
        array of location where cosmics are located 
    medians : numpy array (nframes x nrows x ncols)
        array of median values for each pixel
    science_list: list of str
        science fits file names
    nints   : numpy array  
    jwst    : bool   
    
    """

    # # Make an overwritten directory
    # if not os.path.exists(cleaned_direc):
    #     os.makedirs(cleaned_direc)
    # else:
    #     shutil.rmtree(cleaned_direc)
    #     os.makedirs(cleaned_direc)

    # # Check if the directory is made
    # if os.path.isdir(cleaned_direc):
    #     print ('Directory is made...')

    # Gets number of frames, rows, nd columns
    nframes, nrows, ncols = cosmic_pixels.shape
    total_nints = nints[-1]

    if instrument == 'JWST':

        for i,c in enumerate(cosmic_pixels):

            jwst_fits_counter = np.digitize(i,nints)
            if i == 0 or i in nints:
                fits_file = fits.open(science_list[jwst_fits_counter],memmap=False)
                new_fits_file = copy.deepcopy(fits_file)
                filename = science_list[jwst_fits_counter].split("/")[-1]

            if jwst_fits_counter > 0:
                jwst_index_counter = i-nints[jwst_fits_counter]
            else:
                jwst_index_counter = i

            print("Cleaning integration %d, %s"%(i,filename))

            for row in range(nrows):
                new_fits_file["SCI"].data[jwst_index_counter][row][c[row]] = medians[i][row][c[row]]

            if i in nints-1:
                fits_file.close()
                print("Saving cosmic_cleaned_fits/%s"%(filename))
                file_path = os.path.join(cleaned_direc, filename)
                new_fits_file.writeto(file_path, overwrite=True)
        
        return 
    
    medians = medians.astype('uint16')  
    # array_loop = [1, 3, 5, 7, 9]

    for i in range(nframes):

        # Make copies of fits files
        f = fits.open(science_list[i])
        f_new = copy.deepcopy(f)

        f_new[0].data = f_new[0].data - master_bias

        # Filename from science list
        filename = science_list[i].split("/")[-1]

        print("Cleaning frame {}".format(i+1))

        for row in range(nrows):
            
            f_new[0].data[row][cosmic_pixels[i][row]] = medians[i][row][cosmic_pixels[i][row]]

        file_path = os.path.join(cleaned_direc, filename)

        f_new.writeto(file_path, overwrite=True)
        f.close()
    
    # new_list = input("Do you want to overwrite the cleaned science list for this pixel clip group? [y/n]: ")

    # if new_list == "y":
    
    #     # Make a new file list for the cleaned cosmic fits
    #     print('Making new list for cleaned cosmic fits...')

    #     # Read the global path for all cleaned cosmic fits
    #     file_path = os.path.join(cleaned_direc, '*.fits')
    #     all_files = sorted(glob.glob(file_path))

    #     # Write a new list. This makes sure it's overwriting new files!!!
    #     with open(args.sciencelist + '_cleaned_' + cut_off_name, 'w') as new_list:

    #         for files in all_files:
    #             new_list.write(files+' \n')

    return

# ------------------------------------------------------------------------------------------

def locate_and_correct_cosmic(meta):

    # Read meta database
    instrument = meta.attrs["instrument"]
    sciencelist = meta.attrs["science_list"]
    inputdir = meta.attrs["inputdir_Spock1"]
    save_dir = meta.attrs["outputdir_Spock1"]
    pixelclip = meta.attrs["locate_cosmics_pixelclip"]
    frameclip = meta.attrs["locate_cosmics_frameclip"]
    showfig = meta.attrs["locate_cosmics_showfig"]
    savefig = meta.attrs["locate_cosmics_savefig"]
    cleaned = meta.attrs["locate_cosmics_cleanedfits"]
    badpixel = meta.attrs["locate_cosmics_badpixelmask"]
    groups = meta.attrs["locate_cosmics_groups"]
    gifs = meta.attrs["locate_cosmics_gifs"]

    if badpixel:
        badpixelpath = meta.attrs["bad_pixel_mask"]
        badpixellist = os.path.join(save_dir, "bad_pixel_mask", badpixelpath)

    if groups:
        grouplist = meta.attrs["locate_cosmics_groups"]

    # Loading data
    # ----------------------------------------

    # Load in the list of science file names
    sciencepath = os.path.join(inputdir, sciencelist)
    science_list = np.atleast_1d(np.loadtxt(sciencepath, dtype=str))

    # Check instrument 
    if instrument != 'EFOSC2' and instrument != 'ACAM' and instrument != 'JWST':
        raise ValueError('Currently only set up to deal with EFOSC2/ACAM/JWST')

    # Find master bias
    # Here I do not allow the option to not reduce the image with bias beforehand
    bias_path = os.path.join(save_dir, 'bias', 'master_bias.fits')

    # If bias is not there, raise error
    if not os.path.isfile(bias_path):
        raise FileNotFoundError(f"Bias file not found at: {bias_path}")

    # Load master bias
    master_bias = fits.open(bias_path)[0].data

    # Load in the science data
    print("Loading in data...")
    data = []
    nints = []

    for s in science_list:

        f = fits.open(s,memmap=False)

        if instrument == 'JWST':
            data.append(f["SCI"].data)
            nints.append(f["SCI"].data.shape[0])
        else:
            data.append(f[0].data - master_bias)
            nints.append(f[0].data.shape[0])

        f.close()
    
    nints = np.cumsum(nints)

    # Define data depending if it's JWST data or not 
    if instrument == 'JWST':
        data = np.vstack(data)
    else:
        data = np.array(data)

    # Load bad pixel mask if provided
    if badpixel:

        # Read h5 xarray
        bpm = xr.open_dataset(badpixellist, engine="h5netcdf")

        # Convert xarray to numpy array
        mask = bpm[list(bpm.data_vars)[0]].values

        # Mask the data
        data[:,mask] = np.nan
    
    # Definitons
    # ----------------------------------------

    # Define the cosmic pixel flagged array, initially as an array of zeros matching the dimensions of the input data
    cosmic_pixels = np.zeros_like(data)
    nframes, nrows, ncols = data.shape

    # Define the sigma cut off
    cut_off = pixelclip

    # Make an array of pixel_clip
    cut_off_arr = np.full(nframes, cut_off, dtype=float)

    # print(cut_off_arr)

    # If n_group > 1, then we will ask which frames and which sigma values to use
    # if args.n_group > 1:

    #     for m in range(args.n_group -1):

    #         # Ask for sigma clip
    #         cut_val = float(input("Enter sigma value: "))
    #         frame_list = input("Enter frames number separated by space: ")
    #         frame_list = list(map(int, frame_list.split()))
    #         frame_list = np.array(frame_list) - 1  # Convert to zero-based index

    #         cut_off_arr[frame_list] = cut_val

    # Find unique values in sigma / cut_off
    sigmas = np.unique(cut_off_arr)

    # Join them the unique sigmas to be a string
    sigma_name = "_".join([str(s) for s in sigmas])

    # Make directory for different sigma 


    # Set up working directory
    # -----------------------------------------

    # Define directory for outputs
    cosmic_dir = os.path.join(save_dir, 'cosmic_rays')
    sigma_dir = os.path.join(cosmic_dir, f'sigma_{sigma_name}')
    plot_output = os.path.join(sigma_dir, 'output')
    gif_output = os.path.join(plot_output, "all_frames")
    cleaned_dir = os.path.join(save_dir, 'cleaned_images')

    # Make the directory
    os.makedirs(cosmic_dir, exist_ok=True)
    os.makedirs(sigma_dir, exist_ok=True)
    os.makedirs(plot_output, exist_ok=True)
    os.makedirs(gif_output, exist_ok=True)
    os.makedirs(cleaned_dir, exist_ok=True)

    # Make data arrays for cosmic pixels and median values
    # ----------------------------------------
    
    # Loop through all frames and pixels
    cosmic_pixels = np.zeros_like(data)
    median_values = np.zeros_like(data)

    # Calculate median for each row and find cosmic pixels and median values
    for row in range(nrows):

        print("Calculating medians for row %d of %d"%(row,nrows))

        # Loop through all columns to find cosmic pixel
        for col in range(ncols):

            # Get bad frames and medians
            bad_frames, medians = locate_bad_frames(data, row, col, cut_off_arr, plot_output, 
                                                    showfig=False, savefig=False)

            # Append to the cosmic pixels and median values arrays
            cosmic_pixels[:,row,col][bad_frames] = 1
            median_values[:,row,col] = medians

    # Change the cosmic pixels array to boolean type
    cosmic_pixels = cosmic_pixels.astype(bool)

    # When there is possibility of negative median values
    # if args.bias_frame is not None:

    # Negative median pixels for cosmic pixels
    neg_cosmic_map = cosmic_pixels & (median_values < 0)
    # plot_neg_cosmic_map = np.where(neg_cosmic_map, 1, 0)

    # Negative median in general
    all_array = np.full_like(median_values, True, dtype=bool)
    neg_med_map = all_array & (median_values < 0)
    # plot_neg_med_map = np.where(neg_med_map, 1, 0) 

    # Change median values that is negative to zero
    median_values[median_values < 0] = 0

    # Outputs
    # ----------------------------------------

    if showfig or savefig:

        print("\nPlotting all cosmic-masked pixels...\n")
        plot_cosmic_frames(cosmic_pixels, cut_off_arr, gif_output, showfig, savefig)

        print("\nPlotting frames with high number of cosmics...\n")
        cosmic_pixels = check_cosmic_frames(cosmic_pixels, frameclip, plot_output, showfig, savefig)

    # Save the cosmic masks
    name_file = f'cosmic_pixel_mask_sigma{sigma_name}'
    save_to_xarray(cosmic_pixels, sigma_dir, name_file, sigma_name, bias_path)
    
    
    # For possibility of negative median values
    # if args.bias_frame is not None:
        
    print("\nBecause bias is provided, plotting pixels with negative median values...\n")
    locate_zeroes(nframes, neg_cosmic_map, neg_med_map, plot_output, showfig, savefig)

    # Save cosmic masks where the median values were supposed to be negative
    name_file = f'where_negative_median_is_sigma_{sigma_name}'
    save_to_xarray(neg_cosmic_map, sigma_dir, name_file, sigma_name, bias_path)

    # optionally save new fits files with cosmics replaced by median pixel values
    # note: this doesn't offer much improvement over the interpolation performed in long_slit_science_extraction.py

    # Make gifs
    # if gifs:

    #     print("\nMaking gifs...\n")

    #     # Initialize list
    #     images = list()

    #     # Gif path
    #     gif_path = os.path.join(sigma_dir, "cosmic_frames_combined.gif")
        
    #     # Read all files from folder and append to list
    #     for filename in sorted(os.listdir(gif_output)):

    #         filepath = os.path.join(gif_output, filename)

    #         if filename.endswith(".jpg") or filename.endswith(".png"):
    #             if not os.path.isfile(filepath):
    #                 continue

    #             images.append(iio.imread(filepath))
        
    #     # Make gif
    #     # Duration is in n/1000 seconds per frame. duration=500 means 0.5 seconds per frame
    #     iio.imwrite(gif_path, images, loop=1, duration=500)
    
    # print(median_values[1])
    # print('----')
    # print(cosmic_pixels[1])

    # plt.figure(figsize=[15, 10])
    # plt.imshow(cosmic_pixels[1], cmap='viridis', origin='lower')
    # plt.colorbar()
    # plt.savefig(os.path.join(direc, 'cosmic_pixels_frame_1.png'))
    
    if cleaned:
        replace_cosmics(cosmic_pixels, median_values, science_list, 
                        nints, sigma_name, instrument, cleaned_direc=cleaned_dir,
                        master_bias=master_bias)

# ------------------------------------------------------------------------------------------

if __name__ == "__main__":

    meta = xr.Dataset()

    meta.attrs['instrument'] = 'EFOSC2'
    meta.attrs['science_list'] = 'HAT-P-65_Gr11_27arcsec_list'
    meta.attrs['inputdir_Spock1'] = 'Spock0_calib_output'
    meta.attrs['outputdir_Spock1'] = 'Spock1_pre_processing'
    meta.attrs['locate_cosmics_pixelclip'] = 7
    meta.attrs['locate_cosmics_frameclip'] = 10
    meta.attrs['locate_cosmics_showfig'] = False
    meta.attrs['locate_cosmics_savefig'] = True
    meta.attrs['locate_cosmics_overwrite'] = True 

    meta.attrs['locate_cosmics_badpixelmask'] = True 
    meta.attrs['bad_pixel_mask'] = 'bad_pixel_mask_loose.h5'

    meta.attrs['locate_cosmics_groups'] = False 
    meta.attrs['locate_cosmics_gifs'] = False
    meta.attrs['locate_cosmics_cleanedfits'] = False
 
    locate_and_correct_cosmic(meta)

# ------------------------------------------------------------------------------------------

    

    # Set up the argument parser
    # -----------------------------------------

    # parser = argparse.ArgumentParser()
    # parser.add_argument('sciencelist', help="""Enter list of science .fits file names""")
    # parser.add_argument('-b','--bias_frame',help="""Define the bias frame. Not essential, can be run without a bias frame.""")
    # parser.add_argument('-mask','--bad_pixel_mask',help="""Optionally parse in a bad pixel mask to ignore these pixels from the cosmic flagging""")
    # parser.add_argument('-row','--rows',help="""Optionally define the row location of n test pixels before executing all rows""",type=int,nargs="+")
    # parser.add_argument('-col','--cols',help="""Optionally define the column location of n test pixels before executing all columns""",type=int,nargs="+")
    # parser.add_argument('-pixel_clip','--pixel_clip',help="""Define the outlier rejection threshold/sigma clip. Default = 5""",type=float,default=5.)
    # parser.add_argument('-frame_clip','--frame_clip',help="""Define the multiplicative factor at which a frame's rejection 
    #                                                         detection is deemed to have failed. Default = 3, i.e. if a frame has 3x the median number of cosmics, 
    #                                                         it's deemed to potentially have failed.""",type=float,default=3.)
    # parser.add_argument('-v','--verbose',help="""Display all cosmic pixel masks""",action='store_true')
    # parser.add_argument('-jwst','--jwst',help="""Use this option if we're looking at JWST data as the input fits files have a different format""",action='store_true')
    # parser.add_argument('-t', '--file_type', help="""Define file type for folder naming""",type=str)
    # parser.add_argument('-g', '--gif', help="""Make gifs of cosmic frame plots""",action='store_true')
    # parser.add_argument('-n', '--n_group', help="""Number of groups with different sigma clip. 
    #                                             The -pixel_clip is for the majority of the group.
    #                                             If n_group > 1 then it will prompt a question of how many
    #                                             groups that will be created.""",type=int, default=1)


    # args = parser.parse_args()

    # Optionally load in the master bias
    # if args.bias_frame is not None:
    #     master_bias = fits.open(args.bias_frame)[0].data
    #     bias = True
    # else:
    #     bias = False

    # # If using the verbose option for the median filter, I'm assuming this is a test so I don't run the full script
    # if args.rows is not None:
    #     for r,c in zip(args.rows,args.cols):
    #         locate_bad_frames(data, r, c, cut_off, verbose=True)
    #     raise SystemExit