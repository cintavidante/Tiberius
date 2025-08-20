#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
from scipy.ndimage import median_filter
import pickle
import os
import copy
import glob
import imageio.v3 as iio
import shutil

# Prevent matplotlib plotting frames upside down
plt.rcParams['image.origin'] = 'lower'

# ------------------------------------------------------------------------------------------

def locate_bad_frames(image_data, pixel_row, pixel_col, cut_off_arr, verbose=False):
    """The function that locates the frames/pixels where cosmics are located.

    Inputs:
    image_data - the image data, dimensions of nframes x nrows x ncols
    pixel_row - the row number of the pixel under consideration
    pixel_col - the column number of the pixel under consideration
    cut_off - the sigma cut off / outlier rejection threshold
    verbose - True/False - plot the outlier identififcation for this pixel?

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

    if verbose: # plot output
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
        plt.axhline(cut_off*np.nanstd(residuals),ls='--',color='k')
        plt.axhline(cut_off*-np.nanstd(residuals),ls='--',color='k',label="cut-off")
        plt.plot(np.arange(nframes)[bad_frames],residuals[bad_frames],"rx",label="Flagged outliers")
        plt.legend()
        plt.show()
        plt.savefig('locate_cosmic.png')

    return bad_frames, median

# ---------------------------------------------

def plot_cosmic_frames(cosmic_pixels, cut_off_arr, plotdir):

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

        plt.savefig(filename)
        plt.show(block=False)
        plt.pause(1e-6)
        plt.clf()
    
    return 

# ---------------------------------------------

def check_cosmic_frames(cosmic_pixels, frame_cut_off):
    """A function that plots and optionally resets cosmic pixels for frames where a disproportionate number of pixels have been flagged as cosmics.

    Inputs:
    cosmic_pixels - the array of all cosmic flagged pixels, dimensions of nframes x nrows x ncols

    Returns:
    cosmic_pixels - the new array of all cosmic flagged pixels, taking into the account the user-defined reset frame masks"""

    nframes,nrows,ncols = cosmic_pixels.shape

    ncosmics = []

    for i,c in enumerate(cosmic_pixels):
        ncosmics.append(len(np.where(c)[0]))

    median_cosmics = np.nanmedian(ncosmics)
    if median_cosmics == 0:
        median_cosmics = 1

    print("Median number of cosmics per frame = %d (%.3f%%)"%(median_cosmics,100*median_cosmics/(nrows*ncols)))

    incorrectly_flagged_cosmics = []

    for i,c in enumerate(cosmic_pixels):
        ncosmics = len(np.where(c==1)[0])
        if ncosmics > frame_cut_off*median_cosmics:
            print("Integration %d has %.2fX the median number of cosmics, somethings up"%(i,ncosmics/median_cosmics))
            plt.figure()
            plt.imshow(c,cmap='Greys', interpolation='none',aspect="auto")
            incorrectly_flagged_cosmics.append(i)
            plt.title("Integration %d"%i)
            plt.ylabel("Pixel row")
            plt.xlabel("Pixel column")
            plt.show(block=False)
            # plt.clf()

            reset_mask = input("Reset mask for integration %d? [y/n]: "%i)
            if reset_mask == "y":
                print("...resetting mask\n")
                cosmic_pixels[i] = np.zeros_like(c)

    return cosmic_pixels

# ---------------------------------------------

def replace_cosmics(cosmic_pixels, medians, science_list, nints, cut_off_name, jwst=False, cleaned_direc=os.getcwd()):

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

    if jwst:

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

        # fits_files = [fits.open(s,memmap=False) for s in science_list]
        # new_fits_files = [copy.deepcopy(f) for f in fits_files]
        # filenames = [s.split("/")[-1] for s in science_list]
        # nints = np.cumsum([f["SCI"].data.shape[0] for f in fits_files])
        # total_nints = nints[-1]

        # for s in science_list:
        #     fits_file = fits.open(s,memmap=False)
        #     new_fits_file = copy.deepcopy(fits_file)
        #     filename = s.split("/")[-1]

        #     for i,c in enumerate(cosmic_pixels):

        #         jwst_fits_counter = np.digitize(i,nints)

        #         if jwst_fits_counter > 0:
        #             jwst_index_counter = i-nints[jwst_fits_counter]
        #         else:
        #             jwst_index_counter = i

        #         print("Cleaning integration %d, %s"%(i,filenames[jwst_fits_counter]))

        #         for row in range(nrows):
        #             new_fits_files[jwst_fits_counter]["SCI"].data[jwst_index_counter][row][c[row]] = medians[i][row][c[row]]

        #     for i,nf in enumerate(new_fits_files):
        #         print("Saving cosmic_cleaned_fits/%s"%(filenames[i]))
        #         nf.writeto("cosmic_cleaned_fits/%s"%filenames[i],overwrite=True)

        # return
    
    medians = medians.astype('uint16')  
    # array_loop = [1, 3, 5, 7, 9]

    for i in range(nframes):

        f = fits.open(science_list[i])
        f_new = copy.deepcopy(f)
        filename = science_list[i].split("/")[-1]

        print("Cleaning frame {}".format(i+1))

        for row in range(nrows):
            f_new[0].data[row][cosmic_pixels[i][row]] = medians[i][row][cosmic_pixels[i][row]]

        file_path = os.path.join(cleaned_direc, filename)

        f_new.writeto(file_path, overwrite=True)
        f.close()
    
    # Make a new file list for the cleaned cosmic fits
    print('Making new list for cleaned cosmic fits...')

    # Read the global path for all cleaned cosmic fits
    file_path = os.path.join(cleaned_direc, '*.fits')
    all_files = sorted(glob.glob(file_path))

    # Write a new list. This makes sure it's overwriting new files!!!
    with open(args.sciencelist + '_cleaned_' + cut_off_name, 'w') as new_list:

        for files in all_files:
            new_list.write(files+' \n')

    return

# ------------------------------------------------------------------------------------------

if __name__ == "__main__":

    # Set up the argument parser
    # -----------------------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument('sciencelist', help="""Enter list of science .fits file names""")
    parser.add_argument('-b','--bias_frame',help="""Define the bias frame. Not essential, can be run without a bias frame.""")
    parser.add_argument('-mask','--bad_pixel_mask',help="""Optionally parse in a bad pixel mask to ignore these pixels from the cosmic flagging""")
    parser.add_argument('-row','--rows',help="""Optionally define the row location of n test pixels before executing all rows""",type=int,nargs="+")
    parser.add_argument('-col','--cols',help="""Optionally define the column location of n test pixels before executing all columns""",type=int,nargs="+")
    parser.add_argument('-pixel_clip','--pixel_clip',help="""Define the outlier rejection threshold/sigma clip. Default = 5""",type=float,default=5.)
    parser.add_argument('-frame_clip','--frame_clip',help="""Define the multiplicative factor at which a frame's rejection 
                                                            detection is deemed to have failed. Default = 3, i.e. if a frame has 3x the median number of cosmics, 
                                                            it's deemed to potentially have failed.""",type=float,default=3.)
    parser.add_argument('-v','--verbose',help="""Display all cosmic pixel masks""",action='store_true')
    parser.add_argument('-jwst','--jwst',help="""Use this option if we're looking at JWST data as the input fits files have a different format""",action='store_true')
    parser.add_argument('-t', '--file_type', help="""Define file type for folder naming""",type=str)
    parser.add_argument('-g', '--gif', help="""Make gifs of cosmic frame plots""",action='store_true')
    parser.add_argument('-n', '--n_group', help="""Number of groups with different sigma clip. 
                                                The -pixel_clip is for the majority of the group.
                                                If n_group > 1 then it will prompt a question of how many
                                                groups that will be created.""",type=int, default=1)


    args = parser.parse_args()

    # Loading data
    # ----------------------------------------

    # Load in the list of science file names
    science_list = np.atleast_1d(np.loadtxt(args.sciencelist,dtype=str))

    # Optionally load in the master bias
    if args.bias_frame is not None:
        master_bias = fits.open(args.bias_frame)[0].data
        bias = True
    else:
        bias = False

    # Load in the science data
    print("Loading in data...")
    data = []
    nints = []

    for s in science_list:
        f = fits.open(s,memmap=False)
        if args.jwst:
            data.append(f["SCI"].data)
            nints.append(f["SCI"].data.shape[0])
        else:
            if bias: # subtract the bias if using one
                data.append(f[0].data-master_bias)
            else:
                data.append(f[0].data)
            nints.append(f[0].data.shape[0])
        f.close()

    nints = np.cumsum(nints)

    # Define data depending if it's JWST data or not 
    if args.jwst:
        data = np.vstack(data)
    else:
        data = np.array(data)

    # Load bad pixel mask if provided
    if args.bad_pixel_mask is not None:
        mask = pickle.load(open(args.bad_pixel_mask,"rb"))
        data[:,mask] = np.nan

    # Definitons
    # ----------------------------------------

    # Define the cosmic pixel flagged array, initially as an array of zeros matching the dimensions of the input data
    cosmic_pixels = np.zeros_like(data)
    nframes, nrows, ncols = data.shape

    # Define the sigma cut off
    cut_off = args.pixel_clip

    # Make an array of pixel_clip
    cut_off_arr = np.full(nframes, cut_off, dtype=float)

    # print(cut_off_arr)

    # If n_group > 1, then we will ask which frames and which sigma values to use
    if args.n_group > 1:

        for m in range(args.n_group -1):

            # Ask for sigma clip
            cut_val = float(input("Enter sigma value: "))
            frame_list = input("Enter frames number separated by space: ")
            frame_list = list(map(int, frame_list.split()))
            frame_list = np.array(frame_list) - 1  # Convert to zero-based index

            cut_off_arr[frame_list] = cut_val
    
    # print(cut_off_arr)

    # # Change the first values to be like 8
    # cut_off_arr[:3] = 8
    # cut_off_arr[3:6] = 7
    # cut_off_arr[:-3] = 7

    # Find unique values in sigma / cut_off
    sigmas = np.unique(cut_off_arr)

    # Join them the unique sigmas to be a string
    cut_off_name = "_".join([str(s) for s in sigmas])

    # Set up working directory
    # -----------------------------------------

    # Get the current working directory
    pwd = os.getcwd()

    # Automatically make a new folder
    # The pixel clip here is the majority of the clip for the frames
    direc = os.path.join(pwd, 'cosmic/{}_sigma_{}'.format(args.file_type, cut_off_name))

    # Make sure that directory exists
    os.makedirs(direc, exist_ok=True)

    # # If using the verbose option for the median filter, I'm assuming this is a test so I don't run the full script
    # if args.rows is not None:
    #     for r,c in zip(args.rows,args.cols):
    #         locate_bad_frames(data, r, c, cut_off, verbose=True)
    #     raise SystemExit

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
            bad_frames, medians = locate_bad_frames(data, row, col, cut_off_arr, verbose=False)

            # Append to the cosmic pixels and median values arrays
            cosmic_pixels[:,row,col][bad_frames] = 1
            median_values[:,row,col] = medians

    # Change the cosmic pixels array to boolean type
    cosmic_pixels = cosmic_pixels.astype(bool)

    # Outputs
    # ----------------------------------------

    # Plotting directory
    plot_dir = os.path.join(direc, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Cleaned fits directory
    cleaned_dir = os.path.join(direc, "cleaned_fits")
    os.makedirs(cleaned_dir, exist_ok=True)

    if args.verbose:

        print("\nPlotting all cosmic-masked pixels...\n")
        plot_cosmic_frames(cosmic_pixels, cut_off_arr, plot_dir)

    # double-check the output
    if not args.verbose:

        print("Plotting frames with high number of cosmics...\n")
        cosmic_pixels = check_cosmic_frames(cosmic_pixels,args.frame_clip)

    # Save the cosmic masks
    pickle.dump(cosmic_pixels, 
                open(os.path.join(direc, "cosmic_pixel_mask_{}_sigma_clip.pickle".format(cut_off_name)), "wb"))

    # optionally save new fits files with cosmics replaced by median pixel values
    # note: this doesn't offer much improvement over the interpolation performed in long_slit_science_extraction.py

    # Make gifs
    if args.gif:

        print("\nMaking gifs...\n")

        # Initialize list
        images = list()

        # Gif path
        gif_path = os.path.join(plot_dir, "cosmic_frames_combined.gif")
        
        # Read all files from folder and append to list
        for filename in sorted(os.listdir(plot_dir)):

            filepath = os.path.join(plot_dir, filename)

            if filename.endswith(".jpg") or filename.endswith(".png"):
                if not os.path.isfile(filepath):
                    continue

                images.append(iio.imread(filepath))
        
        # Make gif
        # Duration is in n/1000 seconds per frame. duration=500 means 0.5 seconds per frame
        iio.imwrite(gif_path, images, loop=1, duration=500)
    
    # print(median_values[1])
    # print('----')
    # print(cosmic_pixels[1])

    # plt.figure(figsize=[15, 10])
    # plt.imshow(cosmic_pixels[1], cmap='viridis', origin='lower')
    # plt.colorbar()
    # plt.savefig(os.path.join(direc, 'cosmic_pixels_frame_1.png'))

    # Replace cosmic?
    replace = input("Replace cosmic values with median and save to new fits? [y/n]: ")
    
    if replace == "y":
        replace_cosmics(cosmic_pixels, median_values, science_list, 
                        nints, cut_off_name, args.jwst, cleaned_direc=cleaned_dir)
    
    # print('Max value of master bias: {}'.format(np.max(master_bias)))
    # print('Min value of master bias: {}'.format(np.min(master_bias)))
