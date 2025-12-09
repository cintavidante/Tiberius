#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk

########################################################################
# Master Flat
########################################################################

import os
import glob
import pickle
import argparse

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from astropy.io import fits
from datetime import datetime, timezone
from scipy.ndimage import median_filter, gaussian_filter
from scipy.stats import median_abs_deviation as mad
from global_utils.utils import get_bool_attr
from global_utils.instrument_utils import get_fits_image_extensions_from_config

# Prevent matplotlib plotting frames upside down
plt.rcParams['image.origin'] = 'lower'

# ---------------------------------------------------------------------------
# FLAT PLOTTING FUNCTIONS
# ---------------------------------------------------------------------------

def plot_flat_frame(f, nwin, project, instrument, verbose=False, savefig=False, flat_output=os.getcwd()):

    # Plot master flat
    plt.figure()

    if nwin == 1:
        vmin,vmax = np.nanpercentile(f,[10,90])
        plt.imshow(f,vmin=vmin,vmax=vmax,cmap='hot')
        plt.xlabel('X pixel')
        plt.ylabel('Y pixel')
    else:
        plt.subplot(121)

        vmin,vmax = np.nanpercentile(f[0],[10,90])
        plt.imshow(f[0],vmin=vmin,vmax=vmax,cmap='hot')
        plt.xlabel('X pixel')
        plt.ylabel('Y pixel')

        plt.subplot(122)
        vmin,vmax = np.nanpercentile(f[1],[10,90])
        plt.imshow(f[1],vmin=vmin,vmax=vmax,cmap='hot')

    plt.suptitle(f'Master flat before response fitting ; {project} ; {instrument}')
    plt.colorbar()
    
    if savefig:
        plt.savefig(flat_output + f'/master_flat_before_response_fitting.png')

    if verbose:
        plt.show()


# ---------------------------------------------------------------------------
# FLAT AND BAD MASK HELPERS
# ---------------------------------------------------------------------------

def bias_subtraction(master_bias, data):

    return data - master_bias

# -----------------------------------

def save_fits(data, filename, overwrite=False):

    hdu = fits.PrimaryHDU(data)
    hdu.writeto(filename, overwrite=overwrite)

    return

# -----------------------------------

def save_to_xarray(bad_pixels, bad_pixel_dir, name, method, meta=None, cut_off=None):

    # Make a xarray
    bad_pixels_da = xr.DataArray(bad_pixels, name=name)

    # Add the attributes from the meta file
    if meta is not None:
        for key, value in meta.attrs.items():
            bad_pixels_da.attrs[key] = value

    # Add new attributes
    bad_pixels_da.attrs['method'] = method
    bad_pixels_da.attrs['created'] = datetime.now(timezone.utc).isoformat()

    if cut_off != None:
        bad_pixels_da.attrs['median_cut_off'] = cut_off

    
    # Output as h5 
    bad_pixels_da.to_netcdf(bad_pixel_dir + f"/{name}.h5", engine="h5netcdf")


# ---------------------------------------------------------------------------
# FLAT COMBINATION FUNCTIONS
# ---------------------------------------------------------------------------

def combine_flats_1window(flats_list, master_bias, project, instrument, instr_cfg, idx_list, sat_limit, 
                          verbose=False, savefig=False, flat_output=os.getcwd()):

    idx = idx_list[0]
    id_slice = slice(*instr_cfg.id_slice)
    flat_data = []

    plt.figure()

    for n_frame, f in enumerate(flats_list):

        with fits.open(f) as hdul:
            data_frame = hdul[idx].data

        frame_id = f[id_slice]

        print(f'File {n_frame+1}/{len(flats_list)} ; {frame_id} ; mean={np.mean(data_frame)} ; shape={np.shape(data_frame)} ; max count={np.max(data_frame[:,20:-20])}') # Need to clip extreme edges which can have very high counts
        if np.max(data_frame[:,20:-20]) <= sat_limit:
            flat_data.append(data_frame - master_bias)
        else:
            print(f'--- ingoring potentially saturated frame {frame_id}')

        vmin, vmax = np.nanpercentile(data_frame,[10,90])
        plt.imshow(data_frame, vmin=vmin, vmax=vmax, cmap='hot')
        plt.xlabel("X pixel")
        plt.ylabel("Y pixel")
        plt.colorbar()
        plt.title(f'Frame{n_frame+1:03d} ; {frame_id} ; {project} ; {instrument}')

        if savefig:
            plt.savefig(flat_output + f'/flat_frame_{n_frame+1:03d}.png')

        if verbose:
        
            plt.show(block=False)
            plt.pause(0.5)

        plt.clf()

        hdul.close()

    flat_data = np.array(flat_data)

    mean_combine = np.mean(flat_data,axis=0)
    return mean_combine

# -----------------------------------

def combine_flats_2windows(flats_list, master_bias, project, instrument, instr_cfg, idx_list, sat_limit, 
                           verbose=False, savefig=False, flat_output=os.getcwd()):

    id_slice = slice(*instr_cfg.id_slice)
    flat_data = [[],[]]

    plt.figure()

    for n_frame, f in enumerate(flats_list):

        hdul = fits.open(f)
        frame_id = f[id_slice]

        for w, idx in enumerate(idx_list):

            data_frame = hdul[idx].data

            print(f'File {n_frame+1}/{len(flats_list)} ; {frame_id} ; window {w} ; mean={np.mean(data_frame)} ; shape={np.shape(data_frame)} ; max_count={np.max(data_frame[:,20:-20])}')
            
            if np.max(data_frame[:,20:-20]) <= sat_limit:
                flat_data[w].append(data_frame-master_bias[w])
            else:
                print('--- ingoring potentially saturated frame')
            
            plt.subplot(1,2, w+1)

            vmin,vmax = np.nanpercentile(data_frame,[10,90])
            plt.imshow(data_frame, vmin=vmin, vmax=vmax, cmap='hot')
            if w == 0:
                plt.xlabel("X pixel")
                plt.ylabel("Y pixel")
            else:
                plt.yticks(visible=False)

        plt.colorbar()
        plt.suptitle(f'Frame{n_frame+1:03d} ; {frame_id} ; {project} ; {instrument}')
        if savefig:
            plt.savefig(flat_output + f'/flat_frame_{n_frame+1:03d}.png')
        
        if verbose:
            plt.show(block=False)
            plt.pause(0.5)
            plt.clf()

        hdul.close()

    flat_data = np.array(flat_data)

    mean_combine = [np.mean(flat_data[0],axis=0),np.mean(flat_data[1],axis=0)]

    return mean_combine

# -----------------------------------

def test_smoothing_widths(flat_data, nwindows, project, instrument, showfig=False, savefig=False, flat_output=os.getcwd()):

    """Run a series of median filters with differing box widths to find optimal width."""

    if nwindows > 1:
        nrows,ncols1 = np.shape(flat_data[0])
        _,ncols2 = np.shape(flat_data[1])
        std1 = []
        std2 = []
    else:
        nrows,ncols = np.shape(flat_data)
        std = []

    for i in range(1,nrows//10,2):

        box_width = i

        if nwindows == 1:
            MF = median_filter(flat_data[:,ncols//2],box_width)
            residuals = flat_data[:,ncols//2]/MF
            # rms.append(np.sqrt(np.mean(residuals**2)))
            std.append(np.std(residuals))

        else:
            MF1 = median_filter(flat_data[0][:,ncols1//2],box_width) # just use single column for plotting
            MF2 = median_filter(flat_data[1][:,ncols2//2],box_width)

            residuals1 = flat_data[0][:,ncols1//2]/MF1
            residuals2 = flat_data[1][:,ncols2//2]/MF2

            std1.append(np.std(residuals1))
            std2.append(np.std(residuals2))

    print("Plotting bin width vs. standard deviation - MAKE NOTE OF DESIRED BIN WIDTH!")

    plt.figure()
    plt.xlabel('Bin width')
    plt.ylabel('Standard deviation')
    if nwindows == 1:
        plt.plot(np.arange(1,nrows//10,2),std,'ko')
    else:
        plt.plot(np.arange(1,nrows//10,2),std1,'bo',label='Window 1')
        plt.plot(np.arange(1,nrows//10,2),std2,'ro',label='Window 2')
        plt.legend(loc='upper left')
    
    plt.title(f'bin width vs std ; {project} ; {instrument}')
    if savefig:
        plt.savefig(flat_output + f'/bin_width_vs_std.png')

    if showfig:
        plt.show()

# -----------------------------------

def median_smooth(flat_data, nwindows, box_width, 
                  project, instrument,
                  savefits=False, showfig=False, savefig=False,
                  flat_path=os.getcwd(), flat_output=os.getcwd()):

    """Smooth the flat using a running median, and evaluated at each column individually"""

    if nwindows > 1:
        nrows1,ncols1 = np.shape(flat_data[0])
        nrows2,ncols2 = np.shape(flat_data[1])

        MF1 = median_filter(flat_data[0][:,ncols1//2],box_width) # just use single column for plotting
        MF2 = median_filter(flat_data[1][:,ncols2//2],box_width)
    else:
        nrows,ncols = np.shape(flat_data)
        MF = median_filter(flat_data[:,ncols//2],box_width)

    plt.figure(figsize=(10,8))
    plt.subplot(211)

    if nwindows == 1:
        plt.plot(flat_data[:,ncols//2],lw=4)
        plt.plot(MF,'r',lw=2)
    else:
        plt.plot(flat_data[0][:,ncols1//2],lw=4,color='b',label='window 1')
        plt.plot(MF1,'r',lw=2)
        plt.plot(flat_data[1][:,ncols2//2],lw=4,color='k',label='window 2')
        plt.plot(MF2,'r',lw=2)
        plt.legend(loc='upper left')


    plt.xticks(visible=False)
    plt.ylabel('Counts')

    plt.subplot(212)
    if nwindows == 1:
        plt.plot(flat_data[:,ncols//2]/MF)
    else:
        plt.plot(flat_data[0][:,ncols1//2]/MF1,color='b',label='window 1')
        plt.plot(flat_data[1][:,ncols2//2]/MF2,color='k',label='window 2')
        plt.legend(loc='upper left')

    plt.xlabel('Y pixel')
    plt.ylabel('Residuals')
    plt.ylim(0.951,1.049)
    plt.subplots_adjust(hspace=0)
    plt.suptitle(f'Median filter, column by column ; {project} ; {instrument}')

    if savefig:
        plt.savefig(flat_output + f'/flat_image_vs_median_filter.png')

    if showfig:
        plt.show(block=False)
        plt.pause(5)

    plt.close()

    if nwindows == 1:
        # Now find running median for each column individually
        MF_reshaped =  np.array([median_filter(flat_data[:,i],box_width) for i in range(ncols)]).transpose()

        MF_reshaped[MF_reshaped == 0] = 1 # have to replace 0s which occur near the edges, the bias over subtracts to give negative values which mess up this line

        divided_sky_flat = flat_data/MF_reshaped
        normalised_sky_flat = divided_sky_flat/divided_sky_flat.mean()

    else:
        MF1_reshaped = np.array([median_filter(flat_data[0][:,i],box_width) for i in range(ncols1)]).transpose()
        divided_sky_flat1 = flat_data[0]/MF1_reshaped
        normalised_sky_flat1 = divided_sky_flat1/divided_sky_flat1.mean()

        MF2_reshaped = np.array([median_filter(flat_data[1][:,i],box_width) for i in range(ncols2)]).transpose()
        divided_sky_flat2 = flat_data[1]/MF2_reshaped
        normalised_sky_flat2 = divided_sky_flat2/divided_sky_flat2.mean()

        normalised_sky_flat = np.array([normalised_sky_flat1,normalised_sky_flat2])


    plt.figure(figsize=(10,10))

    if nwindows == 1:
        plt.imshow(normalised_sky_flat,vmin=0.95,vmax=1.05)
    else:
        plt.subplot(121)
        plt.imshow(normalised_sky_flat[0],vmin=0.95,vmax=1.05)
        plt.subplot(122)
        plt.imshow(normalised_sky_flat[1],vmin=0.95,vmax=1.05)

    plt.suptitle(f'Median filter, column by column ; {project} ; {instrument}')
    plt.xlabel('X pixel')
    plt.ylabel('Y pixel')
    plt.colorbar()
    
    if savefig:
        plt.savefig(flat_output + f'/median_filter_column_by_column.png')

    if showfig:
        plt.show(block=False)
        plt.pause(5)

    plt.close()

    if savefits:
        hdu = fits.PrimaryHDU(normalised_sky_flat)
        hdu.writeto(os.path.join(flat_path, f'master_flat_median_filter_column_by_column_box_width_{box_width:d}_FITTED.fits'), 
                overwrite=savefits)

# -----------------------------------

def gaussian_smooth(flat_data, nwindows, project, instrument, 
                    savefits=False, showfig=False, savefig=False, 
                    flat_path=os.getcwd(), flat_output=os.getcwd()):

    """Smooth the flat using a Gaussian filter"""

    if instrument == 'ACAM':
        sigma = 1
    if instrument == 'EFOSC2':
        sigma = 0.4 # this is found through inspection by eye

    if nwindows > 1:
        nrows1,ncols1 = np.shape(flat_data[0])
        nrows2,ncols2 = np.shape(flat_data[1])

        GF1 = gaussian_filter(flat_data[0],sigma)
        GF2 = gaussian_filter(flat_data[1],sigma)
    else:
        nrows,ncols = np.shape(flat_data)
        GF = gaussian_filter(flat_data,sigma)


    # Plot Gaussian filter evaluated at a column in centre of CCD
    plt.figure(figsize=(10,8))
    plt.subplot(211)

    if nwindows == 1:
        plt.plot(flat_data[:,ncols//2],lw=6)
        plt.plot(GF[:,ncols//2],'r',lw=2)
    else:
        plt.plot(flat_data[0][:,ncols1//2],lw=4,color='b',label='window 1')
        plt.plot(GF1[:,ncols1//2],'r',lw=2)
        plt.plot(flat_data[1][:,ncols2//2],lw=4,color='k',label='window 2')
        plt.plot(GF2[:,ncols2//2],'r',lw=2)
        plt.legend(loc='upper left')

    plt.xticks(visible=False)
    plt.ylabel('Counts')

    plt.subplot(212)
    if nwindows == 1:
        plt.plot(flat_data[:,ncols//2]/GF[:,ncols//2])
    else:
        plt.plot(flat_data[0][:,ncols1//2]/GF1[:,ncols1//2],color='b',label='window 1')
        plt.plot(flat_data[1][:,ncols2//2]/GF2[:,ncols2//2],color='k',label='window 2')
        plt.legend(loc='upper left')

    plt.xlabel('Y pixel')
    plt.ylabel('Residuals')
    plt.ylim(0.951,1.049)
    plt.subplots_adjust(hspace=0)
    plt.suptitle(f'Gaussian filter, example column ; {project} ; {instrument}')

    if savefig:
        plt.savefig(flat_output + f'/gaussian_filter_example_column.png')

    if showfig:
        plt.show(block=False)
        plt.pause(5)

    plt.close()


    if nwindows == 1:
        divided_sky_flat = flat_data/GF
        normalised_sky_flat = divided_sky_flat/divided_sky_flat.mean()

    else:
        divided_sky_flat1 = flat_data[0]/GF1
        normalised_sky_flat1 = divided_sky_flat1/divided_sky_flat1.mean()

        divided_sky_flat2 = flat_data[1]/GF2
        normalised_sky_flat2 = divided_sky_flat2/divided_sky_flat2.mean()

        normalised_sky_flat = np.array([normalised_sky_flat1,normalised_sky_flat2])


    plt.figure(figsize=(10,10))
    if nwindows == 1:
        plt.imshow(normalised_sky_flat,vmin=0.95,vmax=1.05)
    else:
        plt.subplot(121)
        plt.imshow(normalised_sky_flat[0],vmin=0.95,vmax=1.05)
        plt.subplot(122)
        plt.imshow(normalised_sky_flat[1],vmin=0.95,vmax=1.05)

    plt.suptitle(f'Gaussian filter ; {project} ; {instrument}')
    plt.xlabel('X pixel')
    plt.ylabel('Y pixel')
    plt.colorbar()
    plt.show(block=False)

    plt.pause(5)
    plt.close()

    if savefits:
        hdu = fits.PrimaryHDU(normalised_sky_flat)
        hdu.writeto(os.path.join(flat_path, f'master_flat_Gaussian_filter_FITTED.fits'), overwrite=savefits)

# ---------------------------------------------------------------------------
# BAD PIXEL MASKING FUNCTIONS
# ---------------------------------------------------------------------------

def pixel_mask_tight(frame, project, instrument, 
                     mask_showfig=False, mask_savefig=False,
                     bad_output=os.getcwd()):
    """A tighter definition of a bad pixel mask using 5 median absolute devitations from the median, 
    computed column by column (along dispersion direction)"""

    medians = np.median(frame,axis=0)
    mads = mad(frame,axis=0,scale='normal')
    good_pixels = []

    for i,row in enumerate(frame):
        good_pixels.append(((row >= medians-5*mads) & (row <= medians+5*mads)))

    good_pixels = np.array(good_pixels)
    bad_pixels = ~good_pixels

    percentage_bad_pixels = 100*len(np.where(bad_pixels)[0])/(frame.shape[0]*frame.shape[1])
    print("Percentage of pixels deemed bad by medians and mads = %.2f%%"%percentage_bad_pixels)

    if mask_savefig or mask_showfig:
        plt.figure()
        plt.imshow(bad_pixels)
        plt.title(f"Bad pixel mask using medians and MADS ; {project} ; {instrument}")

        if mask_savefig:
            plt.savefig(bad_output + '/bad_pixel_mask_medians_MADS.png')

        if mask_showfig:
            plt.show()

    return bad_pixels

# -----------------------------------

def pixel_mask_loose(frame, project, instrument,
                     mask_showfig=False, mask_savefig=False,
                     bad_output=os.getcwd()):
    """A looser definition of the bad pixel mask using 5 standard deviations from the global mean"""

    good_pixels = ((frame >= np.mean(frame)-5*np.std(frame)) & (frame <= np.mean(frame)+5*np.std(frame)))
    bad_pixels = ~good_pixels

    percentage_bad_pixels = 100*len(np.where(bad_pixels)[0])/(frame.shape[0]*frame.shape[1])
    print("Percentage of pixels deemed bad by means and stds = %.2f%%"%percentage_bad_pixels)

    if mask_savefig or mask_showfig:

        plt.figure()
        plt.imshow(bad_pixels)
        plt.title(f"Bad pixel mask using means and stds ; {project} ; {instrument}")

        if mask_savefig:
            plt.savefig(bad_output + '/bad_pixel_mask_means_stds.png')

        if mask_showfig:
            plt.show()

    return bad_pixels

# -----------------------------------

def pixel_mask_medfilt(frame, project, instrument, 
                       mask_showfig=False, mask_savefig=False, mask_cutoff=5,
                       bad_output=os.getcwd()):
    """A function that uses median filters along the rows/cross-dispersion direction to locate outliers. 
    Can be more effective."""

    good_pixels = []

    for i,row in enumerate(frame):
        MF = median_filter(row, mask_cutoff)
        residuals = row-MF
        good_pixels.append(((residuals >= -10*mad(residuals,scale='normal')) & (residuals <= 10*mad(residuals,scale='normal'))))

    good_pixels = np.array(good_pixels)
    bad_pixels = ~good_pixels

    percentage_bad_pixels = 100*len(np.where(bad_pixels)[0])/(frame.shape[0]*frame.shape[1])
    print("Percentage of pixels deemed bad by median filter = %.2f%%"%percentage_bad_pixels)

    if mask_savefig or mask_showfig:

        plt.figure()
        plt.imshow(bad_pixels)
        plt.title(f"Bad pixel mask using median filter ; {project} ; {instrument}")
        plt.xlabel('X pixel')
        plt.ylabel('Y pixel')

        if mask_savefig:
            plt.savefig(bad_output + '/bad_pixel_mask_median_filter.png')

        if mask_showfig:
            plt.show()

    return bad_pixels

# ---------------------------------------------------------------------------
# MAIN CALLABLE FUNCTIONS
# ---------------------------------------------------------------------------

def create_pixel_mask(f, meta, nwin, project, instrument,  mask_savefits, mask_showfig, mask_savefig, mask_cutoff, bad_pixel_dir, bad_output):

    print('\n')
    print('Testing different pixel masks.')
    # Make bad pixel masks
    if nwin == 1:
        pmt = pixel_mask_tight(f, project, instrument, mask_showfig, mask_savefig, bad_output)
        pml = pixel_mask_loose(f, project, instrument, mask_showfig, mask_savefig, bad_output)
        pmmf = pixel_mask_medfilt(f, project, instrument, mask_showfig, mask_savefig, mask_cutoff, bad_output)
    else:
        pmt = np.array([pixel_mask_tight(f[0], project, instrument, mask_showfig, mask_savefig),
                        pixel_mask_tight(f[1], project, instrument, mask_showfig, mask_savefig)])

        pml = np.array([pixel_mask_loose(f[0], project, instrument, mask_showfig, mask_savefig),
                        pixel_mask_loose(f[1], project, instrument, mask_showfig, mask_savefig)])

        pmmf = np.array([pixel_mask_medfilt(f[0], project, instrument, mask_showfig, mask_savefig, mask_cutoff),
                         pixel_mask_medfilt(f[1], project, instrument, mask_showfig, mask_savefig, mask_cutoff)])

    if mask_savefits:
        #pickle.dump(pmt,open("bad_pixel_mask_tight.pickle","wb"))
        save_to_xarray(pmt,
                    bad_pixel_dir,
                    name='bad_pixel_mask_tight',
                    method='tighter definition using 5 std from median',
                    meta=meta)
            
        #pickle.dump(pml,open("bad_pixel_mask_loose.pickle","wb"))
        save_to_xarray(pml,
                    bad_pixel_dir,
                    name='bad_pixel_mask_loose',
                    method='looser def, 5 std from the global mean',
                    meta=meta)
            
        #pickle.dump(pmmf,open("bad_pixel_mask_medfilt.pickle","wb"))
        save_to_xarray(pmmf,
                    bad_pixel_dir,
                    name=f'bad_pixel_mask_medfilt_{mask_cutoff}',
                    method='median filter',
                    meta=meta,
                    cut_off=mask_cutoff)

# -----------------------------------

def create_master_flat_pixel_mask(meta, instr_cfg):

    # Read meta database

    instrument = meta.attrs["instrument"]
    project = meta.attrs["project_name"]
    savedir = meta.attrs["outputdir_Spock1"]
    
    flatlist = meta.attrs["flat_list"]
    savefits = get_bool_attr(meta.attrs["flat_savefits"])
    showfig = get_bool_attr(meta.attrs["flat_showfig"])
    savefig = get_bool_attr(meta.attrs["flat_savefig"])
    satlimit = float(meta.attrs["flat_saturation_limit"])
    boxwidth = int(meta.attrs["flat_box_width"])

    mask = get_bool_attr(meta.attrs["bad_pixel_mask"])
    mask_savefits = get_bool_attr(meta.attrs["bad_pixel_mask_savefits"])
    mask_savefig = get_bool_attr(meta.attrs["bad_pixel_mask_savefig"])
    mask_showfig = get_bool_attr(meta.attrs["bad_pixel_mask_showfig"])
    mask_medianfilter_cutoff = int(meta.attrs["bad_pixel_mask_medianfilter_cutoff"])


    # Define directory for master flat and bad pixel mask
    flat_dir = os.path.join(savedir, 'flat', os.path.basename(flatlist))
    flat_output = os.path.join(flat_dir, 'plots')
    bad_pixel_dir = os.path.join(savedir, 'bad_pixel_mask', os.path.basename(flatlist))
    bad_output = os.path.join(bad_pixel_dir, 'plots')

    # Make the directory
    os.makedirs(flat_dir, exist_ok=True)
    os.makedirs(flat_output, exist_ok=True)
    if mask:
        os.makedirs(bad_pixel_dir, exist_ok=True)
        os.makedirs(bad_output, exist_ok=True)

    # Check instrument 
    if instrument != 'EFOSC2' and instrument != 'ACAM':
        raise ValueError('Currently only set up / tested to deal with ACAM or EFOSC data')
    
    # Load flat list
    flats_files = np.loadtxt(flatlist, str)

    # Find master bias
    bias_path = glob.glob(os.path.join(savedir, 'bias', 'master_bias.fits'))[0]

    # If bias is not there, raise error
    if not os.path.isfile(bias_path):
        raise FileNotFoundError(f"Bias file not found at: {bias_path}")

    # Load master bias
    master_bias_data = fits.open(bias_path)[0].data

    # Find out how many windows we're dealing with
    test_file = flats_files[0]
    idx_list = get_fits_image_extensions_from_config(instr_cfg, test_file)
    nwin = len(idx_list)
    
    # Create master flat for different windows
    if nwin == 1:
        f = combine_flats_1window(flats_files, master_bias_data, project, instrument, instr_cfg, idx_list, 
                                  satlimit, showfig, savefig, flat_output)
    else:
        f = combine_flats_2windows(flats_files, master_bias_data, project, instrument, instr_cfg, idx_list, 
                                   satlimit, showfig, savefig, flat_output)


    # Plot Master Flat before Response Fitting
    if showfig or savefig:
        plot_flat_frame(f, nwin, project, instrument, showfig, savefig, flat_output)

    # Save flats
    if savefits:
        save_fits(f/np.median(f), os.path.join(flat_dir, f'master_flat_{project}_{instrument}.fits'), overwrite=savefits)

    # Smoothing width
    test_smoothing_widths(np.array(f), nwin, project, instrument, showfig, savefig, flat_output)

    # Median smooth
    median_smooth(np.array(f), nwin, boxwidth, project, instrument, savefits, 
                  showfig, savefig, flat_dir, flat_output)
    
    # Gaussing smooth
    # ## Note: Gaussian smooth not currently used, despite good perfomance, as it removes features in x as well as y
    gaussian_smooth(np.array(f), nwin, project, instrument, savefits,
                    showfig, savefig, flat_dir, flat_output)

    # Create bad pixel mask
    if mask:
        create_pixel_mask(f, meta, nwin, project, instrument,
                           mask_savefits, mask_showfig, mask_savefig,
                             mask_medianfilter_cutoff, bad_pixel_dir, bad_output)

    return

# ---------------------------------------------------------------------------
# Optional: Testing the script
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    meta = xr.Dataset()
    meta.attrs['instrument'] = 'EFOSC2'
    meta.attrs['projet_name'] = 'Planet'
    meta.attrs['inputdir_Spock1'] = 'Spock0_calib_output'
    meta.attrs['outputdir_Spock1'] = 'Spock1_pre_processing'
    
    meta.attrs['flat_list'] = '/Users/cividante/Tiberius/hatp65/Spock0_calib_output/SKY_FLAT_Gr11_27arcsec_list'
    meta.attrs['flat_savefits'] = True
    meta.attrs['flat_showfig'] = True
    meta.attrs['flat_savefig'] = True
    meta.attrs['flat_saturation_limit'] = 55000
    meta.attrs['flat_box_width'] = 9

    meta.attrs['bad_pixel_mask'] = True
    meta.attrs['bad_pixel_mask_savefits'] = True
    meta.attrs['bad_pixel_mask_showfig'] = True
    meta.attrs['bad_pixel_mask_savefig'] = True
    meta.attrs['bad_pixel_mask_medianfilter_cutoff'] = 5


    create_master_flat_pixel_mask(meta, None)
