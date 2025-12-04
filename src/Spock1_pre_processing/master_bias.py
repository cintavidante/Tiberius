#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk

#############################################################################
# Master Bias
#############################################################################

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import sys
import glob
from PIL import Image
from pathlib import Path
from global_utils.utils import get_bool_attr
from global_utils.instrument_utils import get_fits_image_extensions_from_config

# Prevent matplotlib flipping images
plt.rcParams['image.origin'] = 'lower'

# ---------------------------------------------------------------------------
#  PLOTTING HELPERS
# ---------------------------------------------------------------------------

def plot_bias_frame(data, title=None, save_as=None, verbose=False, savefig=False):
    """Plot the single bias frame or the multiple bias frames for one science file"""
    plt.figure(figsize=(10, 6))

    if isinstance(data, list):
        # Multiple windows (for example ACAM)
        for i, img in enumerate(data):
            plt.subplot(1, len(data), i+1)
            plt.imshow(
                img,
                vmin=np.median(img) * 0.99,
                vmax=np.median(img) * 1.01,
                cmap="viridis"
                )
            plt.title(f"Window {i+1}")
            plt.colorbar(label="ADU")
            plt.xlabel("X pixel")
            plt.ylabel("Y pixel")
        
        if title:
            plt.suptitle(title)
    
    else: # Single window
        plt.imshow(
            data,
            vmin=np.median(data) * 0.99,
            vmax=np.median(data) * 1.01,
            cmap="viridis"
            )
        plt.colorbar(label="ADU")
        plt.xlabel("X pixel")
        plt.ylabel("Y pixel")
        
        if title:
            plt.title(title)

    if savefig:
        plt.savefig(save_as, dpi=120)
    if verbose:
        plt.show()

    plt.close()


def plot_master_bias(master_bias, instrument, save_dir):
    """Plot the resulting master bias"""
    plt.figure(figsize=(11, 7))
    plt.suptitle(f"Master Bias – {instrument}", size=16)

    if master_bias.ndim == 2:
        plt.imshow(
            master_bias,
            vmin=np.median(master_bias)*0.99,
            vmax=np.median(master_bias)*1.01,
            cmap="viridis"
        )
        plt.colorbar()
        plt.xlabel("X pixel")
        plt.ylabel("Y pixel")
    else:
        for i in range(len(master_bias)):
            plt.subplot(1, 2, i + 1)
            plt.title(f"Window {i+1}")
            plt.imshow(master_bias[i], cmap="gray")
            plt.colorbar()
            plt.xlabel("X pixel")
            plt.ylabel("Y pixel")

    Path(save_dir).mkdir(exist_ok=True)
    plt.savefig(f"{save_dir}/master_bias_{instrument}.png", dpi=120)
    plt.savefig(f"{save_dir}/master_bias_{instrument}.pdf")
    plt.close()


def create_bias_gif(instrument, save_dir):
    """Create animated GIF from saved bias PNG frames"""
    plot_dir = save_dir + '/plots'
    pattern = f"{plot_dir}/{instrument}_*.png"
    png_files = sorted(glob.glob(pattern))

    if not png_files:
        print("No PNG files found to create GIF.")
        return

    frames = [Image.open(f) for f in png_files]
    gif_name = f"{save_dir}/master_bias_{instrument}.gif"

    frames[0].save(
        gif_name,
        save_all=True,
        append_images=frames[1:],
        duration=400,
        loop=0
    )

    print(f"GIF created → {gif_name}")


# ---------------------------------------------------------------------------
#  BIAS COMBINATION FUNCTIONS
# ---------------------------------------------------------------------------

def combine_biases_1window(bias_files, instrument, instr_cfg, idx_list, save_dir, verbose=False, savefig=False, eyeball=False):
    """Calculate the combined biases for 1 window"""
    idx = idx_list[0]
    id_slice = slice(*instr_cfg.id_slice)

    bias_data = []
    good_frames, bad_frames = [], []

    for n, f in enumerate(bias_files):
        with fits.open(f) as hdul:
            data_frame = hdul[idx].data

        frame_id = f[id_slice]
        mean, var = np.mean(data_frame), np.var(data_frame)
        print(f"File {n+1}/{len(bias_files)} ; {frame_id} ; mean={mean:.1f}, var={var:.1f}, var/mean={var/mean:.2f}")

        if verbose or savefig:
            plot_dir = save_dir + '/plots'
            Path(plot_dir).mkdir(exist_ok=True)
            png_name = f"{plot_dir}/{instrument}_{n+1:03d}.png"
            plot_bias_frame(data_frame, title=f"{frame_id}", save_as=png_name, verbose=verbose, savefig=savefig)

        if eyeball:
            ans = input("good/bad? [g/b]: ").lower()
            (good_frames if ans == "g" else bad_frames).append(f)

        bias_data.append(data_frame)

    if eyeball:   # plot each frame and have user select what is good and bad.  NEEDS IMPLEMENTING FOR 2 WINDOWS.
        good_path = Path(save_dir) / "bias_good.txt"
        bad_path  = Path(save_dir) / "bias_bad.txt"

        np.savetxt(good_path, good_frames, fmt="%s")
        np.savetxt(bad_path,  bad_frames,  fmt="%s")
        sys.exit()

    return np.median(np.array(bias_data), axis=0)


def combine_biases_2windows(bias_files, instrument, instr_cfg, idx_list, save_dir, verbose=False, savefig=False): 
    """Calculate the combined biases for 2 windows seperately"""
    id_slice = slice(*instr_cfg.id_slice)


    bias_data = [[] for _ in idx_list]

    for n, f in enumerate(bias_files):
        with fits.open(f) as hdul:
            frame_id = f[id_slice]
            
            window_frames = []
            for i, idx in enumerate(idx_list):
                data_frame = hdul[idx].data
                bias_data[i].append(data_frame)
                window_frames.append(data_frame)

                mean, var = np.mean(data_frame), np.var(data_frame)
                print(f"File {n+1}/{len(bias_files)} ; window {idx} ; {frame_id} ; mean={mean:.1f}, var={var:.1f}, var/mean={var/mean:.2f}")

            if verbose or savefig:
                plot_dir = save_dir + '/plots'
                Path(plot_dir).mkdir(exist_ok=True)
                png_name = f"{plot_dir}/{instrument}_{n+1:03d}.png"
                plot_bias_frame(window_frames, title=f"{frame_id}", save_as=png_name, verbose=verbose, savefig=savefig)

    return np.array([np.median(bias_data[i], axis=0) for i in range(len(idx_list))])


# ---------------------------------------------------------------------------
#  MAIN CALLABLE FUNCTION
# ---------------------------------------------------------------------------

def create_master_bias(meta, instr_cfg):
    """Main function to create a master bias"""
    attr = getattr(meta, "ds", meta).attrs

    biaslist = attr["bias_list"]
    instrument = attr["instrument"]
    save_dir = attr["outputdir_Spock1"]
    save_dir = './' + save_dir + '/bias'
    verbose = get_bool_attr(attr["bias_verbose"])
    savefig = get_bool_attr(attr["bias_savefig"])
    savefits = get_bool_attr(attr["bias_clobber"])
    eyeball = get_bool_attr(attr["bias_eyeball"])

    Path(save_dir).mkdir(exist_ok=True)
    
    bias_files = np.loadtxt(biaslist, str)
    
    test_file = bias_files[0]
    idx_list = get_fits_image_extensions_from_config(instr_cfg, test_file)
    nwin = len(idx_list)

    # Combine frames
    if nwin == 1:
        master_bias = combine_biases_1window(
            bias_files, instrument, instr_cfg, idx_list, save_dir, verbose, savefig, eyeball
        )
    else:
        master_bias = combine_biases_2windows(
            bias_files, instrument, instr_cfg, idx_list, save_dir, verbose, savefig
        )

    # Save FITS
    if savefits:
        fits.PrimaryHDU(master_bias).writeto(f"{save_dir}/master_bias.fits", overwrite=True)

    # Create GIF if requested
    if verbose:
        create_bias_gif(instrument, save_dir)

    # Plot master bias
    plot_master_bias(master_bias, instrument, save_dir)

    print("\nMaster bias generation complete.\n")

    return master_bias


# ---------------------------------------------------------------------------
# Optional: Testing the script
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import xarray as xr 
    from global_utils.InstrumentDataClass import load_instrument_config

    meta = xr.Dataset()
    meta.attrs['instrument'] = 'EFOSC2'
    meta.attrs['inputdir_Spock1'] = 'Spock0_calib_output'
    meta.attrs['outputdir_Spock1'] = 'Spock1_pre_processing'

    meta.attrs['bias_list'] = './Spock0_calib_output/BIAS_Free_27arcsec_list'
    meta.attrs['bias_verbose'] = True
    meta.attrs['bias_clobber'] = True
    meta.attrs['bias_savefig'] = True
    meta.attrs['bias_eyeball'] = False

    instr_cfg = load_instrument_config(meta.attrs["instrument"].strip("'"))
    create_master_bias(meta, instr_cfg)
