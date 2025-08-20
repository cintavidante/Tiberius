import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
import os
import glob

# Run this python script in calib_files directory

# ------------------------------------------------------------------------------------------

if __name__ == "__main__":

    # Define parser
    parser = argparse.ArgumentParser()

    parser.add_argument('sciencelist', help="""Enter list of science .fits file names""")
    parser.add_argument('-n','--n_sep',help="""Number of separated list.""", type=int, default=2)

    args = parser.parse_args()

    # Define dictionary
    frame_groups = {}

    # Loop based on n_sep
    for i in range(args.n_sep):

        # As long as i is not the last iteration, ask for name and frames
        if (i != args.n_sep - 1):

            name = input("Enter name for the separated frames: ")
            bad_list = input("Enter frames number separated by space: ")
            my_list = list(map(int, bad_list.split()))
            frame_groups.update({name: my_list})

        # Otherwise ask for the name for the rest of the frames
        else:
            name = input("Enter name for the rest of the frames: ")
            frame_groups.update({name: []})
    
    # Read science list
    science_list_list = np.atleast_1d(np.loadtxt(args.sciencelist,dtype=str))

    # Separate frames group
    sep_groups = []

    # Loop dictionary
    for j, (name, frames) in enumerate(frame_groups.items()):

        # Create new list name
        new_list_name = args.sciencelist + '_' + name

        # Open new list file
        with open(new_list_name, 'w') as new_list:

            # The first few dictionary items
            if j != args.n_sep:

                # Append the separated frames to a list
                sep_groups.extend(frames)

                for i, file in enumerate(science_list_list):
                    
                    # If frame is in items, write to new list
                    if (i+1) in frames:
                        new_list.write(file + '\n')
            
            # Last dictionary item
            if j == args.n_sep - 1:
    
                for i, file in enumerate(science_list_list):

                    # Write the remaining frames
                    if (i+1) not in sep_groups:
                        new_list.write(file + '\n')



                






