#### Author of this code: James Kirk
#### Contact: jameskirk@live.co.uk 

## Cinta's note:
## Adjust so that if I run this python script in date_folder/raw, the list will be 
#$ saved in date_folder/calib_files

from astropy.io import fits
import glob
import os
import argparse

parser = argparse.ArgumentParser(description='Load all fits files within a directory and produced lists of file types')
parser.add_argument("-c","--clobber",help="Overwrite previously saved lists if they exit. Particularly useful when doing reductions in real-time",action='store_true')
args = parser.parse_args()

all_files = sorted(glob.glob('*.fits'))

# Current working directory
pwd = os.getcwd()

# Parent working directory (one directory above)
parent_dir = os.path.dirname(pwd)

if args.clobber:
     preexisting = [open(i,'w') for i in glob.glob("*_list")]
     [i.close() for i in preexisting]
else:
	pass

def split_list(file_names,pwd):

    for file_number,i in enumerate(file_names):
        f = fits.open(i)
        hdr = f[0].header
        
        try:
            grism = hdr['HIERARCH ESO INS GRIS1 NAME']
            filt = hdr['HIERARCH ESO INS FILT1 NAME']
            slit = hdr['HIERARCH ESO INS SLIT1 NAME']
            obj = hdr['OBJECT']
        except:
            continue
        
        if obj == 'SKY,FLAT':
            obj = 'SKY_FLAT'
        
        if obj == 'WAVE':
            obj = 'ARC'
        
        if obj == 'FOCUS' or obj == 'OTHER':
            continue
        
        if filt == 'Free':
            filt = ''
        else:
            filt = '_'+filt

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
        
        print(i, grism, filt, slit, obj)

        # Create the filename
        filename = obj + '_' + grism + '_' + slit + filt + '_list'
        
        # Output the file one directory above the Raw directory
        filepath = os.path.join(parent_dir, 'calib_files/file_lists', filename)
        
        try:
            file_list = open(filepath,'a')
        except:
            file_list = open(filepath,'w')
        
        file_list.write(pwd+'/'+i+' \n')
        file_list.close()
        
        f.close()

split_list(all_files,pwd)
