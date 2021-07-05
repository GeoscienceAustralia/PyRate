import numpy as np
from matplotlib import pyplot as plt
import glob
import re
import math
import rasterio as rio
import argparse
import os

"""
This script plots the original interferogram, the corresponding correction file, 
and the resulting corrected interferogram from a PyRate directory with already,
processed data. directories are given as user arguments to the script, and the
number of plots is determined by a number range given by user. 

Usage: python3 plot_correction_files.py <IFG_DIR> <CORRECTION_DIR> <CORRECTED_DIR> <SAVE_DIR> <FIRST_IFG> <LAST_IFG>

Command-line arguments:
           IFG_DIR       - full path to uncorrected interferograms in PyRate.
    CORRECTION_DIR       - full path to correction files in PyRate.
     CORRECTED_DIR       - full path to corrected interferograms in PyRate.
          SAVE_DIR       - full path to directory where images will be saved (needs to exist).
         FIRST_IFG       - first IFG in range of IFGs to plot (e.g. 1 to start plotting at 1st IFG in directory).
          LAST_IFG       - last IFG in range of IFGs to plot (e.g. 37 will plot up until the 37th IFG in directory).
"""

# Arguments
parser = argparse.ArgumentParser(description="Script to plot correction files with uncorrected and corrected interferogram")
parser.add_argument("IFG_DIR", type=str, help="full path to uncorrected interferograms in PyRate")
parser.add_argument("CORRECTION_FILE_DIR", type=str, help="full path to correction files in PyRate")
parser.add_argument("CORRECTED_IFG_DIR", type=str, help="full path to corrected interferograms in PyRate")
parser.add_argument("SAVE_DIR", type=str, help="full path to directory where images will be saved")
parser.add_argument("FIRST_IFG", type=int, help="first IFG in range of IFGs to plot (e.g. 1 to start plotting at first IFG in directory)")
parser.add_argument("LAST_IFG", type=int, help="last IFG in range of IFGs to plot (e.g. 37 will plot up until the 37th IFG in directory)")
args = parser.parse_args()



# Directories and variable from user agruments
ifg_dir = os.path.abspath(args.IFG_DIR)
corr_dir = os.path.abspath(args.CORRECTION_FILE_DIR)
tempml_dir = os.path.abspath(args.CORRECTED_IFG_DIR)
save_dir = os.path.abspath(args.SAVE_DIR)

first_ifg_num = args.FIRST_IFG - 1
last_ifg_num = args.LAST_IFG - 1


# Create Lists
ifg_list = []
for file in glob.glob(f'{ifg_dir}/*ifg.tif'):
    ifg_list.append(file)

corr_list = []
for file in glob.glob(f'{corr_dir}/*.npy'):
    corr_list.append(file)

tempml_list = []
for file in glob.glob(f'{tempml_dir}/*ifg.tif'):
    tempml_list.append(file)

# Sort
ifg_list.sort()
corr_list.sort()
tempml_list.sort()

for i in range(first_ifg_num, last_ifg_num + 1):    
    
    # Read data
    with rio.open(ifg_list[i]) as src:
        ifg = src.read(1)
        ifg[ifg==0.0] = np.nan
        mask = np.isnan(ifg)
    
        # convert to mm 
        ifg = ifg * 1000 * (0.0562356424 / (4 * math.pi))
    
    corr =  np.load(corr_list[i])
    corr[mask] = np.nan
    

    with rio.open(tempml_list[i]) as src:
        ifg_corr = src.read(1)
    
    # Identify Date Pair
    date_pair_list = re.findall(r'\d{8}-\d{8}', ifg_list[i])
    date_pair_string = date_pair_list[0]
    
    print(f'\nPlotting for {date_pair_string}...\n')
    
    # Plot
    climit = 100
    fig, ax = plt.subplots(1,3, figsize=(6, 3))
     
    # IFG
    s0 = ax[0].imshow(ifg-np.nanmedian(ifg), cmap='bwr', clim=(-1*climit, climit))
    ax[0].set_axis_off()

    # CORRECTION FILE
    s1 =ax[1].imshow(corr-np.nanmedian(corr), cmap='bwr', clim=(-1*climit, climit))
    ax[1].set_axis_off()

    # IFG CORRECTED
    s2 = ax[2].imshow(ifg_corr-np.nanmedian(ifg_corr), cmap='bwr', clim=(-1*climit,climit))
    ax[2].set_axis_off()

    # Extra
    fig.colorbar(s0, ax=ax[0], location='bottom', label='mm')
    fig.colorbar(s1, ax=ax[1], location='bottom', label='mm')
    fig.colorbar(s2, ax=ax[2], location='bottom', label='mm')

    fig.set_facecolor('grey')
    
    fig.tight_layout()

    # Title
    ax[1].set_title(f'{date_pair_string}', fontsize=10, fontweight='bold')

    plt.savefig(f'{save_dir}/{date_pair_string}.png', dpi=300)
    
    plt.close()

    i = i + 1

