#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This script plots the original interferogram, the corresponding correction file, 
and the final corrected interferogram from a PyRate directory with already,
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
         NORMALISE       - Switch to subtract median to from figures (0=YES, 1=No, default is 0). 
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import glob
import re
import math
import argparse
import os
from pyrate.core.shared import Ifg

# Arguments
parser = argparse.ArgumentParser(description="Script to plot correction files with uncorrected and corrected interferogram")
parser.add_argument("IFG_DIR", type=str, help="full path to uncorrected interferograms in PyRate")
parser.add_argument("CORRECTION_FILE_DIR", type=str, help="full path to correction files in PyRate")
parser.add_argument("CORRECTED_IFG_DIR", type=str, help="full path to corrected interferograms in PyRate")
parser.add_argument("SAVE_DIR", type=str, help="full path to directory where images will be saved")
parser.add_argument("FIRST_IFG", type=int, help="first IFG in range of IFGs to plot (e.g. 1 to start plotting at first IFG in directory)")
parser.add_argument("LAST_IFG", type=int, help="last IFG in range of IFGs to plot (e.g. 37 will plot up until the 37th IFG in directory)")
parser.add_argument("-sm", "--subtract_median", type=int, default=1, help="Switch to subtract median to from figures (1=YES, 0=No, default is 1)")
args = parser.parse_args()

# Directories and variable from user agruments
ifg_dir = os.path.abspath(args.IFG_DIR)
corr_dir = os.path.abspath(args.CORRECTION_FILE_DIR)
tempml_dir = os.path.abspath(args.CORRECTED_IFG_DIR)
save_dir = os.path.abspath(args.SAVE_DIR)

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

# define colour map
cmap = cm.Spectral_r
cmap.set_bad(color='grey')

# loop over each ifg in turn
for i in range(args.FIRST_IFG - 1, args.LAST_IFG):    
    # Read data
    orig_ifg = Ifg(ifg_list[i])
    orig_ifg.open()
    orig_ifg.convert_to_mm() # force mm conversion
    ifg = orig_ifg.phase_data
    orig_ifg.close()

    corr =  np.load(corr_list[i])
    
    corr_ifg = Ifg(tempml_list[i])
    corr_ifg.open()
    corr_ifg.convert_to_mm() # force mm conversion
    ifg_corr = corr_ifg.phase_data
    corr_ifg.close()
    
    # Identify Date Pair
    date_pair_list_ifg = re.findall(r'\d{8}-\d{8}', ifg_list[i])
    date_pair_string_ifg = date_pair_list_ifg[0]

    date_pair_list_corr = re.findall(r'\d{8}-\d{8}', corr_list[i])
    date_pair_string_corr = date_pair_list_corr[0]

    date_pair_list_ifgcorr = re.findall(r'\d{8}-\d{8}', tempml_list[i])
    date_pair_string_ifgcorr = date_pair_list_ifgcorr[0]

    # Check the Date-pairs are the same in case of mismatched files saved into the directories
    if date_pair_string_ifg == date_pair_string_corr and date_pair_string_ifg == date_pair_string_ifgcorr:
        print(f'\nPlotting for {date_pair_string_ifg}...\n')
        pass
    else:
        print(f'\nERROR: Interferogram datepair mismatch at {date_pair_string_ifg}, check that directories have the same interferograms.\n')
        break
    
    # PLOTTING
#    ifg_quant = np.nanquantile(ifg, [0.05, 0.95])
#    climit = np.nanmax(np.abs(ifg_quant))
    climit = 30
    
    # subtract median of image
    if args.subtract_median == 1:
        ifg -= np.nanmedian(ifg)
        corr -= np.nanmedian(corr)
        ifg_corr -= np.nanmedian(ifg_corr)

    # three sub-plots
    fig, ax = plt.subplots(1,3, figsize=(6, 3))
     
    # ORIGINAL IFG
    s0 = ax[0].imshow(ifg, cmap=cmap, clim=(-1*climit, climit))
    ax[0].set_title('Original Ifg', fontsize=8)
    ax[0].set_axis_off()

    # CORRECTION FILE
    s1 =ax[1].imshow(corr, cmap=cmap, clim=(-1*climit, climit))
    ax[1].set_title('Correction', fontsize=8)
    ax[1].set_axis_off()

    # CORRECTED IFG
    s2 = ax[2].imshow(ifg_corr, cmap=cmap, clim=(-1*climit,climit))
    ax[2].set_title('Corrected Ifg', fontsize=8)
    ax[2].set_axis_off()

    # Colorbar
    cbar = fig.colorbar(s1, ax=ax.ravel().tolist(), orientation='horizontal', label='mm', \
            ticks=[-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])
    cbar.ax.tick_params(labelsize=8)

    # set figure window colour to grey
    fig.set_facecolor('grey')
    
    # Title
    fig.suptitle(f'{date_pair_string_ifg}', y=0.75, fontsize=10, fontweight='bold')

    plt.savefig(f'{save_dir}/{date_pair_string_ifg}.png', dpi=300)
    
    plt.close()

