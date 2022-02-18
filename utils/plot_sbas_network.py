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
This python script can be used to generate a baseline-time plot of the 
interferograms used in the PyRate SBAS (Small Baseline Subset) network.
The functions 'plot_baseline_time_sbas' and 'epoch_baselines' are copies of the corresponding
functions in 'calc_baselines_functions.py' in GA's gamma-insar repository, see
https://github.com/GeoscienceAustralia/gamma_insar/blob/develop/calc_baselines_functions.py

Usage: python3 utils/plot_sbas_network.py <path to PyRate outdir> 
"""
import rasterio
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os, sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime, timedelta


print('')
if len(sys.argv) != 2:
    print('Exiting: Provide path to <PyRate outdir> as command line argument')
    print('')
    print('Usage: python3 utils/plot_time_series.py <path to PyRate outdir>')
    exit()
else:
    path = sys.argv[1]
    print(f"Looking for PyRate products in: {path}")


def readtif(tifname: str):
    """
    wrapper for rasterio tif reading
    """
    print(f"Reading file: {tifname}")
    with rasterio.open(tifname) as src:
        md = src.tags()

    return md


def plot_baseline_time_sbas(epochs, Bperps, epoch1, epoch2, filename):
    """
    Make a baseline time plot including IFG connections and save to disk
    """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    divider = make_axes_locatable(ax1)

    # plot interferograms as lines
    for n, m in zip(epoch1, epoch2):
        #print n, m
        x = [epochs[n], epochs[m]]
        y = [Bperps[n], Bperps[m]]
        #  baselines[x]
        ax1.plot_date(x, y, xdate=True, ydate=False, linestyle='-',
        color = 'r', linewidth=1.0)

    # plot epochs as filled circles
    ax1.plot_date(epochs, Bperps, xdate=True, ydate=False, marker="o",
                  markersize=14, markerfacecolor="black", linestyle="None")

    # plot epoch numbers as symbols
    labels = [i+1 for i in range(len(Bperps))]
    for a, b, c in zip(epochs, Bperps, labels):
        ax1.text(a, b, c, color="white", ha="center", va="center", size=9,
                 weight="bold")

    #format the time axis ticks
    years   = mdates.MonthLocator(bymonth=[1, 7])   # every 0.5 year
    months   = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter("%Y-%m-%d")
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_major_formatter(yearsFmt)
    ax1.xaxis.set_minor_locator(months)

    #set the time axis range
    date_min = epochs.min()
    date_max = epochs.max()
    date_range = date_max - date_min
    date_add = date_range.days/15
    ax1.set_xlim(date_min - timedelta(days=date_add), date_max +
                 timedelta(days=date_add))

    # set the Bperp axis range
    Bperp_min = min(Bperps)
    Bperp_max = max(Bperps)
    Bperp_range = Bperp_max - Bperp_min
    ax1.set_ylim(Bperp_min - Bperp_range/15, Bperp_max + Bperp_range/15)

    #set axis titles
    ax1.set_xlabel("Date (YYYY-MM-DD)")
    ax1.set_ylabel("Perpendicular Baseline (m)")
    ax1.grid(True)

    #rotates and right aligns the date labels
    fig.autofmt_xdate()

    # Save plot to PNG file
    plt.savefig(filename, orientation="landscape", transparent=False,
                format="png")
    return


def epoch_baselines(epochs, bperp, masidx, slvidx, supermaster):
    '''
    Determine relative perpendicular baselines of epochs from
    interferometric baselines

    INPUT:
    epochs    list of epoch dates
    bperp    list of interferogram absolute perpendicular baselines
    masidx    list of master indices from get_index()
    slvidx    list of slave indices from get_index()
    supermaster    epoch to set relative bperp to zero (integer)

    OUTPUT:
    epochbperp    list of epoch relative perpendicular baselines
    '''

    # Count number of ifgs and epochs
    nifgs = len(bperp)
    nepochs = len(epochs)
    print(nifgs, "interferograms and", nepochs, "epochs in the network.")

    # Initialise design matrix 'A'
    A = np.zeros((nifgs+1,nepochs))

    # assign super-master epoch to constrain relative baselines
    A[0,supermaster] = 1
    b = np.zeros(nifgs+1)
    b[1:nifgs+1] = bperp

    # Construct design matrix
    for i in range(nifgs):
        imas = masidx[i]
        islv = slvidx[i]
        A[i+1,imas] = -1
        A[i+1,islv] = 1

    # Do overdetermined linear inversion x=A\b
    x = np.linalg.lstsq(A, b, rcond=None)
    return x[:][0]


###########
# Main code
print('')
print('Generating the baseline-time plot from tif-files in temp_mlooked_dir')

# reading metadata from tif file
path_to_tif = os.path.join(path, 'temp_mlooked_dir/*.tif')
# some empty lists
Bperps_ifg = []
epoch1 = []
epoch2 = []
for tif_file in glob.glob(path_to_tif):
    md = readtif(tif_file)
    # look for perpendicular baseline in geotiff metadata; skip ifg if not there
    if 'BASELINE_PERP_METRES' in md:
        Bperps_ifg.append(float(md['BASELINE_PERP_METRES']))
        epoch1.append(md['FIRST_DATE'])
        epoch2.append(md['SECOND_DATE'])
print('')

# Quit if no ifg has a perpendicular baseline value
if len(Bperps_ifg) == 0:
    print('No perpendicular baseline values in ifg metadata. ' +
          'First run the DEM error correction to calculate baselines.')
    quit()

# create date vector containing all epochs in the network
epochs= list(set(epoch1+epoch2))
epochs.sort()

# get epoch indices for all interferograms
epoch1_ix = []
for e in epoch1:
    epoch1_ix.append(epochs.index(e))
epoch2_ix = []
for e in epoch2:
    epoch2_ix.append(epochs.index(e))
# convert epochs to datetime object
epochs[:] = [datetime.strptime(epo, "%Y-%m-%d") for epo in epochs]

# convert IFG Bperps into single-master Bperps referenced to the first epoch
Bperps_epoch = epoch_baselines(epochs,Bperps_ifg,epoch1_ix,epoch2_ix, 0)

# filename of output image
filename = os.path.join(path, 'temp_mlooked_dir/baseline_time_plot.png')
# call the function to create the plot
plot_baseline_time_sbas(np.array(epochs), Bperps_epoch, epoch1_ix, epoch2_ix, filename)
print('Network plot saved to ' + filename)
print('')

