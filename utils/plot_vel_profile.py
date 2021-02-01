#!/usr/bin/python3

'''
This script has been copied from LiCSBAS package and modified for PyRate package.
See: https://github.com/yumorishita/LiCSBAS/blob/master/bin/LiCSBAS_plot_ts.py
'''

# plotting PyRate velocity and ts files
import rasterio
import matplotlib.pyplot as plt
import matplotlib.backend_bases
import numpy as np
import os, sys, re
import xarray as xr
from datetime import datetime as dt
from pylab import plot, ginput, show, axis # for velocity profile

if len(sys.argv) != 2:
    print('Exiting: Provide abs path to <outdir> as command line argument')
    exit()
else:
    path = sys.argv[1]
    print(f"Looking for PyRate products in: {path}")

# path = "/g/data/dg9/INSAR_ANALYSIS/EROMANGA/S1/PYRATE/out_15mlk/result_cc06_quadfit_indep_orbit_15mlk_without_DEMerr/"


# ----- Reading linear velocity file --------------------------------
with rasterio.open(os.path.join(path, 'linear_rate.tif')) as src2:
    vel = src2.read()
    bounds2 = src2.bounds
    x_coord2 = np.linspace(bounds2[0], bounds2[2], src2.width)
    y_coord2 = np.linspace(bounds2[1], bounds2[3], src2.height)
    # grab list of time series dates from metadata
    ed = src2.tags()['EPOCH_DATE']

# convert metadata string to list of strings
date_str = re.findall(r'\'(.+?)\'', ed)
imdates_dt = [dt.strptime(x, '%Y-%m-%d') for x in date_str]

# make velocity xarray 
dac2 = xr.DataArray(vel[0,:,:], coords={'lon': x_coord2, 'lat': y_coord2}, dims=['lat', 'lon'])
vs = xr.Dataset()
vs['vel'] = dac2
longitude = vs.coords['lon']
latitude = vs.coords['lat']



# ------ Masking the velocity file using linear sample --------------
# # linear resample *tif file and mask out velocity
src3 = rasterio.open(os.path.join(path, 'linear_samples.tif'))
lsample = src3.read()
lsamp_val = lsample[0,:,:]
[mask_x,mask_y] = np.where((lsamp_val >= (len(imdates_dt)-2)))
vel_test = np.empty((vel.shape[1],vel.shape[2],)) * np.nan
vel_test[mask_x,mask_y] = vs.vel.data[mask_x,mask_y]



# -------- Plot velocity profile across two points -------------------
def vel_profile(line):
    num = 500 
    x, y = np.linspace(line[0][0], line[1][0], num), np.linspace(line[0][1], line[1][1], num)
    ## Extract the values along the line, using nearest-neighbor interpolation
    zi = vel_test[y.astype(int), x.astype(int)]
    return x, y, zi

# vmin = -50; vmax = 50
refx1 = int(len(x_coord2)/ 2)
refx2 = int(len(x_coord2)/ 2) + 1
refy1 = int(len(y_coord2)/ 2)
refy2 = int(len(y_coord2)/ 2) + 1

auto_crange: float = 99.8
refvalue_vel = np.nanmean(vel[0, refy1:refy2 + 1, refx1:refx2 + 1])
vmin_auto = np.nanpercentile(vel[0, :, :], 100 - auto_crange)
vmax_auto = np.nanpercentile(vel[0, :, :], auto_crange)
vmin = vmin_auto - refvalue_vel
vmax = vmax_auto - refvalue_vel

cmap = matplotlib.cm.bwr_r
cmap.set_bad('grey',1.)
fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]}, figsize=(7,10))
cax = axes[0].imshow(vel_test, clim=[vmin, vmax], cmap = cmap)
cbr = fig.colorbar(cax,ax=axes[0], orientation='vertical')
cbr.set_label('Velocity [mm/yr]')

print("Please click any two points")
ig = 1
while ig!=2:
    pts = ginput(2) # it will wait for two clicks
    [x,y,zi] = vel_profile(pts)
    if axes[0].lines:
        del axes[0].lines[0]
        axes[1].cla()
    axes[0].plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], 'ro-') # axes[0].plot([x0, x1], [y0, y1], 'ro-')
    axes[1].plot(x, zi, 'gray')
    axes[1].plot(x, zi, 'r.')
    # axes[1].set_ylim(-5,5)
    axes[1].set_ylabel('LOS velocity [mm/yr]')
    axes[1].set_xlabel('x-axis')
    axes[1].grid(zorder=0)
    plt.pause(0.1)
    fig.canvas.draw()
    fig.canvas.flush_events()

plt.show()
