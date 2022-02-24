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
This python script can be used to make a profile through the PyRate linear rate product.

Usage: python3 utils/plot_linear_rate_profile.py <path to PyRate outdir>
"""
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

# ----- Reading linear velocity file --------------------------------
with rasterio.open(os.path.join(path, 'velocity_dir', 'linear_rate.tif')) as src2:
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
src3 = rasterio.open(os.path.join(path, 'velocity_dir', 'linear_samples.tif'))
lsample = src3.read()
lsamp_val = lsample[0,:,:]
[mask_x,mask_y] = np.where((lsamp_val >= (len(imdates_dt)-2)))
vel_masked = np.empty((vel.shape[1],vel.shape[2],)) * np.nan
vel_masked[mask_x,mask_y] = vs.vel.data[mask_x,mask_y]

# -------- Plot velocity profile across two points -------------------
def get_profile(pts, img):
    '''
    Extract values from image at points along the profile line, using nearest-neighbor interpolation
    '''
    num = 500 # hardcoded 500 points per profile
    # generate 500 points between the provided points
    x, y = np.linspace(pts[0][0], pts[1][0], num), np.linspace(pts[0][1], pts[1][1], num)
    # pull z value from nearest whole pixel
    z = img[y.astype(int), x.astype(int)]
    return x, y, z

# vmin = -50; vmax = 50
#refx1 = int(len(x_coord2)/ 2)
#refx2 = int(len(x_coord2)/ 2) + 1
#refy1 = int(len(y_coord2)/ 2)
#refy2 = int(len(y_coord2)/ 2) + 1
#refvalue_vel = np.nanmean(vel[0, refy1:refy2 + 1, refx1:refx2 + 1])

auto_crange: float = 100
vmin_auto = np.nanpercentile(vel[0, :, :], 100 - auto_crange)
vmax_auto = np.nanpercentile(vel[0, :, :], auto_crange)
# find the absolute max displacement; round to nearest 10 units, for colour bar limits
lim = np.round(np.amax(np.array([np.abs(vmin_auto), np.abs(vmax_auto)])), decimals=-1)
#vmin = vmin_auto - refvalue_vel
#vmax = vmax_auto - refvalue_vel

cmap = matplotlib.cm.Spectral_r
cmap.set_bad('grey',1.)
fig = plt.figure('PyRate Linear Rate Profile', figsize=(7,10))
axes = fig.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]})
cax = axes[0].imshow(vel_masked, clim=[-lim, lim], cmap = cmap)
cbr = fig.colorbar(cax,ax=axes[0], orientation='vertical')
cbr.set_label('LOS velocity [mm/yr]')

print("Please click any two points")
ig = 1
while ig!=2:
    pts = ginput(2) # it will wait for two clicks
    [x,y,z] = get_profile(pts, vel_masked)
    if axes[0].lines:
        del axes[0].lines[0]
        axes[1].cla()
    axes[0].plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], 'ro-') # axes[0].plot([x0, x1], [y0, y1], 'ro-')
    axes[1].plot(x, z, 'gray')
    axes[1].plot(x, z, 'r.')
    # axes[1].set_ylim(-5,5)
    axes[1].set_ylabel('LOS velocity [mm/yr]')
    axes[1].set_xlabel('X coordinate')
    axes[1].grid(zorder=0)
    plt.pause(0.1)
    fig.canvas.draw()
    fig.canvas.flush_events()

plt.show()
