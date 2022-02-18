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
This python script can be used to make an animated gif of PyRate cumulative time series products.

Usage: python3 utils/make_tscuml_animation.py <path to PyRate outdir>
"""
import rasterio
import matplotlib.pyplot as plt
import matplotlib.backend_bases
import numpy as np
import os, sys, re
import xarray as xr
from datetime import datetime as dt
import matplotlib.animation as animation

if len(sys.argv) != 2:
    print('Exiting: Provide abs path to <outdir> as command line argument')
    exit()
else:
    path = sys.argv[1]
    print(f"Looking for PyRate products in: {path}")
#################################

# Reading velocity data
with rasterio.open(os.path.join(path, 'velocity_dir', 'linear_rate.tif')) as src2:
    vel = src2.read()
    bounds2 = src2.bounds
    x_coord2 = np.linspace(bounds2[0], bounds2[2], src2.width)
    y_coord2 = np.linspace(bounds2[1], bounds2[3], src2.height)
    # grab list of time series dates from metadata
    ed = src2.tags()['EPOCH_DATE']

# convert metadata string to list of strings
date_str = re.findall(r'\'(.+?)\'', ed)

# make velocity xarray (also apply deramp to the velocity map directly)
dac2 = xr.DataArray(vel[0,:,:], coords={'lon': x_coord2, 'lat': y_coord2}, dims=['lat', 'lon'])
vs = xr.Dataset()
vs['vel'] = dac2
longitude = vs.coords['lon']
latitude = vs.coords['lat']

# pre-allocate a 3D numpy array to read the tscuml tiff raster bands to
# include first 'zero' time slice, which PyRate does not save to disk
tscuml = np.zeros((len(date_str), vel.shape[1], vel.shape[2])) # commented Chandra

# reading *tif files and generate cumulative variable
print('Reading tscuml files:')
for i, d in enumerate(date_str[1:]):
    print(i+1, 'tscuml_' + d + '.tif')
    with rasterio.open(os.path.join(path, 'timeseries_dir', 'tscuml_' + d + '.tif')) as src:
        data = src.read()
        bounds = src.bounds
        x_coord = np.linspace(bounds[0], bounds[2], src.width)
        y_coord = np.linspace(bounds[1], bounds[3], src.height)
    tscuml[i+1, :, :] = np.squeeze(data, axis=(0,)) # commented chandra

# copy Nans in first time slice to zero epoch
zeroepoch = np.zeros((vel.shape[1], vel.shape[2]))
zeroepoch[np.isnan(tscuml[1,:,:])] = np.nan
tscuml[0,:,:] = zeroepoch

# convert date strings to datetime objects
imdates_dt = [dt.strptime(x, '%Y-%m-%d') for x in date_str]

# make tscuml xarray
dac = xr.DataArray(tscuml, coords={'time': imdates_dt, 'lon': x_coord, 'lat': y_coord}, dims=['time', 'lat', 'lon'])
ds = xr.Dataset()
ds['tscuml'] = dac
n_im, length, width = tscuml.shape

# Add max and min displacement range
#refx1 = int(len(x_coord2)/ 2)
#refx2 = int(len(x_coord2)/ 2) + 1
#refy1 = int(len(y_coord2)/ 2)
#refy2 = int(len(y_coord2)/ 2) + 1
#refvalue_lastepoch = np.nanmean(tscuml[-1, refy1:refy2, refx1:refx2]) # reference values

auto_crange: float = 100
dmin_auto = np.nanpercentile(tscuml, 100 - auto_crange)
dmax_auto = np.nanpercentile(tscuml, auto_crange)
# find the absolute max displacement; round to nearest 10 units, for colour bar limits
lim = np.round(np.amax(np.array([np.abs(dmin_auto), np.abs(dmax_auto)])), decimals=-1)
#dmin = dmin_auto - refvalue_lastepoch
#dmax = dmax_auto - refvalue_lastepoch

####  ANIMATION of LOS Cumulative displacement TS
fig = plt.figure('PyRate Cumulative Displacement Animation', figsize=(5,5))
faxv = fig.add_axes([0.15,0.15,0.75,0.75])
cmap = matplotlib.cm.Spectral_r
cmap.set_bad('grey',1.) # filled grey color to nan value
ims = [] # pre-allocate list for appending time slices

# loop over all timeslices, including zero epoch
for ii in range(0,len(imdates_dt)):
    im = faxv.imshow(ds.tscuml[ii], cmap=cmap, alpha=1, origin='upper',extent=[ds.coords['lon'].min(),
         ds.coords['lon'].max(), ds.coords['lat'].min(), ds.coords['lat'].max()], clim=[-lim, lim])
    title = fig.text(0.40, 0.90, "Date: {}".format(imdates_dt[ii].date()), fontsize=12, va='bottom' )
    ims.append([im, title])

fcbr = fig.colorbar(im, orientation='horizontal')
fcbr.set_label('LOS Displacement [mm]')
ani = animation.ArtistAnimation(fig, ims, interval=500, blit=False)
#plt.show()
file = path + '/timeseries_dir/' +'tscuml_animation.gif'
ani.save(file, writer='imagemagick', fps=10, dpi=100)
print('Animation saved to ' + file)
