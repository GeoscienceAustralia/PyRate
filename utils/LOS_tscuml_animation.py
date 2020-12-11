#!/usr/bin/python3

'''
This script is for generating animated *gif file for the cumulative LOS displacement time-series

'''

# plotting PyRate velocity and ts files
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
with rasterio.open(os.path.join(path, 'linear_rate.tif')) as src2:
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
    with rasterio.open(os.path.join(path, 'tscuml_' + d + '.tif')) as src:
        data = src.read()
        bounds = src.bounds
        x_coord = np.linspace(bounds[0], bounds[2], src.width)
        y_coord = np.linspace(bounds[1], bounds[3], src.height)
    tscuml[i+1, :, :] = np.squeeze(data, axis=(0,)) # commented chandra

# convert date strings to datetime objects
imdates_dt = [dt.strptime(x, '%Y-%m-%d') for x in date_str]

# make tscuml xarray
dac = xr.DataArray(tscuml, coords={'time': imdates_dt, 'lon': x_coord, 'lat': y_coord}, dims=['time', 'lat', 'lon'])
ds = xr.Dataset()
ds['tscuml'] = dac
n_im, length, width = tscuml.shape

# choose final time slice
time_slice = len(imdates_dt)-1

####  ANIMATION of LOS Cumulative displacement TS
fig = plt.figure('Animation Displacement', figsize=(5,5))
faxv = fig.add_axes([0.15,0.15,0.75,0.75])
cmap = matplotlib.cm.bwr_r #
cmap.set_bad('grey',1.) # filled grey color to nan value
ims = []
for ii in range(1,time_slice+1):
    print(ii)
    im = faxv.imshow(ds.tscuml[ii], cmap=cmap, alpha=1, origin='upper',extent=[ds.coords['lon'].min(), ds.coords['lon'].max(), ds.coords['lat'].min(), ds.coords['lat'].max()], clim=(-50, 50)) #for displacement
    title = fig.text(0.55, 0.90, "Date: {}".format(imdates_dt[ii].date()), fontsize=8, va='bottom' )
    ims.append([im, title])
fcbr = fig.colorbar(im, orientation='horizontal')
fcbr.set_label('LOS Displacement [mm]')
ani = animation.ArtistAnimation(fig, ims, interval=2, blit=False)
plt.show()
ani.save(path + 'Animation.gif', writer='imagemagick', fps=10, dpi=100)
print('Animation Done!')
#######
