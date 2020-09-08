#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
This python script can be used to plot PyRate linear-rate and cumulative time series
products.
This script is a version of a script in the LiCSBAS package; see
https://github.com/yumorishita/LiCSBAS/blob/master/bin/LiCSBAS_plot_ts.py

Usage: python3 utils/plot_time_series.py <path to PyRate outdir> 
"""
import rasterio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, RectangleSelector, CheckButtons
import matplotlib.dates as mdates
import matplotlib.backend_bases
import numpy as np
import fnmatch
import os, sys, re
import statsmodels.api as sm
import xarray as xr
from datetime import datetime as dt
import warnings

if len(sys.argv) != 2:
    print('Exiting: Provide path to <PyRate outdir> as command line argument')
    print('')
    print('Usage: python3 utils/plot_time_series.py <path to PyRate outdir>')
    exit()
else:
    path = sys.argv[1]
    print(f"Looking for PyRate products in: {path}")


###############################
def readtif(tifname: str):
    """
    wrapper for rasterio tif reading
    """
    print(f"Reading file: {tifname}")
    with rasterio.open(tifname) as src:
        img = src.read()
        bounds = src.bounds
        md = src.tags()

    x_coord = np.linspace(bounds[0], bounds[2], src.width)
    y_coord = np.linspace(bounds[1], bounds[3], src.height)

    return img, x_coord, y_coord, md
###############################

# reading velocity data from linear_rate product
vel, x_coord, y_coord, md = readtif(os.path.join(path, 'linear_rate.tif'))

# read regression intercept from linear_intercept product
intercept, _, _, _ = readtif(os.path.join(path, 'linear_intercept.tif'))

# convert time series dates from metadata string to list of strings
date_str = re.findall(r'\'(.+?)\'', md['EPOCH_DATE'])

# convert date strings to datetime objects
imdates_dt = [dt.strptime(x, '%Y-%m-%d') for x in date_str]
imdates_ordinal = [x.toordinal() for x in imdates_dt]

# make velocity xarray
dac2 = xr.DataArray(vel[0,:,:], coords={'lon': x_coord, 'lat': y_coord}, dims=['lat', 'lon'])
vs = xr.Dataset()
vs['vel'] = dac2
longitude = vs.coords['lon']
latitude = vs.coords['lat']

# pre-allocate a 3D numpy array to read the tscuml tiff raster bands to
# include first 'zero' time slice, which PyRate does not save to disk
tscuml = np.zeros((len(date_str), vel.shape[1], vel.shape[2]))

# reading tscuml*tif files and add to tscuml variable
for i, d in enumerate(date_str[1:]):
    data, x_coord, y_coord, _ = readtif(os.path.join(path, 'tscuml_' + d + '.tif'))
    tscuml[i+1, :, :] = np.squeeze(data, axis=(0,))

# make tscuml xarray
dac = xr.DataArray(tscuml, coords={'time': imdates_dt, 'lon': x_coord, 'lat': y_coord}, dims=['time', 'lat', 'lon'])
ds = xr.Dataset()
ds['tscuml'] = dac
n_im, length, width = tscuml.shape

# choose final time slice
time_slice = len(imdates_dt)-1

# set reference area and initial point (scene centre)
refx1 = int(len(x_coord)/ 2)
refx2 = int(len(x_coord)/ 2) + 1
refy1 = int(len(y_coord)/ 2)
refy2 = int(len(y_coord)/ 2) + 1
refarea = (refx1,refx2,refy1,refy2)
point_x = refx1
point_y = refy1


def get_range(arr, refarea):
    """
    determine plot scale range
    """
    auto_crange: float = 99.8
    refvalue = np.nanmean(arr[refarea[0]:refarea[1], refarea[2]:refarea[3]]) # reference values
    if str(refvalue) == 'nan':
        refvalue = 0
    dmin_auto = np.nanpercentile((tscuml[-1, :, :]), 100 - auto_crange)
    dmax_auto = np.nanpercentile((tscuml[-1, :, :]), auto_crange)
    dmin = dmin_auto - refvalue
    dmax = dmax_auto - refvalue

    return dmin, dmax

# range from last tscuml epoch
dmin, dmax = get_range(tscuml[-1, :, :], refarea)

# range from velocity
vmin, vmax = get_range(vel[0, :, :], refarea)

# Plot figure of Velocity and Cumulative displacement
figsize = (7,7)
pv = plt.figure('PyRate: linear_rate / tscuml map viewer', figsize)
axv = pv.add_axes([0.15,0.15,0.75,0.83])
axv.set_title('linear_rate')
axt2 = pv.text(0.01, 0.99, 'Left-doubleclick:\n Plot time series\nRight-drag:\n Change ref area', fontsize=8, va='top')
axt = pv.text(0.01, 0.78, 'Ref area:\n X {}:{}\n Y {}:{}\n (start from 0)'.format(refx1, refx2, refy1, refy2),
              fontsize=8, va='bottom')

# create masked array for NaNs
mvel = np.ma.array(vs.vel, mask=np.isnan(vs.vel))
cmap = matplotlib.cm.bwr_r
cmap.set_bad('grey')
cax = axv.imshow(mvel, clim=[vmin, vmax], cmap=cmap)
cbr = pv.colorbar(cax, orientation='vertical')
cbr.set_label('mm/yr')

# Radio buttom for velocity selection
mapdict_data = {}
mapdict_unit = {}
mapdict_vel = {'Vel': vel}
mapdict_unit.update([('Vel', 'mm/yr')])
mapdict_data = mapdict_vel

axrad_vel = pv.add_axes([0.01, 0.3, 0.13, len(mapdict_data)*0.025+0.04])
# Radio buttons
radio_vel = RadioButtons(axrad_vel, tuple(mapdict_data.keys()))
for label in radio_vel.labels:
    label.set_fontsize(10)

tscuml_disp_flag = False
climauto = True

def line_select_callback(eclick, erelease):
    """
    Set ref function
    """
    global refx1, refx2, refy1, refy2, dmin, dmax   ## global cannot change existing values... why?

    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    if x1 <= x2:
        refx1, refx2 = [int(np.round(x1)), int(np.round(x2))]
    elif x1 > x2:
        refx1, refx2 = [int(np.round(x1)), int(np.round(x2))]
    if y1 <= y2:
        refy1, refy2 = [int(np.round(y1)), int(np.round(y2))]
    elif y1 > y2:
        refy1, refy2 = [int(np.round(y1)), int(np.round(y2))]
    refarea = (refx1,refx2,refy1,refy2)

    axt.set_text('Ref area:\n X {}:{}\n Y {}:{}\n (start from 0)'.format(refx1, refx2, refy1, refy2))
    pv.canvas.draw()

    ### Change clim
    if climauto:  ## auto
        dmin, dmax = get_range(tscuml[-1, :, :], refarea)

    ### Update draw
    if not tscuml_disp_flag:  ## vel or noise indice # Chandra
        val_selected = radio_vel.value_selected
        val_ind = list(mapdict_data.keys()).index(val_selected)
        radio_vel.set_active(val_ind)
    else:  ## cumulative displacement
        time_selected = tslider.val
        tslider.set_val(time_selected)

    if lastevent:  ## Time series plot
        printcoords(lastevent)

RS = RectangleSelector(axv, line_select_callback, drawtype='box', useblit=True, button=[3], spancoords='pixels', interactive=False)

plt.connect('key_press_event', RS)


vlimauto = True
def show_vel(val_ind):
    global vmin, vmax, tscuml_disp_flag
    tscuml_disp_flag = False

    if 'Vel' in val_ind:  ## Velocity
        data = mapdict_data[val_ind]
        if vlimauto:  ## auto
            vmin = np.nanpercentile(data, 100 - auto_crange)
            vmax = np.nanpercentile(data, auto_crange)
        cax.set_cmap(cmap)
        cax.set_clim(vmin, vmax)
        cbr.set_label('mm/yr')

    cbr.set_label(mapdict_unit[val_ind])
    cax.set_data(data)
    axv.set_title(val_ind)
    #cax.set_clim(-100, 100)

    pv.canvas.draw()

radio_vel.on_clicked(show_vel)

# Slider for cumulative displacement
axtim = pv.add_axes([0.1, 0.08, 0.8, 0.05], yticks=[])
mdt = mdates.date2num(imdates_dt)
tslider = Slider(axtim, 'year', mdates.date2num(imdates_dt[0]), mdates.date2num(imdates_dt[-1]))
tslider.ax.bar(mdt, np.ones(len(mdt)), facecolor='black', width=4)
tslider.ax.bar(mdt[0], 1, facecolor='red', width=12)

loc_tslider = tslider.ax.xaxis.set_major_locator(mdates.AutoDateLocator())
try:  # Only support from Matplotlib 3.1!
    tslider.ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc_tslider))
    for label in tslider.ax.get_xticklabels():
        label.set_rotation(20)
        label.set_horizontalalignment('right')
except:
    tslider.ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    for label in tslider.ax.get_xticklabels():
        label.set_rotation(20)
        label.set_horizontalalignment('right')

dstr_ref = imdates_dt[0].strftime('%Y/%m/%d') # reference image date
def tim_slidupdate(val):
    global tscuml_disp_flag
    timein = tslider.val
    timenearest = np.argmin(np.abs(mdates.date2num(imdates_dt) - timein))

    dstr = imdates_dt[timenearest].strftime('%Y/%m/%d')
    axv.set_title('tscuml: %s (Ref: %s)' % (dstr, dstr_ref))

    newv = (tscuml[timenearest, :, :])
    cax.set_data(newv)
    cax.set_cmap(cmap)
    cax.set_clim(dmin, dmax)
    cbr.set_label('mm')
    tscuml_disp_flag = True

    pv.canvas.draw()

tslider.on_changed(tim_slidupdate)

##### Plot figure of time-series at a point
pts = plt.figure('PyRate: pixel time-series graph', (9,5))
axts = pts.add_axes([0.12, 0.14, 0.7, 0.8])

axts.scatter(imdates_dt, np.zeros(len(imdates_dt)), c='b', alpha=0.6)
axts.grid()

axts.set_xlabel('Time [Year]')
axts.set_ylabel('Displacement [mm]')

loc_ts = axts.xaxis.set_major_locator(mdates.AutoDateLocator())
try:  # Only support from Matplotlib 3.1
    axts.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc_ts))
except:
    axts.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    for label in axts.get_xticklabels():
        label.set_rotation(20)
        label.set_horizontalalignment('right')

### Ref info at side
axtref = pts.text(0.83, 0.95, 'Ref area:\n X {}:{}\n Y {}:{}\n (start from 0)\nRef date:\n {}'.format(refx1, refx2, refy1, refy2, dstr_ref), fontsize=8, va='top')

### Fit function for time series
fitbox = pts.add_axes([0.83, 0.10, 0.16, 0.25])
models = ['PyRate linear rate', 'Annual+L model', 'Quad model', 'Annual+Q model']
visibilities = [True, False, False, False]
fitcheck = CheckButtons(fitbox, models, visibilities)
for label in fitcheck.labels:
    label.set_fontsize(8)

def fitfunc(label):
    index = models.index(label)
    visibilities[index] = not visibilities[index]
    lines1[index].set_visible(not lines1[index].get_visible())

    pts.canvas.draw()

fitcheck.on_clicked(fitfunc)

### First show of selected point in image window
pax, = axv.plot([point_y], [point_x], 'k', linewidth=3)
pax2, = axv.plot([point_y], [point_x], 'Pk')

### Plot time series at clicked point
lastevent = []
geocod_flag = False
label1 = 'tscuml'
label2 = 'PyRate linear_rate'

ylen = []


def calc_model(dph, imdates_ordinal, xvalues, model, vel1p, intercept1p):
    """
    Function to calculate model to fit cumulative time series data for one pixel
    """
    imdates_years = np.divide(imdates_ordinal,365.25) # imdates_ordinal / 365.25
    xvalues_years = xvalues / 365.25

    A = sm.add_constant(imdates_years)  # [1, t]
    An = sm.add_constant(xvalues_years)  # [1, t]
    if model == 0:  # PyRate linear fit results
        # print(" M = ", vel1p, " & C = ", intercept1p )
        A = np.dot(A, vel1p) + intercept1p
        An = np.dot(An, vel1p) + intercept1p
        # pass
    if model == 1:  # Annual+L
        sin = np.sin(2 * np.pi * imdates_years)
        cos = np.cos(2 * np.pi * imdates_years)
        A = np.concatenate((A, sin[:, np.newaxis], cos[:, np.newaxis]), axis=1)
        sin = np.sin(2 * np.pi * xvalues_years)
        cos = np.cos(2 * np.pi * xvalues_years)
        An = np.concatenate((An, sin[:, np.newaxis], cos[:, np.newaxis]), axis=1)
    if model == 2:  # Quad
        A = np.concatenate((A, (imdates_years ** 2)[:, np.newaxis]), axis=1)
        An = np.concatenate((An, (xvalues_years ** 2)[:, np.newaxis]), axis=1)
    if model == 3:  # Annual+Q
        sin = np.sin(2 * np.pi * imdates_years)
        cos = np.cos(2 * np.pi * imdates_years)
        A = np.concatenate((A, (imdates_years ** 2)[:, np.newaxis], sin[:, np.newaxis], cos[:, np.newaxis]), axis=1)
        sin = np.sin(2 * np.pi * xvalues_years)
        cos = np.cos(2 * np.pi * xvalues_years)
        An = np.concatenate((An, (xvalues_years ** 2)[:, np.newaxis], sin[:, np.newaxis], cos[:, np.newaxis]), axis=1)

    result = sm.OLS(dph, A, missing='drop').fit()
    return result.predict(An)
#################################

def printcoords(event):
    global dph, lines1, lines2, lastevent, imdates_ordinal, imdates_dt
    # outputting x and y coords to console
    if event.inaxes != axv:
        return
    elif event.button != 1:
        return
    elif not event.dblclick:  ## Only double click
        return
    else:
        lastevent = event  ## Update last event

    ii = np.int(np.round(event.ydata))
    jj = np.int(np.round(event.xdata))

    ### Plot on image window
    ii1h = ii - 0.5; ii2h = ii + 1 - 0.5
    jj1h = jj - 0.5; jj2h = jj + 1 - 0.5

    pax.set_data([jj1h, jj2h, jj2h, jj1h, jj1h], [ii1h, ii1h, ii2h, ii2h, ii1h])
    pax2.set_data(jj, ii)
    pv.canvas.draw()

    axts.cla()
    axts.grid(zorder=0)
    axts.set_axisbelow(True)
    axts.set_xlabel('Date')
    axts.set_ylabel('Cumulative Displacement (mm)')

    ### Get values of noise indices and incidence angle
    noisetxt = ''
    for key in mapdict_data:
        # val_temp = mapdict_data[key]
        val = mapdict_data[key][0, ii, jj]
        # val = val_temp[0,ii,jj]
        unit = mapdict_unit[key]
        if key.startswith('Vel'):
            continue
        elif key.startswith('n_') or key == 'mask':
            noisetxt = noisetxt + '{}: {:d} {}\n'.format(key, int(val), unit)
        else:
            noisetxt = noisetxt + '{}: {:.2f} {}\n'.format(key, int(val), unit)

    try:  # Only support from Matplotlib 3.1!
        axts.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc_ts))
    except:
        axts.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
        for label in axts.get_xticklabels():
            label.set_rotation(20)
            label.set_horizontalalignment('right')

    ### If not masked
    ### tscuml file
    vel1p = vel[0, ii, jj]
    intercept1p = intercept[0,ii,jj]
    dph = tscuml[:,ii,jj]

    ## fit function
    lines1 = [0, 0, 0, 0]
    xvalues = np.arange(imdates_ordinal[0], imdates_ordinal[-1], 10)
    xdates = [dt.fromordinal(pp) for pp in xvalues]
    # Mask to exclude nan elements
    mask = ~np.isnan(dph)
    # remove nan elements from both arrays
    imo = np.asarray(imdates_ordinal)
    imo = imo[mask]
    imt = np.asarray(imdates_dt)
    imt = imt[mask]
    dph = dph[mask]

    for model, vis in enumerate(visibilities):
        if len(dph) > 1:
            yvalues = calc_model(dph, imo, xvalues, model, vel1p, intercept1p)
            if model == 0:
                lines1[model], = axts.plot(xdates, yvalues, 'r-', label=label2, visible=vis, alpha=0.6, zorder=3)
                axts.legend()
            else:
                lines1[model], = axts.plot(xdates, yvalues, 'g-', visible=vis, alpha=0.6, zorder=3)

    axts.scatter(imt, dph, label=label1, c='b', alpha=0.6, zorder=5)
    axts.set_title('Velocity = {:.1f} mm/yr @({}, {})'.format(vel1p, jj, ii), fontsize=10)

    ### Y axis
    if ylen:
        vlim = [np.nanmedian(dph) - ylen / 2, np.nanmedian(dph) + ylen / 2]
        axts.set_ylim(vlim)

    ### Legend
    axts.legend()

    pts.canvas.draw()


#%% First show of time series window
event = matplotlib.backend_bases.LocationEvent
event.xdata = point_x
event.ydata = point_y
event.inaxes = axv
event.button = 1
event.dblclick = True
lastevent = event
printcoords(lastevent)


#%% Final linking of the canvas to the plots.
cid = pv.canvas.mpl_connect('button_press_event', printcoords)
with warnings.catch_warnings(): ## To silence user warning
    warnings.simplefilter('ignore', UserWarning)
    plt.show()
pv.canvas.mpl_disconnect(cid)

