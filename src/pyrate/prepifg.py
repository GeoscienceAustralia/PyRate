"""
Created on 23/10/2012
@author: Ben Davies
"""

import os, sys
from os.path import join, splitext
from math import modf
from numbers import Number
from tempfile import mkstemp
from itertools import product
from subprocess import check_call

from scipy.stats.stats import nanmean
from numpy import array, where, nan, isnan, mean, float32, zeros

from shared import Ifg, DEM
from roipac import filename_pair, write_roipac_header
from config import parse_namelist
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST
from config import IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST, DEM_FILE
from config import PROJECTION_FLAG, ORBITAL_FIT_LOOKS_X, ORBITAL_FIT_LOOKS_Y
from ifgconstants import X_FIRST, Y_FIRST, X_LAST, Y_LAST, X_STEP, Y_STEP
from ifgconstants import WIDTH, FILE_LENGTH, WAVELENGTH
from ifgconstants import Z_OFFSET, Z_SCALE, PROJECTION, DATUM

# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP]

GRID_TOL = 1e-6



def prepare_ifgs(params, thresh=0.5, use_exceptions=False, verbose=False):
	"""Produces multilooked/resampled data files for PyRate analysis.
	params - dict of named values from pyrate config file
	threshhold - 0.0->1.0 controls NaN handling when resampling to coarser grids,
	             value is proportion above which the number of NaNs in an area is
	             considered invalid. threshhold=0 resamples to NaN if 1 or more
	             contributing cells are NaNs. At 0.25, it resamples to NaN if 1/4
	             or more contributing cells are NaNs. At 1.0, areas are resampled
	             to NaN only if all area cells are NaNs.
	use_exceptions - TODO: PROB REMOVE
	verbose - controls level of gdalwarp output
	"""
	# validate config file settings
	crop_opt = params[IFG_CROP_OPT]
	if crop_opt not in CROP_OPTIONS:
		raise PreprocessingException("Unrecognised crop option: %s" % crop_opt)

	check_looks(params)
	paths = [join(params[OBS_DIR], p) for p in parse_namelist(params[IFG_FILE_LIST])]
	ifgs = [Ifg(p) for p in paths]

	# treat DEM as Ifg as core API is equivalent
	if params.has_key(DEM_FILE):
		ifgs.append(DEM(params[DEM_FILE]))

	check_resolution(ifgs)
	for i in ifgs:
		i.open()

	# calculate crop option for passing to gdalwarp
	if crop_opt == MINIMUM_CROP:
		xmin, ymin, xmax, ymax = min_bounds(ifgs)
	elif crop_opt == MAXIMUM_CROP:
		xmin, ymin, xmax, ymax = max_bounds(ifgs)

	elif crop_opt == CUSTOM_CROP:
		xmin, xmax = params[IFG_XFIRST], params[IFG_XLAST]
		ymin, ymax = params[IFG_YLAST], params[IFG_YFIRST]
		check_crop_coords(ifgs, xmin, xmax, ymin, ymax, use_exceptions)

	# calculate args for gdalwarp reprojection
	extents = [str(s) for s in (xmin, ymin, xmax, ymax) ]
	resolution = None
	if params[IFG_LKSX] > 1 or params[IFG_LKSY] > 1:
		resolution = [params[IFG_LKSX] * i.X_STEP, params[IFG_LKSY] * i.Y_STEP ]

	# Generate gdalwarp args for each ifg, some hackery required to interface with
	# cmd line interface. For resampling, gdalwarp is called 2x, once to subset
	# the source data for Pirate's form of averaging/resampling, the second to
	# generate the final dataset with correct extents/shape/cell count. Without
	# resampling, gdalwarp is only needed to cut out the required segment.

	orb = _do_orbital_multilooking(params)
	xl, yl = params[IFG_LKSX], params[IFG_LKSY]
	for i in ifgs:
		lyr = warp(i, xl, yl, extents, resolution, thresh, verbose)

		# handle 2nd round resampling for orbital correction data
		if(orb):
			xl2, yl2 = params[ORBITAL_FIT_LOOKS_X], params[ORBITAL_FIT_LOOKS_Y]
			xres = resolution[0] * xl2
			yres = resolution[1] * yl2
			warp(lyr, xl2, yl2, extents, [xres, yres], thresh, verbose)


def _do_orbital_multilooking(params):
	'''Returns True if params dict has keys/values for orbital multilooking'''
	keys = [ORBITAL_FIT_LOOKS_X, ORBITAL_FIT_LOOKS_Y]
	if all([params.has_key(k) for k in keys]):
		if all([params[k] > 1 for k in keys]):
			return True
	return False


def _file_ext(raster):
	'''Returns file ext string based on type of raster.'''
	if isinstance(raster, Ifg):
		return "tif"
	elif hasattr(raster, DATUM):
		return "dem"
	else:
		raise NotImplementedError("TODO: additional raster types?")


def _resample_ifg(ifg, cmd, x_looks, y_looks, thresh):
	'''Convenience function to resample data from a given Ifg (more coarse).'''

	# HACK: create tmp ifg, extract data array for manual resampling as gdalwarp
	# lacks Pirate's averaging method
	tmp_path = mkstemp()[1]
	check_call(cmd + [ifg.data_path, tmp_path])
	tmp = type(ifg)(tmp_path, ifg.hdr_path) # dynamically handle Ifgs & Rasters
	tmp.open()

	if isinstance(ifg, Ifg):
		# TODO: resample amplitude band too?
		data = tmp.phase_band.ReadAsArray()
		data = where(data == 0, nan, data) # flag incoherent cells as NaNs
	elif isinstance(ifg, DEM):
		data = tmp.height_band.ReadAsArray()
	else:
		raise NotImplementedError("TODO: other raster types to handle?")

	del tmp # manual close
	os.remove(tmp_path)
	return resample(data, x_looks, y_looks, thresh)


def warp(ifg, x_looks, y_looks, extents, resolution, thresh, verbose):
	'''TODO
	xlooks - integer factor to scale X axis by, 5 is 5x smaller, 1 is no change.
	ylooks - as xlooks, but for Y axis
	extents - TODO
	resolution - [xres, yres] or None. Sets resolution output Ifg metadata. Use
	             None if raster size is not being changed.
	thresh - see thresh in prepare_ifgs().
	verbose - True = print gdalwarp output to stdout
	'''

	# dynamically build command for call to gdalwarp
	cmd = ["gdalwarp", "-overwrite", "-srcnodata", "None", "-te"] + extents
	if not verbose: cmd.append("-q")

	# HACK: if resampling, cut segment with gdalwarp & manually average tiles
	data = None
	if resolution:
		data = _resample_ifg(ifg, cmd, x_looks, y_looks, thresh)
		cmd += ["-tr"] + [str(r) for r in resolution] # change res of final output

	# use GDAL to cut (and resample) the final output layers
	s = splitext(ifg.data_path)
	ext = _file_ext(ifg)
	looks_path = s[0] + "_%srlks.%s" % (y_looks, ext)
	cmd += [ifg.data_path, looks_path]
	check_call(cmd)

	# Add missing/updated metadata to resampled ifg/DEM
	new_lyr = type(ifg)(looks_path, ifg.hdr_path)
	new_lyr.open(readonly=False)
	new_hdr_path = filename_pair(looks_path)[1]
	_create_new_roipac_header(ifg, new_lyr)

	# for non-DEMs, phase bands need extra metadata & conversions
	if hasattr(new_lyr, "phase_band"):
		new_lyr.phase_band.SetNoDataValue(nan)

		if data is None: # data wasn't resampled, so flag incoherent cells
			data = new_lyr.phase_band.ReadAsArray()
			data = where(data == 0, nan, data)

		# TODO: projection
		#if params.has_key(PROJECTION_FLAG):
		#	reproject()

		# tricky: write either resampled or the basic cropped data to new layer
		new_lyr.phase_band.WriteArray(data)
		new_lyr.nan_converted = True

	return new_lyr


def resample(data, xscale, yscale, threshold):
	"""Resamples/averages 'data' to return an array from the averaging of blocks
	of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

	data - source array to resample to different size
	xscale - number of cells to average along X axis
	yscale - number of Y axis cells to average
	threshold - minimum allowable proportion of NaN cells (range from 0.0-1.0),
	eg. 0.25 = 1/4 or more as NaNs results in a NaN value for the output cell.
	"""
	if threshold < 0 or threshold > 1:
		raise ValueError("threshold must be >= 0 and <= 1")

	# TODO: check scaling factors are ints

	ysize, xsize = data.shape
	xres, yres = (xsize / xscale), (ysize / yscale)
	dest = zeros((yres, xres), dtype=float32) * nan
	tile_cell_count = xscale * yscale

	# calc mean without nans (fractional threshold ignores tiles with excess NaNs)
	for y,x in product(xrange(yres), xrange(xres)):
		tile = data[y * yscale : (y+1) * yscale, x * xscale : (x+1) * xscale]
		non_nans = [ tile[crd] for crd in product(xrange(yscale), xrange(xscale))
		                  if not isnan(tile[crd])]
		nan_fraction = (tile_cell_count - len(non_nans)) / float(tile_cell_count)

		if nan_fraction < threshold or (nan_fraction == 0 and threshold == 0):
			dest[y,x] = mean(non_nans)

	return dest


def reproject():
	raise NotImplementedError("TODO: Reprojection LOS/Horiz/Vert")


def _create_new_roipac_header(src_ifg, new_ifg, dest=None):
	"""Translates an existing ROIPAC Ifg header to new values following cropping
	and resampling."""

	geotrans = new_ifg.dataset.GetGeoTransform()
	newpars = { WIDTH: new_ifg.dataset.RasterXSize,
							FILE_LENGTH: new_ifg.dataset.RasterYSize,
							X_FIRST: geotrans[0],
							X_STEP: geotrans[1],
							Y_FIRST: geotrans[3],
							Y_STEP: geotrans[5] }

	if isinstance(src_ifg, Ifg):
		newpars[WAVELENGTH] = src_ifg.WAVELENGTH
	elif isinstance(src_ifg, DEM):
		for a in [Z_OFFSET, Z_SCALE, PROJECTION, DATUM]:
			newpars[a] = getattr(src_ifg, a)

	if dest is None:
		dest = filename_pair(new_ifg.data_path)[1]
	write_roipac_header(newpars, dest)


def check_resolution(ifgs):
	"""Verifies Ifg resolutions are equal for the given grids"""
	for var in [X_STEP, Y_STEP]:
		values = array([getattr(i, var) for i in ifgs])
		if not (values == values[0]).all():
			msg = "Grid resolution does not match for %s" % var
			raise PreprocessingException(msg)


def check_looks(params):
	"""Verifies looks parameters are valid"""
	xscale = params[IFG_LKSX]
	yscale = params[IFG_LKSY]
	if not (isinstance(xscale, Number) and isinstance(yscale, Number)):
		msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xscale, yscale)
		raise PreprocessingException(msg)

	if not (xscale > 0 and yscale > 0):
		msg = "Invalid looks parameter(s), x: %s, y: %s" % (xscale, yscale)
		raise PreprocessingException(msg)


def min_bounds(ifgs):
	'''Returns bounds for overlapping area of the given interferograms.'''
	xmin = max([i.X_FIRST for i in ifgs])
	ymax = min([i.Y_FIRST for i in ifgs])
	xmax = min([i.X_LAST for i in ifgs])
	ymin = max([i.Y_LAST for i in ifgs])
	return xmin, ymin, xmax, ymax


def max_bounds(ifgs):
	'''Returns bounds for the total area covered by the given interferograms.'''
	xmin = min([i.X_FIRST for i in ifgs])
	ymax = max([i.Y_FIRST for i in ifgs])
	xmax = max([i.X_LAST for i in ifgs])
	ymin = min([i.Y_LAST for i in ifgs])
	return xmin, ymin, xmax, ymax


def check_crop_coords(ifgs, xmin, xmax, ymin, ymax, use_exceptions=False):
	'''Ensures cropping coords line up with grid system within tolerance.'''
	# NB: assumption is the first Ifg is correct, so only test against it
	i = ifgs[0]
	for par, crop, step in zip([X_FIRST, X_LAST, Y_FIRST, Y_LAST],
														[xmin, xmax, ymax, ymin],
														[i.X_STEP, i.X_STEP, i.Y_STEP, i.Y_STEP]):

		# is diff of the given extent from grid a multiple of X|Y_STEP ?
		param = getattr(i, par)
		diff = abs(crop - param)
		remainder = abs(modf(diff / step)[0])

		# handle cases where division gives remainder near zero, or just < 1
		if remainder > GRID_TOL and remainder < (1 - GRID_TOL):
			msg = "%s crop extent not within %s of grid coordinate" % (par, GRID_TOL)
			if use_exceptions:
				raise PreprocessingException(msg)
			sys.stderr.write("WARN: %s\n" % msg)


class PreprocessingException(Exception):
	pass
