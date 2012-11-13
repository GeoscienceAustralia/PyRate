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
from ifgconstants import X_FIRST, Y_FIRST, X_LAST, Y_LAST, X_STEP, Y_STEP
from ifgconstants import WIDTH, FILE_LENGTH, WAVELENGTH
from ifgconstants import Z_OFFSET, Z_SCALE, PROJECTION, DATUM

# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP]

GRID_TOL = 1e-6



def prepare_ifgs(params, threshold=0.5, use_exceptions=False, verbose=False):
	"""Produces multilooked/resampled data files for PyRate analysis. Params is a
	dict of named values from pyrate config file. Threshhold controls resampling
	to coarser grids, ranging from 0.0 -> 1.0. Threshold is the proportion at and
	above which the proportion of NaNs in a segment is considered invalid. Setting
	threshhold=0 resamples to NaN if 1+ contributing cells are NaNs. At 0.25, it
	resamples to NaN if 1/4 or more contributing cells are NaNs. At 1.0, segments
	are resampled to NaN only if all contributing cells are NaNs."""
	
	# validate config file settings
	crop_opt = params[IFG_CROP_OPT]
	if not crop_opt in CROP_OPTIONS:
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

		# check cropping coords line up with grid system within tolerance
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

	# calculate args for gdalwarp reprojection
	extents = [str(s) for s in (xmin, ymin, xmax, ymax) ]
	resolution = None
	if params[IFG_LKSX] > 1 or params[IFG_LKSY] > 1:
		resolution = [params[IFG_LKSX] * i.X_STEP, params[IFG_LKSY] * i.Y_STEP ]

	# Generate gdalwarp args for each interferogram. Some trickery is required to
	# interface with gdalwarp. For resampling, gdalwarp is called 2x, once to subset
	# the source data for averaging, the second to generate the final dataset with
	# correct extents/shape/cell count. Without resampling, gdalwarp is needed 1x.
	for i in ifgs:
		s = splitext(i.data_path)
		if isinstance(i, Ifg):
			ext = "tif"
		elif hasattr(i, DATUM):
			ext = "dem"
		else:
			raise NotImplementedError("TODO: additional raster types?")

		looks_path = s[0] + "_%srlks.%s" % (params[IFG_LKSY], ext)
		cmd = ["gdalwarp", "-overwrite", "-srcnodata", "None", "-te"] + extents
		if not verbose: cmd.append("-q")

		# HACK: if resampling, cut segment with gdalwarp & manually do tile averages
		# resampling (gdalwarp lacks Pirate averaging method)
		data = None
		if resolution:
			tmp_path = mkstemp()[1]
			check_call(cmd + [i.data_path, tmp_path])

			i_tmp = type(i)(tmp_path, i.hdr_path) # dynamically handle Ifgs & Rasters
			i_tmp.open()

			if isinstance(i, Ifg):
				# TODO: resample amplitude band too?
				data = i_tmp.phase_band.ReadAsArray()
				data = where(data == 0, nan, data) # flag incoherent cells as NaNs
			elif isinstance(i, DEM):
				data = i_tmp.height_band.ReadAsArray()
			else:
				raise NotImplementedError("TODO: other raster types to handle?")

			data = resample(data, params[IFG_LKSX], params[IFG_LKSY], threshold)
			del i_tmp
			os.remove(tmp_path)

			# args to change resolution for final outout
			cmd += ["-tr"] + [str(r) for r in resolution]

		# use GDAL to cut (and resample) the final output layers
		cmd += [i.data_path, looks_path]
		check_call(cmd)

		# Add missing/updated metadata to resampled interferograms/rasters
		new_lyr = type(i)(looks_path, i.hdr_path) # dynamically handle Ifgs/Rasters
		new_lyr.open(readonly=False)
		new_hdr_path = filename_pair(looks_path)[1]
		_create_new_roipac_header(i, new_lyr) # create header file of revised values

		# for non-DEMs, phase bands need extra metadata & conversions
		if hasattr(new_lyr, "phase_band"):
			new_lyr.phase_band.SetNoDataValue(nan)

			if data is None:
				# data not resampled, so data for conversion must be read from new dataset
				data = new_lyr.phase_band.ReadAsArray()
				data = where(data == 0, nan, data)

			# tricky: write either resampled or the basic cropped data to new layer
			new_lyr.phase_band.WriteArray(data)


# TODO: refactor out to another module?
def resample(data, xscale, yscale, threshold):
	"""Resamples/averages 'data' from tile size given by scaling factors. Assumes
	incoherent cells have been converted to NaNs. threshold is the minimum allowable
	proportion of NaN cells (range from 0-1), eg. 0.25 = 1/4 or more as NaNs results
	in a NaN value for the output cell."""

	if threshold < 0 or threshold > 1:
		raise ValueError("threshold must be >= 0 and <= 1")

	ysize, xsize = data.shape
	xres, yres = (xsize / xscale), (ysize / yscale)
	dest = zeros((yres, xres), dtype=float32) * nan
	tile_cell_count = xscale * yscale

	# calc mean without nans (fractional threshold ignores tiles with excess NaNs)
	for y,x in product(xrange(yres), xrange(xres)):
		tile = data[y * yscale : (y+1) * yscale, x * xscale : (x+1) * xscale]
		non_nans = [ tile[crd] for crd in product(xrange(yscale), xrange(xscale)) if not isnan(tile[crd])]
		nan_fraction = (tile_cell_count - len(non_nans)) / float(tile_cell_count)

		if nan_fraction < threshold or (nan_fraction == 0 and threshold == 0):
			dest[y,x] = mean(non_nans)

	return dest


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


class PreprocessingException(Exception):
	pass
