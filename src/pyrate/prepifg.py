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

from shared import Ifg
from roipac import filename_pair, write_roipac_header
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST
from config import IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST
from ifgconstants import X_FIRST, Y_FIRST, X_LAST, Y_LAST, X_STEP, Y_STEP
from ifgconstants import WIDTH, FILE_LENGTH, WAVELENGTH

# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
NO_CROP = 4
GRID_TOL = 1e-6



def prepare_ifgs(params, use_exceptions=False, verbose=False):
	"""Produces multilooked/resampled data files for PyRate analysis"""
	
	check_looks(params)
	with open(params[IFG_FILE_LIST]) as f:
		filelist = f.readlines()
	
	paths = [join(params[OBS_DIR], p) for p in filelist ]
	ifgs = [Ifg(p) for p in paths]
	check_resolution(ifgs)	
	for i in ifgs:
		i.open()
		
	# calculate crop option for passing to gdalwarp
	crop_opt = params[IFG_CROP_OPT]
	if crop_opt == MINIMUM_CROP:
		xmin = max([i.X_FIRST for i in ifgs])
		ymax = min([i.Y_FIRST for i in ifgs])
		xmax = min([i.X_LAST for i in ifgs])
		ymin = max([i.Y_LAST for i in ifgs])
		
	elif crop_opt == MAXIMUM_CROP:
		xmin = min([i.X_FIRST for i in ifgs])
		ymax = max([i.Y_FIRST for i in ifgs])
		xmax = max([i.X_LAST for i in ifgs])
		ymin = min([i.Y_LAST for i in ifgs])

	elif crop_opt == CUSTOM_CROP:
		xmin = params[IFG_XFIRST]
		ymax = params[IFG_YFIRST]
		xmax = params[IFG_XLAST]
		ymin = params[IFG_YLAST]
		
		# check cropping coords line up with grid system within tolerance
		# NB: only tests against the first Ifg
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
		
	elif crop_opt == NO_CROP:
		raise NotImplementedError("Same as maximum for IFGs of same size")
		
	else:
		raise PreprocessingException("Unrecognised crop option: %s" % crop_opt)
	
	# calculate args for gdalwarp reprojection
	extents = [str(s) for s in (xmin, ymin, xmax, ymax) ]	
	resolution = None
	if params[IFG_LKSX] > 1 or params[IFG_LKSY] > 1:
		resolution = [params[IFG_LKSX] * i.X_STEP, params[IFG_LKSY] * i.Y_STEP ]
	
	# generate args for gdalwarp for each interferogram
	for i in ifgs:
		s = splitext(i.data_path)
		looks_path = s[0] + "_%srlks.tif" % params[IFG_LKSY] # NB: hardcoded to .tif  
		cmd = ["gdalwarp", "-overwrite", "-srcnodata", "None", "-te"] + extents
		if not verbose: cmd.append("-q")		
				
		# HACK: if resampling, cut original seg with gdalwarp & perform tile averaging
		# resampling method (gdalwarp lacks Pirate averaging method)
		data = None
		if resolution:			
			tmp_path = mkstemp()[1]
			check_call(cmd + [i.data_path, tmp_path])
			ifg_tmp = Ifg(tmp_path, i.hdr_path)
			ifg_tmp.open()
			data = ifg_tmp.phase_band.ReadAsArray()
			data = where(data == 0, nan, data) # flag incoherent cells as NaN
			data = resample(data, params[IFG_LKSX], params[IFG_LKSY])
			del ifg_tmp
			os.remove(tmp_path)
									
			# args to change resolution for final outout
			cmd += ["-tr"] + [str(r) for r in resolution]
		
		cmd += [i.data_path, looks_path]
		check_call(cmd)		
		
		# Add missing metadata to new interferograms
		ifg = Ifg(looks_path, i.hdr_path)
		ifg.open(readonly=False)
		ifg.phase_band.SetNoDataValue(nan) # phase == 0 is incoherent
		
		new_hdr_path = filename_pair(looks_path)[1]
		_create_new_roipac_header(i, ifg)
		
		# data is only None if there is no resampling
		if data is None:
			data = ifg.phase_band.ReadAsArray()
			data = where(data == 0, nan, data) # flag incoherent cells as NaN
		
		ifg.phase_band.WriteArray(data)


def resample(data, xscale, yscale):
	"""Resamples/averages 'data' from tile size given by scaling factors. Assumes
	incoherent cells have been converted to NaNs."""
	
	ysize, xsize = data.shape
	xres, yres = (xsize / xscale), (ysize / yscale)
	dest = zeros((yres, xres), dtype=float32) * nan
	
	# alternate mean without nans
	# TODO: threshold for # NaN cells?
	for y,x in product(xrange(yres), xrange(xres)):
		tile = data[y * yscale : (y+1) * yscale, x * xscale : (x+1) * xscale]
		non_nans = [ tile[crd] for crd in product(xrange(yscale), xrange(xscale)) if not isnan(tile[crd])]
		
		if non_nans:
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
		Y_STEP: geotrans[5],
		WAVELENGTH: src_ifg.WAVELENGTH }
	
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


class PreprocessingException(Exception):
	pass
