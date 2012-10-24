"""
Created on 23/10/2012
@author: Ben Davies
"""

import os
from os.path import join, splitext
from subprocess import check_call

from shared import Ifg
from config import OBS_DIR, IFG_CROP_OPT, IFG_LKSX, IFG_LKSY, IFG_FILE_LIST
from config import IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST
from ifgconstants import X_FIRST, Y_FIRST

# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
NO_CROP = 4



def prepare_ifgs(params):
	with open(params[IFG_FILE_LIST]) as f:
		filelist = f.readlines()
	
	paths = [join(params[OBS_DIR], p) for p in filelist ]
	ifgs = [Ifg(p) for p in paths]
	for i in ifgs:
		i.open() # creates header files if needed
		
	# TODO: for CROP, LOOKING, RESCALE, use gdalwarp  eg:  gdalwarp -te 1 -2 5 2 -tr 2 2 -r bilinear
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
		# TODO: warn if the extents don't match up with the grid
		xmin = params[IFG_XFIRST]
		ymax = params[IFG_YFIRST]
		xmax = params[IFG_XLAST]
		ymin = params[IFG_YLAST]
		
	elif crop_opt == NO_CROP:
		raise NotImplementedError("Same as maximum for IFGs of same size")
		
	else:
		raise PreprocessingException("Uncrecognised crop option: %s" % crop_opt)
	
	# use GDAL to reproject dataset
	extents = [str(s) for s in (xmin, ymin, xmax, ymax) ]
	resolution = None
	
	for i in ifgs:
		s = splitext(i.data_path)
		looks_path = s[0] + "_TODO" + ".tif" 
		args = ["gdalwarp", "-te"] + extents + [i.data_path, looks_path]
		check_call(args)	


class PreprocessingException(Exception):
	pass
