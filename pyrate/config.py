'''
Utilities to parse pyrate.conf config files. Includes numerous general PyRate
constants relating to options in config files.

Created on 17/09/2012
@author: bpd900
'''

# TODO: add regex column to check if some values are within bounds? Potential
# problem with the checking being done in the middle of the runs, as bad values
# could cause crashes & destroying some of the results.

from orbital import NETWORK_METHOD, QUADRATIC


# general constants
NO_MULTILOOKING = 1

# constants for lookups
NUMBER_OF_SETS = 'nsets'
IFG_FILE_LIST = 'ifgfilelist'
OBS_DIR = 'obsdir'
OUT_DIR = 'outdir'
NUM_SETS = 'nsets'
SIM_DIR = 'simdir'
DEM_FILE = 'demfile'
AMPLITUDE_FLAG = 'ampflag'
PERP_BASELINE_FLAG = 'basepflag'
PROJECTION_FLAG = 'prjflag'

IFG_CROP_OPT = 'ifgcropopt' # 1=min, 2=max, 3=custom, 4=all ifgs already same size
IFG_LKSX = 'ifglksx'
IFG_LKSY = 'ifglksy'

IFG_XFIRST = 'ifgxfirst'
IFG_XLAST = 'ifgxlast'
IFG_YFIRST = 'ifgyfirst'
IFG_YLAST = 'ifgylast'

# reference pixel parameters
REFX = 'refx'
REFY = 'refy'
REFNX = "refnx"
REFNY = "refny"
REF_CHIP_SIZE = 'refchipsize'
REF_MIN_FRAC = 'refminfrac'

# orbital error correction/parameters
ORBITAL_FIT = 'orbfit' # numerical BOOL (1/0), do/don't do orbital correction
ORBITAL_FIT_METHOD = 'orbfitmethod'  # BOOL (1/2) 1: ifg by ifg/independent, 2: epoch by epoch
ORBITAL_FIT_DEGREE = 'orbfitdegrees' # BOOL (1/2) 1=planar, 2=quadratic
ORBITAL_FIT_LOOKS_X = 'orbfitlksx' # int of 1+, X multi looking factor
ORBITAL_FIT_LOOKS_Y = 'orbfitlksy' # int of 1+, Y multi looking factor
# ORBITAL_FIT_orbrefest:     1 BOOLEAN (1/0) # remove reference phase
# ORBITAL_FIT_ orbmaskflag:   1 BOOLEAN (1/0) # mask some patches for orbital correction


# Lookup to help convert args to correct type/defaults
# format is    key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = { OBS_DIR : (None, "obs"),
					IFG_FILE_LIST : (None, "ifg.list"),
					OUT_DIR : (None, "out"),
					DEM_FILE : (str, None),
					PERP_BASELINE_FLAG : (bool, True),
					AMPLITUDE_FLAG : (bool, False),
					NUM_SETS : (int, 1),
					IFG_CROP_OPT : (int, None),
					IFG_LKSX : (int, NO_MULTILOOKING),
					IFG_LKSY : (int, NO_MULTILOOKING),
					IFG_XFIRST : (float, None),
					IFG_XLAST : (float, None),
					IFG_YFIRST : (float, None),
					IFG_YLAST : (float, None),
					PROJECTION_FLAG : (int, 3),

					REFX : (int, 0),
					REFY : (int, 0),
					REFNX : (int, None), # was 50 in original Pirate code
					REFNY : (int, None), # was 50 in original Pirate code
					REF_CHIP_SIZE : (int, None), # defaults to 21 in orig
					REF_MIN_FRAC : (float, 0.8), # uses Pirate default

					ORBITAL_FIT : (bool, True),
					ORBITAL_FIT_METHOD : (int, NETWORK_METHOD),
					ORBITAL_FIT_DEGREE : (int, QUADRATIC),
					ORBITAL_FIT_LOOKS_X : (int, NO_MULTILOOKING),
					ORBITAL_FIT_LOOKS_Y : (int, NO_MULTILOOKING)
				}


def parse_conf_file(conf_file):
	"""Returns a dict for the key:value pairs from the .conf file"""
	with open(conf_file) as f:
		txt = f.read().splitlines()
		lines = [line.split() for line in txt if line != "" and line[0] not in "%#"]
		lines = [(e[0].rstrip(":"), e[1]) for e in lines] # strip colons from keys
		parameters = dict(lines)
		_parse_pars(parameters)
		return parameters


def _parse_pars(pars):
	"""Parses and converts config file params from text"""
	for k in PARAM_CONVERSION.keys():
		if pars.has_key(k):
			conversion_func = PARAM_CONVERSION[k][0]
			if conversion_func:
				pars[k] = conversion_func(pars[k])
		else:
			# revert empty options to default value
			pars[k] = PARAM_CONVERSION[k][1]

	return pars


def parse_namelist(nml):
	"""Parses name list file into array of paths"""
	with open(nml) as f:
		return [ln.strip() for ln in f.readlines() if ln != ""]


class ConfigException(Exception):
	'''Default exception class for configuration errors.'''
	pass
