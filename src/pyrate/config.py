'''
Utilities to parse pyrate.conf and PyRate general config files

Created on 17/09/2012
@author: bpd900
'''

import numpy

from ifgconstants import X_FIRST, Y_FIRST, WIDTH, FILE_LENGTH, X_STEP, Y_STEP


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

IFG_CROP_OPT = 'ifgcropopt' # 1: minimum, 2: maximum, 3: customize, 4: all ifms already same size
IFG_LKSX = 'ifglksx'
IFG_LKSY = 'ifglksy'

IFG_XFIRST = 'ifgxfirst'
IFG_XLAST = 'ifgxlast'
IFG_YFIRST = 'ifgyfirst'
IFG_YLAST = 'ifgylast'

REFX = 'refx'
REFY = 'refy'
REFNX = "refnx"
REFNY = "refny"
REF_CHIP_SIZE = 'refchipsize'
REF_MIN_FRAC = 'refminfrac'



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
					IFG_LKSX : (int, 0),
					IFG_LKSY : (int, 0),
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
