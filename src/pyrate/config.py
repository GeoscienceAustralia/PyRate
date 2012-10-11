'''
Utilities to parse pyrate.conf and PyRate general config files

Created on 17/09/2012
@author: bpd900
'''

import numpy
from numpy import unique, reshape, histogram

from shared import Ifg, IfgException


# constants for lookups
NUMBER_OF_SETS = 'nsets'
IFG_FILE_LIST = 'ifgfilelist'
OBS_DIR = 'obsdir'
OUT_DIR = 'outdir'
NUM_SETS = 'nsets'
SIM_DIR = 'simdir'
AMPLITUDE_FLAG = 'ampflag'
PERP_BASELINE_FLAG = 'basepflag'

IFG_CROP_OPT = 'ifgcropopt' # 1: minimum, 2: maximum, 3: customize, 4: all ifms already same size
IFG_LKSX = 'ifglksx'  # INT
IFG_LKSY = 'ifglksy' # INT


# Lookup to help convert args to correct type/defaults
# format is    key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = { OBS_DIR : (None, "obs"),
					IFG_FILE_LIST : (None, "ifg.list"),
					OUT_DIR : (None, "out"),
					PERP_BASELINE_FLAG : (bool, True),
					AMPLITUDE_FLAG : (bool, False),
					NUM_SETS : (int, 1),
					IFG_CROP_OPT : (int, None),
					IFG_LKSX : (int, 0),
					IFG_LKSY : (int, 0),
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



class EpochList(object):
	'''TODO'''

	def __init__(self, date=None, repeat=None, span=None):
		self.date = date
		self.repeat = repeat
		self.span = span


def get_epochs(ifgs):
	masters = [i.MASTER for i in ifgs]
	slaves = [i.SLAVE for i in ifgs]

	combined = masters + slaves
	dates, n = unique(combined, False, True)
	repeat, _ = histogram(n, bins=len(set(n)))

	# absolute span for each date from the zero/start point
	span = [ (dates[i] - dates[0]).days / 365.25 for i in range(len(dates)) ]
	return EpochList(dates, repeat, span)


def prepare_ifgs(ifgs, params, conversion=None, amplitude=None, projection=None):

	# TODO: initial port of the ugly prepifg.m code
	res = _check_xy_steps(ifgs)
	if res is False:
		msg = "Cell sizes unequal for supplied interferograms"
		raise IfgException(msg)

	# TODO: handle multilooking
	# TODO: does the multilooking need to be in squares? ie. SX==SY?
	if params[IFG_LKSX] > 0 and params[IFG_LKSY] > 0:
		raise NotImplementedError
	else:
		raise NotImplementedError


def _check_xy_steps(ifgs):
	'''Validates X_STEP and Y_STEP for given list of interferograms. Returns True
	if the values for X_STEP match (and for Y_STEP)'''
	xsteps = numpy.array([i.X_STEP for i in ifgs])
	if not numpy.all(xsteps == xsteps[0]):
		return False

	ysteps = numpy.array([i.Y_STEP for i in ifgs])
	return numpy.all(ysteps == ysteps[0])






