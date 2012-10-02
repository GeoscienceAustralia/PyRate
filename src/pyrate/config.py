'''
Utilities to parse pyrate.conf and PyRate general config files

Created on 17/09/2012
@author: bpd900
'''

from numpy import unique, histogram 

from ifgconstants import DATE, DATE12


# constants for lookups
NUMBER_OF_SETS = 'nsets'
IFG_FILE_LIST = 'ifgfilelist'
OBS_DIR = 'obsdir'
OUT_DIR = 'outdir'
NUM_SETS = 'nsets'
SIM_DIR = 'simdir'
AMPLITUDE_FLAG = 'ampflag'
PERP_BASELINE_FLAG = 'basepflag'


# TODO

# Lookup to help convert args to correct type/defaults
# format is    key : (conversion, default value)
# None = no conversion
PARAM_CONVERSION = { OBS_DIR : (None, "obs"),
											IFG_FILE_LIST : (None, "ifg.list"),
											OUT_DIR : (None, "out"),
											PERP_BASELINE_FLAG : (bool, True),
											AMPLITUDE_FLAG : (bool, False),
											NUM_SETS : (int, 1),
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
	
	def __init__(self, date=None, repeat=None, span=None):
		self.date = date
		self.repeat = repeat
		self.span = span


def get_epochs(ifgs):
	masters = [i.header[DATE] for i in ifgs]
	slaves = [i.header[DATE12][-1] for i in ifgs]
		
	combined = masters + slaves
	#date_set, _, n = unique(combined)  # TODO: need newer numpy
	date_set = unique(combined)
	
	hist = None # TODO: histogram([d.toordinal() for d in date_set], bins=len(date_set))
	span = None
	
	return EpochList(date_set, hist, span)
