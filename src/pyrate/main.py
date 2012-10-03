'''
Created on 17/09/2012
@author: bpd900
'''

import os
import roipac, config
from shared import Ifg

from config import IFG_FILE_LIST, OBS_DIR


def main(cfgfile='pyrate.conf', verbose=True):
	"""TODO: pirate workflow"""

	if verbose: print 'Reading configuration file ...'
	params = config.parse_conf_file(cfgfile)
	#preprocess_parameters(params)

	# get interferogram list
	if verbose:
		print 'Making ifg and epoch namelist ...'
	ifg_namelist = config.parse_namelist(params[IFG_FILE_LIST])
	ifg_namelist = [os.path.join(params[OBS_DIR], p) for p in ifg_namelist]

	ifgs = [Ifg(p) for p in ifg_namelist]
	epochlist = config.get_epochs(ifgs)
