'''
Main workflow script for PyRate

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import config

from shared import Ifg
import algorithm


# metadata constants
META_UNITS = 'PHASE_UNITS'
MILLIMETRES = 'MILLIMETRES'


# TODO: add basic logging statements 
def main(cfgfile='pyrate.conf', verbose=True):
	"""TODO: pirate workflow"""
	raise NotImplementedError

	# TODO: add parameter error checking to fail fast, before number crunching
	#params = config.get_config_params(cfgfile)

	# TODO: get list of mlooked/cropped ifgs
	# NB: keep source files intact, should be run after prepifg code
	#ifglist = config.parse_namelist(params[config.IFG_FILE_LIST])
	#ifg_namelist = [os.path.join(params[config.OBS_DIR], p) for p in ifglist]
	#ifgs = [Ifg(p) for p in ifg_namelist]


def process_ifgs(ifgs, params):
	'''
	Perform correction steps on the given ifgs
	ifgs: sequence of Ifg objs (unopened)
	params: dict of run config params
	'''
	for i in ifgs:
		i.open()
		convert_wavelength(i)

	# final close
	while ifgs:
		i = ifgs.pop()
		i.dataset.FlushCache()
		i = None # force close TODO: may need to implement close()


def convert_wavelength(ifg):
	ifg.data = algorithm.wavelength_radians_to_mm(ifg.phase_data, ifg.wavelength)
	ifg.dataset.SetMetadataItem(META_UNITS, MILLIMETRES)
