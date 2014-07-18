'''
Main workflow script for PyRate

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import config
import logging
import datetime

from shared import Ifg
import algorithm, mst, refpixel


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
	High level function to perform correction steps on supplied ifgs
	ifgs: sequence of Ifg objs (unopened)
	params: dict of run config params
	'''
	init_logging(logging.DEBUG)
	for i in ifgs:
		i.open()
		convert_wavelength(i)

	# TODO: epochs = algorithm.get_epochs(ifgs)
	# TODO: mst_grid = mst.mst_matrix(ifgs, epochs)
	# TODO: refy, refx = refpixel.ref_pixel(params, ifgs)

	# final close
	while ifgs:
		i = ifgs.pop()
		i.dataset.FlushCache()
		i = None # force close TODO: may need to implement close()

	logging.debug('End PyRate processing\n')

# TODO: write to alternate file if log exists
def init_logging(level):
	t = datetime.datetime.now()
	path = 'pyrate_%s_%02d_%02d.log' % (t.year, t.month, t.day)
	fmt = '%(asctime)s %(message)s'
	datefmt = '%d/%m/%Y %I:%M:%S %p'
	logging.basicConfig(filename=path, format=fmt, datefmt=datefmt, level=level)
	logging.debug('Log started')


def convert_wavelength(ifg):
	if ifg.dataset.GetMetadataItem(META_UNITS) == MILLIMETRES:
		msg = '%s: previous wavelength conversion detected'
		logging.debug(msg % ifg.data_path)
		return

	ifg.data = algorithm.wavelength_radians_to_mm(ifg.phase_data, ifg.wavelength)
	ifg.dataset.SetMetadataItem(META_UNITS, MILLIMETRES)
	msg = '%s: converted wavelength to millimetres'
	logging.debug(msg % ifg.data_path)


# function template
#
# add check for pre-existing metadata flag / skip if required
# perform calculation
# optionally save modified data to disk if required
# optionally save correction component to disk (more useful for debugging) 
# set flag in dataset for correction
# write to log file
