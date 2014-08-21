'''
Main workflow script for PyRate

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import os
import logging
import datetime

import config as cf
from shared import Ifg

import mst
import prepifg
import algorithm
import refpixel
import orbital


# constants for metadata flags
META_UNITS = 'PHASE_UNITS'
MILLIMETRES = 'MILLIMETRES'
META_ORBITAL = 'ORBITAL_ERROR'
META_REMOVED = 'REMOVED'


# TODO: does verbose need to be included? Simpler to dump everything into a log
# TODO: add clean exception handling
def main(cfgfile='pyrate.conf', verbose=True):
	"""TODO: pirate workflow"""

	# TODO: add parameter error checking to fail fast, before number crunching
	try:
		params = cf.get_config_params(cfgfile)
	except IOError as err:
		msg = 'Config file error: %s "%s"' % (err.strerror, err.filename)
		logging.debug(msg)
		print msg
		return err.errno

	# NB: keep source files intact, should be run after prepifg code
	ifglist = cf.parse_namelist(params[cf.IFG_FILE_LIST])

	if warp_required(params[cf.IFG_LKSX],
						params[cf.IFG_LKSY],
						params[cf.IFG_CROP_OPT]):

		# TODO check for mlooked files
		mlooked_exists = None

		if mlooked_exists:
			logging.debug('Loading transformed interferograms...')
			raise NotImplementedError
		else:
			# TODO: get list of mlooked/cropped ifgs if they exist
			# otherwise warp raw ifgs, and return these as 'ifgs' (sans DEM)
			logging.debug('Transforming interferograms...')
			raise NotImplementedError
	else:
		# TODO: edits the ifgs in place, not the 'out' ones
		ifg_paths = [os.path.join(params[cf.OBS_DIR], p) for p in ifglist]
		ifgs = [Ifg(p) for p in ifg_paths]

	process_ifgs(ifgs, params)


def process_ifgs(ifgs, params):
	'''
	High level function to perform correction steps on supplied ifgs
	ifgs: sequence of Ifg objs (unopened)
	params: dict of run config params
	'''
	for i in ifgs:
		i.open(readonly=False)
		convert_wavelength(i)

	remove_orbital_error(ifgs, params)

	mst_grid = mst.mst_matrix_ifgs_only(ifgs)
	refpx, refpy = find_reference_pixel(ifgs, params)

	# TODO: missing
	# VCM
	# Linear Rate
	# Time series

	# final cleanup
	while ifgs:
		i = ifgs.pop()
		i.write_phase()
		i.dataset.FlushCache()
		i = None # force close    TODO: may need to implement close()

	logging.debug('End PyRate processing\n')


def warp_required(xlooks, ylooks, crop):
	"""Returns True if params show rasters need to be cropped and/or resized."""
	if xlooks > 1 or ylooks > 1:
		return True

	if crop is None or crop == prepifg.ALREADY_SAME_SIZE:
		return False

	return True


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
		msg = '%s: ignored as previous wavelength conversion detected'
		logging.debug(msg % ifg.data_path)
		return

	ifg.data = algorithm.wavelength_radians_to_mm(ifg.phase_data, ifg.wavelength)
	ifg.dataset.SetMetadataItem(META_UNITS, MILLIMETRES)
	msg = '%s: converted wavelength to millimetres'
	logging.debug(msg % ifg.data_path)


def remove_orbital_error(ifgs, params):
	if not params[cf.ORBITAL_FIT]:
		logging.debug('Orbital correction not required.')
		return

	# perform some general error/sanity checks
	flags = [i.dataset.GetMetadataItem(META_ORBITAL) for i in ifgs]

	if all(flags):
		msg = 'Skipped orbital correction, ifgs already corrected'
		logging.debug(msg)
		return
	else:
		check_orbital_ifgs(ifgs, flags)

	mlooked = None
	if params[cf.ORBITAL_FIT_LOOKS_X] > 1 or params[cf.ORBITAL_FIT_LOOKS_Y] > 1:
		# resampling here to use all prior corrections to orig data
		# TODO: avoid writing mlooked to disk by using mock ifgs/in mem arrays?
		mlooked = prepifg.prepare_ifgs(ifgs,
									crop_opt=prepifg.ALREADY_SAME_SIZE,
									xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
									ylooks=params[cf.ORBITAL_FIT_LOOKS_Y])

	orbital.orbital_correction(ifgs,
							degree=params[cf.ORBITAL_FIT_DEGREE],
							method=params[cf.ORBITAL_FIT_METHOD],
							mlooked=mlooked)

	# mlooked layers discarded as not used elsewhere
	if mlooked:
		for path in [m.data_path for m in mlooked]:
			msg = '%s: deleted (multilooked orbital correction file)'
			logging.debug(msg % path)
			os.remove(path)

	for i in ifgs:
		i.dataset.SetMetadataItem(META_ORBITAL, META_REMOVED)
		logging.debug('%s: orbital error removed' % i.data_path)


def check_orbital_ifgs(ifgs, flags):
	count = sum([f == META_REMOVED for f in flags])
	if count < len(flags) and count > 0:
		msg = 'Detected mix of corrected and uncorrected orbital error in ifgs'
		logging.debug(msg)

		for i, flag in zip(ifgs, flags):
			if flag:
				msg = '%s: prior orbital error correction detected'
			else:
				msg = '%s: no orbital correction detected'
			logging.debug(msg % i.data_path)

		raise orbital.OrbitalError(msg)


def find_reference_pixel(ifgs, params):
	# unlikely, but possible the refpixel can be (0,0)
	# check if there is a pre-specified reference pixel coord
	refx = params.get(params[cf.REFX])
	if refx > ifgs[0].ncols - 1:
		raise ValueError("Invalid reference pixel X coordinate: %s" % refx)

	refy = params.get(params[cf.REFY])
	if refy > ifgs[0].nrows - 1:
		raise ValueError("Invalid reference pixel Y coordinate: %s" % refy)

	if refx >= 0 and refy >= 0:
		msg = 'Reusing config file reference pixel (%s, %s)'
		logging.debug(msg % (refx, refy))
		return (refy, refx) # reuse preset ref pixel

	refy, refx = refpixel.ref_pixel(ifgs,
									params[cf.REFNX],
									params[cf.REFNY],
									params[cf.REF_CHIP_SIZE],
									params[cf.REF_MIN_FRAC])

	logging.debug('Reference pixel coordinate: (%s, %s)' % (refx, refy))
	return refx, refy


# general function template
#
# add check for pre-existing metadata flag / skip if required
# perform calculation
# optionally save modified data to disk if required
# optionally save correction component to disk (more useful for debugging)
# set flag in dataset for correction
# write to log file


if __name__ == "__main__":
	from optparse import OptionParser
	parser = OptionParser()
	# TODO: add options as they arise
	options, args = parser.parse_args()

	init_logging(logging.DEBUG)

	if args:
		if len(args) != 1:
			parser.error('Too many args')
		main(cfgfile=args[0])
	else:
		main()
