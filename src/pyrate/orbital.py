'''
Module for calculating orbital correction for interferograms.

Created on 31/3/13
@author: Ben Davies
'''

from itertools import product
from numpy import sum, isnan, reshape, zeros, float32, vstack, squeeze
from scipy.linalg import lstsq
from numpy.linalg import pinv

import algorithm
from algorithm import ifg_date_lookup
from mst import default_mst

# Orbital correction tasks
#
# 4) forward calculation (orbfwd.m)
#		create 2D orbital correctionl layer from the model params


# TODO: options for multilooking
# 1) do the 2nd stage mlook at prepifg.py/generate up front, then delete in workflow after
# 2) refactor prep_ifgs() call to take input filenames and params & generate mlooked versions from that
#    this needs to be more generic to call at any point in the runtime.


# constants
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2

PLANAR = 1
QUADRATIC = 2


# TODO: do MST as a filter step outside the main func? More MODULAR?
# FIXME: operate on files

def orbital_correction(ifgs, degree, method, mlooked=None, offset=True):
	'''
	TODO: Top level method for correcting orbital error in the given Ifgs. It is
	assumed the given ifgs have been reduced using MST for the network method.

	ifgs - list of Ifg objs to correct
	degree - PLANAR or QUADRATIC
	method - INDEPENDENT_METHOD or NETWORK_METHOD
	mlooked - sequence of multilooked Ifgs
	offset = True/False to include the constant/offset component
	'''
	# TODO: save corrected layers to new file or use intermediate arrays?

	if degree not in [PLANAR, QUADRATIC]:
		msg = "Invalid degree of %s for orbital correction" % degree
		raise OrbitalError(msg)

	if method not in [INDEPENDENT_METHOD, NETWORK_METHOD]:
		msg = "Unknown method '%s', need INDEPENDENT or NETWORK method" % method
		raise OrbitalError(msg)

	if method == NETWORK_METHOD:
		# FIXME: net correction has to have some kind of handling for mlooked or not

		if mlooked:
			_validate_mlooked(mlooked, ifgs)
			return _get_net_correction(sub_mlooked, degree, offset) # TODO: fwd corr
		else:
			return _get_net_correction(ifgs, degree, offset) # TODO: fwd corr

	elif method == INDEPENDENT_METHOD:
		# FIXME: determine how to work this into the ifgs. Generate new Ifgs? Update
		# the old ifgs and flag the corrections in metadata?
		#
		#for i in ifgs:
		#	i.phase_data -= _get_ind_correction(i, degree, offset)
		return [_get_ind_correction(i, degree, offset) for i in ifgs]


def _validate_mlooked(mlooked, ifgs):
	'''Basic sanity checking of the multilooked ifgs.'''

	if len(mlooked) != len(ifgs):
		msg = "Mismatching # ifgs and # multilooked ifgs"
		raise OrbitalError(msg)

	tmp = [hasattr(i, 'phase_data') for i in mlooked]
	if all(tmp) is False:
		msg = "Mismatching types in multilooked ifgs arg:\n%s" % mlooked
		raise OrbitalError(msg)


def _get_ind_correction(ifg, degree, offset):
	'''Calculates and returns orbital correction array for an ifg'''

	vphase = reshape(ifg.phase_data, ifg.num_cells) # vectorised, contains NODATA
	dm = get_design_matrix(ifg, degree, offset)
	assert len(vphase) == len(dm)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model, _, rank, _ = lstsq(tmp, fd)
	exp_len = 2 if degree == PLANAR else 5
	if offset:
		exp_len += 1
	assert len(model) == exp_len

	# calculate forward model & morph back to 2D
	tmp = sum(dm * model, axis=1) # d = ax + by
	correction = reshape(tmp, ifg.phase_data.shape)
	return correction


def _get_net_correction(ifgs, degree, offset):
	'''
	Returns the TODO
	ifgs - assumed to be Ifgs from a prior MST step
	degree - PLANAR or QUADRATIC
	offset - True/False for including TODO
	'''

	# FIXME: correction needs to apply to *all* ifgs - can still get the correction
	#        as there are model params for each epoch (+ ifgs relate to all those epochs).
	# TODO: multilooking (do seperately/prior to this as a batch job?)

	# get DM / clear out the NaNs based on obs
	tmp = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in ifgs])
	vphase = squeeze(tmp)
	dm = get_network_design_matrix(ifgs, degree, offset)
	assert len(vphase) == len(dm)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model = pinv(tmp, 1e-6) * fd

	# TODO forward correction

	return model


def get_design_matrix(ifg, degree, offset):
	'''Returns design matrix with 2 columns for linear model parameters'''

	nparams = 2 if degree == PLANAR else 5
	if offset:
		nparams += 1  # eg. y = mx + offset

	# init design matrix
	shape = (ifg.num_cells, nparams)
	data = zeros(shape, dtype=float32)
	rows = iter(data)

	dmfun = _planar_dm if degree == PLANAR else _quadratic_dm
	dmfun(ifg, rows, offset)
	return data


# TODO: can this be refactored under one get_design_matrix() func?
def get_network_design_matrix(ifgs, degree, offset):
	'''Returns a larger format design matrix for networked error correction.'''

	if degree not in [PLANAR, QUADRATIC]:
		raise OrbitalError("Invalid degree argument")

	num_ifgs = len(ifgs)
	if num_ifgs < 1:
		# can feasibly do correction on a single Ifg/2 epochs
		raise OrbitalError("Invalid number of Ifgs")

	# TODO: refactor to prevent duplication with get_design_matrix()?
	nparams = 2 if degree == PLANAR else 5

	# sort out master and slave date IDs
	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	ids = algorithm.master_slave_ids(dates)
	num_epochs = max(ids.values()) + 1 # convert from zero indexed ID
	np = nparams * num_epochs

	# init design matrix
	shape = [ifgs[0].num_cells * num_ifgs, np]
	if offset:
		shape[1] += num_ifgs # add extra cols for offset component
	data = zeros(shape, dtype=float32)

	# paste in individual design matrices
	for i, ifg in enumerate(ifgs):
		tmp = get_design_matrix(ifg, degree, False) # network method doesn't have extra col
		rs = i * ifg.num_cells
		rf = rs + ifg.num_cells

		# generate column indices into data based on master position
		mascs = ids[ifg.MASTER] * nparams
		mascf = mascs + nparams
		data[rs:rf, mascs:mascf] = -tmp

		# then for slave
		slvcs =	ids[ifg.SLAVE] * nparams
		slvcf = slvcs + nparams
		data[rs:rf, slvcs:slvcf] = tmp

		if offset:
			data[rs:rf, np + i] = 1  # add offset cols

	return data


def _planar_dm(ifg, rows, offset):
	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)

	# TODO: optimise with meshgrid calls?
	# TODO: make more efficient by pre generating xranges and doing array ops?
	# TODO: coordinates generator for Ifgs?

	if offset:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next() # TODO: make faster with vstack?
			row[:] = [x * ifg.X_SIZE, y * ifg.Y_SIZE, 1]
	else:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next() # TODO: make faster with vstack?
			row[:] = [x * ifg.X_SIZE, y * ifg.Y_SIZE]


def _quadratic_dm(ifg, rows, offset):
	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)
	yst, xst = ifg.Y_SIZE, ifg.X_SIZE

	# TODO: refactor, use ones +/- final col and paste these values over it
	if offset:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next()
			y2 = y * yst
			x2 = x * xst
			row[:] = [x2**2, y2**2, x2*y2, x2, y2, 1]
	else:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next()
			y2 = y * yst
			x2 = x * xst
			row[:] = [x2**2, y2**2, x2*y2, x2, y2]



class OrbitalError(Exception):
	'''Generic class for errors in orbital correction'''
	pass
