'''
TODO

Created on 31/3/13
@author: Ben Davies
'''

from itertools import product
from numpy import sum, where, nan, isnan, reshape, zeros, float32
from numpy.ma import masked_array
from scipy.linalg import lstsq

import algorithm


# Orbital correction
# 0) Config file stuff:
#		Add the orbital params
#		orbfit = 0: off, 1=do it
#		orbfitmethod = 1 or 2 (see conf file)  (AKA independent and networked methods)
#		orbrefest = orb_ref_est, 0=off, 1=on (IGNORE for the time being - something to do with removing the constant param 'c' from y = mx+c)
#		orbmaskflag = IGNORE for now. Used to mask out some patches, eg. if there is an earthquake signal somewhere.
#
# 1) resample/multilook the data to have less pixels (for comp efficiency)
# optional - vu can handle the larger arrays (new sys=?) so not mandatory (orbcorrect.m)
#
# 2) design matrix (orbdesign.m)
#
# 3) linear regression/inversion - get params (ala orbcorrect.m)
#
# 4) forward calculation (orbfwd.m)
#		create 2D orbital correctionl layer from the model params


# constants
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2

PLANAR = 1
QUADRATIC = 2


def orbital_correction(ifgs, degree, method, offset=True):
	'''TODO'''

	# TODO: save corrected layers to new file or use intermediate arrays?
	# TODO: offsets

	if degree not in [PLANAR, QUADRATIC]:
		msg = "Invalid degree of %s for orbital correction" % degree
		raise OrbitalCorrectionError(msg)

	if method == NETWORK_METHOD:
		raise NotImplementedError

	return [_get_correction(i, degree, offset) for i in ifgs]


def _get_correction(ifg, degree, offset):
	'''Calculates and returns orbital correction array for an ifg'''

	# precalculate matrix shapes
	orig = ifg.phase_data.shape
	vsh = orig[0] * orig[1] # full shape of vectorised form of data
	vphase = reshape(ifg.phase_data, vsh) # vectorised data, contains NODATA

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
	correction = reshape(tmp, orig)
	return correction


def get_design_matrix(ifg, degree, offset):
	'''Returns design matrix with 2 columns for linear model parameters'''

	nparams = 2 if degree == PLANAR else 5
	if offset:
		nparams += 1  # eg. y = mx + offset

	# init design matrix
	shape = ((ifg.WIDTH * ifg.FILE_LENGTH), nparams)
	data = zeros(shape, dtype=float32)
	rows = iter(data)

	dmfun = _planar_dm if degree == PLANAR else _quadratic_dm
	dmfun(ifg, rows, offset)
	return data


# TODO: can this be refactored under one get_design_matrix() func?
def get_network_design_matrix(ifgs, degree, offset):
	'''TODO'''

	if degree not in [PLANAR, QUADRATIC]:
		raise OrbitalCorrectionError("Invalid degree argument")

	num_ifgs = len(ifgs)
	if num_ifgs < 2:
		raise OrbitalCorrectionError("Invalid number of Ifgs")

	num_epochs = num_ifgs + 1

	# TODO: refactor to prevent duplication here
	nparams = 2 if degree == PLANAR else 5
	if offset:
		nparams += 1  # eg. b/offset in (y = mx + b) is an extra param

	# sort out master and slave date IDs
	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	ids = algorithm.master_slave_ids(dates)

	# init design matrix
	nrows, ncols = ifgs[0].FILE_LENGTH, ifgs[0].WIDTH
	shape = (ifg.num_cells * num_ifgs, nparams * num_epochs)
	data = zeros(shape, dtype=float32)

	# paste in individual design matrices
	for i, ifg in enumerate(ifgs):
		tmp = get_design_matrix(ifg, degree, offset)
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
			row[:] = [x * ifg.X_STEP, y * ifg.Y_STEP, 1] # FIXME: change to (Y|X)SIZE
	else:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next() # TODO: make faster with vstack?
			row[:] = [x * ifg.X_STEP, y * ifg.Y_STEP] # FIXME: change to (Y|X)SIZE, needs proj4


def _quadratic_dm(ifg, rows, offset):
	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)
	yst, xst = ifg.Y_STEP, ifg.X_STEP # FIXME: sizes

	if offset:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next()
			y2 = y * yst
			x2 = x * xst
			row[:] = [x2**2, y2**2, x2*y2, x2, y2, 1] # FIXME: change to (Y|X)SIZE
	else:
		for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
			row = rows.next()
			y2 = y * yst
			x2 = x * xst
			row[:] = [x2**2, y2**2, x2*y2, x2, y2] # FIXME: change to (Y|X)SIZE, needs proj4



class OrbitalCorrectionError(Exception):
	pass
