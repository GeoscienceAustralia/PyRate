'''
TODO

Created on 31/3/13
@author: Ben Davies
'''

from itertools import product
from numpy import sum, where, nan, isnan, reshape, zeros, float32
from numpy.ma import masked_array

#from numpy.linalg import lstsq
from scipy.linalg import lstsq


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


def orbital_correction(ifgs, degree, method):

	# TODO: save corrected layers to new file or use intermediate arrays?
	# TODO: offsets

	if degree not in [PLANAR, QUADRATIC]:
		msg = "Invalid degree of %s for orbital correction" % degree
		raise OrbitalCorrectionError(msg)

	if method == NETWORK_METHOD:
		raise NotImplementedError

	return [_get_correction(i, degree) for i in ifgs]


def _get_correction(ifg, degree):
	'''Calculates and returns orbital correction array for an ifg'''

	# precalculate matrix shapes
	orig = ifg.phase_data.shape
	vsh = orig[0] * orig[1] # full shape of vectorised form of data

	# vectorise data
	vphase = reshape(ifg.phase_data, vsh) # contains NODATA

	# TODO: refactor DM into single func?
	if degree == PLANAR:
		dm = get_design_matrix(ifg)
	else:
		dm = get_design_matrix_quadratic(ifg)
	assert len(vphase) == len(dm)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model, _, rank, _ = lstsq(tmp, fd)

	lm = len(model)
	if degree == PLANAR:
		assert lm == 2
	else:
		assert lm == 5

	# TODO: assert rank == EXP, "Got rank of %s" % rank

	# calculate forward model & morph back to 2D
	tmp = sum(dm * model, axis=1) # d = ax + by
	correction = reshape(tmp, orig)
	return correction


def get_design_matrix(ifg):
	'''Returns design matrix with 2 columns for linear model parameters'''

	# init design matrix
	shape = ((ifg.WIDTH * ifg.FILE_LENGTH), 2)
	data = zeros(shape, dtype=float32)
	idata = iter(data)

	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)
	for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
		row = idata.next()
		row[:] = [y * ifg.Y_STEP, x * ifg.X_STEP] # FIXME: change to (Y|X)SIZE, needs proj4

	return data


def get_design_matrix_quadratic(ifg):
	# init design matrix
	shape = ((ifg.WIDTH * ifg.FILE_LENGTH), 5)
	data = zeros(shape, dtype=float32)
	idata = iter(data)

	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)
	yst, xst = ifg.Y_STEP, ifg.X_STEP
	for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
		row = idata.next()
		y2 = y * yst
		x2 = x * xst
		row[:] = [x2**2, y2**2, x2*y2, x2, y2] # FIXME: change to (Y|X)SIZE, needs proj4

	return data



class OrbitalCorrectionError(Exception):
	pass
