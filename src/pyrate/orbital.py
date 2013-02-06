'''
TODO

Created on 31/3/13
@author: Ben Davies
'''

from itertools import product
from numpy import sum, reshape, zeros, float32
from numpy.linalg import lstsq


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


def orbital_correction(ifgs, degree, method):
	# TODO: get the orbital corrections
	# TODO: save corrected layers to new file or use intermediate arrays?
	raise NotImplementedError


# TODO: hide this function?
def orbital_correction_np(ifgs, degree, method):
	'''Calculates orbital error and returns layers of the forward correction.'''

	# TODO: offsets

	if method == NETWORK_METHOD:
		raise NotImplementedError

	# grab head ifg to prepare some data
	head = ifgs[0]
	orig = (head.FILE_LENGTH, head.WIDTH)

	sh = head.FILE_LENGTH * head.WIDTH # new shape for vectorised data
	dm = independent_design_matrix(head)

	# TODO: performance - add filtering for NODATA?
	models = [lstsq(dm, reshape(i.phase_data, sh))[0] for i in ifgs]
	corrections = [sum(dm * m, axis=1) for m in models] # d = ax + by
	return [reshape(c, orig) for c in corrections]


def independent_design_matrix(ifg):
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


def forward_calculation(design_mat, model):
	raise NotImplementedError
