'''
TODO

Created on 31/3/13
@author: Ben Davies
'''

from numpy import zeros, float32
from itertools import product


# Orbital correction
# 0) Config file stuff:
#		Add the orbital params
#		orbfit = 0: off, 1=do it
#		orbfitmethod = 1 or 2 (see conf file)  (AKA independent and networked methods)
#		orbrefest = orb_ref_est, 0=off, 1=on (IGNORE for the time being - something to do with removing the constant param 'c' from y = mx+c)
#		orbmaskflag = IGNORE for now. Used to mask out some patches, eg. if there is an earthquake signal somewhere.
#
# 1) resample/multilook the data to have less pixels (for comp efficiency)
# optional - vu can handle the larger arrays (new sys=?) so not mandatory
#
# 2) design matrix
#
# 3) linear regression - get params
#
# 4) forward calculation
#		create 2D orbital correctionl layer from the model params


def independent_design_matrix(ifg):
	'''Returns design matrix with 2 columns for linear model parameters'''

	# init design matrix
	shape = ((ifg.WIDTH * ifg.FILE_LENGTH), 2)
	data = zeros(shape, dtype=float32)
	idata = iter(data)

	# apply positional parameter values
	for y,x in product(xrange(ifg.FILE_LENGTH), xrange(ifg.WIDTH)):
		row = idata.next()
		row[:] = [y * ifg.Y_STEP, x * ifg.X_STEP]

	return data


def forward_calculation(design_mat, model):
	raise NotImplementedError
