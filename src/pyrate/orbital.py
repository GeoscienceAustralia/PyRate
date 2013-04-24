'''
Module for calculating orbital correction for interferograms.

Created on 31/3/13
@author: Ben Davies
'''

from numpy import empty, isnan, reshape, float32, squeeze
from numpy import dot, vstack, zeros, meshgrid
from scipy.linalg import lstsq
from numpy.linalg import pinv

from algorithm import master_slave_ids, get_all_epochs


# Orbital correction tasks
#
# 3) FIXME: REDO/REFACTOR TESTING:
#  Manually make small DMs with corect offsets cols etc and test against these!
#
# 4) forward calculation (orbfwd.m)
#		* create 2D orbital correction layer from the model params
#		* add handling for mlooked ifgs/correcting the full size ones
#		* correct ifgs in place + save


# TODO: options for multilooking
# 1) do the 2nd stage mlook at prepifg.py/generate up front, then delete in workflow after
# 2) refactor prep_ifgs() call to take input filenames and params & generate mlooked versions from that
#    this needs to be more generic to call at any point in the runtime.


# constants
INDEPENDENT_METHOD = 1
NETWORK_METHOD = 2

PLANAR = 1
QUADRATIC = 2


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
		# FIXME: handle mlooked and straight Ifgs

		if mlooked is not None:
			_validate_mlooked(mlooked, ifgs)
			return _get_net_correction(sub_mlooked, degree, offset) # TODO: fwd corr
		else:
			return _get_net_correction(ifgs, degree, offset) # TODO: fwd corr

	elif method == INDEPENDENT_METHOD:
		# FIXME: determine how to work this into the ifgs. Generate new Ifgs?
		# Update passed in ifgs and flag the corrections in metadata?
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


def get_num_params(degree, offset):
	'''Returns number of model parameters'''
	nparams = 2 if degree == PLANAR else 5
	if offset is True:
		nparams += 1  # eg. y = mx + offset
	return nparams


def _get_ind_correction(ifg, degree, offset):
	'''Calculates and returns orbital correction array for an ifg'''

	vphase = reshape(ifg.phase_data, ifg.num_cells) # vectorised, contains NODATA
	dm = get_design_matrix(ifg, degree, offset)
	assert len(vphase) == len(dm)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model, _, rank, _ = lstsq(tmp, fd)

	# calculate forward model & morph back to 2D
	correction = reshape(dot(dm, model), ifg.phase_data.shape)
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

	# get DM / clear out the NaNs based on obs
	tmp = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in ifgs])
	vphase = squeeze(tmp)
	dm = get_network_design_matrix(ifgs, degree, offset)
	assert len(vphase) == len(dm)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model = dot(pinv(tmp, 1e-6), fd)
	ncoef = get_num_params(degree, offset)
	coefs = [model[i:i+ncoef] for i in range(len(model))[::ncoef]]
	ids = master_slave_ids(get_all_epochs(ifgs))

	dm = get_design_matrix(ifgs[0], degree, offset)
	shp = ifgs[0].shape

	corrections = []
	for i in ifgs:
		par = coefs[ids[i.SLAVE]] - coefs[ids[i.MASTER]]
		corrections.append(dot(dm, par))  # TODO: plus offset?

	corrections = [e.reshape(shp) for e in corrections]

	# TODO apply corrections to Ifgs
	return corrections


def get_design_matrix(ifg, degree, offset):
	'''
	Returns design matrix with columns for model parameters.
	ifg - interferogram to base the DM on
	degree - PLANAR or QUADRATIC
	offset - True to include offset cols, otherwise False.
	'''
	# apply positional parameter values, multiply pixel coordinate by cell size to
	# get distance (a coord by itself doesn't tell us distance from origin)

	# init design matrix
	data = empty((ifg.num_cells, get_num_params(degree, offset)), dtype=float32)
	x, y = meshgrid(range(ifg.WIDTH), range(ifg.FILE_LENGTH))
	x = x.reshape(ifg.num_cells) * ifg.X_SIZE
	y = y.reshape(ifg.num_cells) * ifg.Y_SIZE

	if degree == PLANAR:
		data[:, 0] = x
		data[:, 1] = y
	elif degree == QUADRATIC:
		data[:, 0] = x**2
		data[:, 1] = y**2
		data[:, 2] = x * y
		data[:, 3] = x
		data[:, 4] = y

	if offset is True:
		data[:, -1] = 1

	return data


# TODO: should this be refactored under one get_design_matrix() func?
def get_network_design_matrix(ifgs, degree, offset):
	'''Returns a larger format design matrix for networked error correction.'''

	if degree not in [PLANAR, QUADRATIC]:
		raise OrbitalError("Invalid degree argument")

	num_ifgs = len(ifgs)
	if num_ifgs < 1:
		# can feasibly do correction on a single Ifg/2 epochs
		raise OrbitalError("Invalid number of Ifgs")

	# init design matrix
	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	num_epochs = len(set(dates)) # TODO: test this
	nparams = get_num_params(degree, offset)
	shape = [ifgs[0].num_cells * num_ifgs, nparams * num_epochs] # FIXME: has 1 too many cols for OFFSETS for PLANAR AND QUADRATIC
	data = zeros(shape, dtype=float32)

	# individual design matrices
	ids = master_slave_ids(dates)
	ncoef = get_num_params(degree, False) # only base level of coefficients
	offset_col = num_epochs * ncoef # base offset for the offset cols

	for i, ifg in enumerate(ifgs):
		tmp = get_design_matrix(ifg, degree, False) # DMs within full DM don't have extra col
		rs = i * ifg.num_cells # starting row
		m = ids[ifg.MASTER] * ncoef  # start col for master
		s = ids[ifg.SLAVE] * ncoef  # start col for slave
		data[rs:rs + ifg.num_cells, m:m + ncoef] = -tmp
		data[rs:rs + ifg.num_cells, s:s + ncoef] = tmp

		if offset:
			data[rs:rs + ifg.num_cells, offset_col + i] = 1  # init offset cols

	return data



class OrbitalError(Exception):
	'''Generic class for errors in orbital correction'''
	pass
