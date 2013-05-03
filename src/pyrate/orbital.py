'''
Module for calculating orbital correction for interferograms.

Created on 31/3/13
@author: Ben Davies
'''

from numpy import empty, isnan, reshape, float32, squeeze
from numpy import dot, vstack, zeros, median, meshgrid
from scipy.linalg import lstsq
from numpy.linalg import pinv

from algorithm import master_slave_ids, get_all_epochs, get_epoch_count


# Orbital correction tasks
#
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
	NB: this modifies the Ifgs in place.

	ifgs - list of Ifg objs to correct
	degree - PLANAR or QUADRATIC
	method - INDEPENDENT_METHOD or NETWORK_METHOD
	mlooked - sequence of multilooked Ifgs
	offset = True/False to include the constant/offset component
	'''

	if degree not in [PLANAR, QUADRATIC]:
		msg = "Invalid degree of %s for orbital correction" % degree
		raise OrbitalError(msg)

	if method not in [INDEPENDENT_METHOD, NETWORK_METHOD]:
		msg = "Unknown method '%s', need INDEPENDENT or NETWORK method" % method
		raise OrbitalError(msg)

	if method == NETWORK_METHOD:
		if mlooked is None:
			return _network_correction(ifgs, degree, offset)
		else:
			_validate_mlooked(mlooked, ifgs)
			return _network_correction(ifgs, degree, offset, mlooked)

	elif method == INDEPENDENT_METHOD:
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


def get_num_params(degree, offset=None):
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
	model = lstsq(tmp, fd)[0] # first arg is the model params

	# calculate forward model & morph back to 2D
	correction = reshape(dot(dm, model), ifg.phase_data.shape)
	return correction


def _network_correction(ifgs, degree, offset, m_ifgs=None):
	'''
	Returns the TODO
	ifgs - assumed to be Ifgs from a prior MST step
	degree - PLANAR or QUADRATIC
	offset - True/False for including TODO
	'''
	# get DM & filter out NaNs
	src_ifgs = ifgs if m_ifgs is None else m_ifgs
	vphase = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in src_ifgs])
	vphase = squeeze(vphase)
	dm = get_network_design_matrix(src_ifgs, degree, offset)

	# filter NaNs out before getting model
	tmp = dm[~isnan(vphase)]
	fd = vphase[~isnan(vphase)]
	model = dot(pinv(tmp, 1e-6), fd)

	ncoef = get_num_params(degree)
	ids = master_slave_ids(get_all_epochs(ifgs))
	coefs = [model[i:i+ncoef] for i in range(0, len(set(ids)) * ncoef, ncoef)]

	# create DM to expand into surface from params
	dm = get_design_matrix(ifgs[0], degree, offset=False)

	for i in ifgs:
		orb = dot(dm, coefs[ids[i.SLAVE]] - coefs[ids[i.MASTER]])
		orb = orb.reshape(ifgs[0].shape)

		# estimate offsets
		if offset:
			tmp = i.phase_data - orb
			i.phase_data -= (orb + median(tmp[~isnan(tmp)]))
		else:
			i.phase_data -= orb


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

	# TODO: performance test this vs np.concatenate (n by 1 cols)

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


def get_network_design_matrix(ifgs, degree, offset):
	'''
	Returns a larger format design matrix for networked error correction. This is
	a fullsize DM, including rows which relate to those of NaN cells.
	ifgs - sequence of interferograms
	degree - PLANAR or QUADRATIC mode
	offset - True to include offset cols, otherwise False.
	'''
	if degree not in [PLANAR, QUADRATIC]:
		raise OrbitalError("Invalid degree argument")

	nifgs = len(ifgs)
	if nifgs < 1:
		# can feasibly do correction on a single Ifg/2 epochs
		raise OrbitalError("Invalid number of Ifgs")

	# init design matrix
	nepochs = get_epoch_count(ifgs)
	ncoef = get_num_params(degree)
	shape = [ifgs[0].num_cells * nifgs, ncoef * nepochs]
	if offset:
		shape[1] += nifgs # add extra offset cols

	ndm = zeros(shape, dtype=float32)

	# individual design matrices
	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	ids = master_slave_ids(dates)
	offset_col = nepochs * ncoef # base offset for the offset cols

	tmp = get_design_matrix(ifgs[0], degree, False)
	for i, ifg in enumerate(ifgs):
		rs = i * ifg.num_cells # starting row
		m = ids[ifg.MASTER] * ncoef  # start col for master
		s = ids[ifg.SLAVE] * ncoef  # start col for slave
		ndm[rs:rs + ifg.num_cells, m:m + ncoef] = -tmp
		ndm[rs:rs + ifg.num_cells, s:s + ncoef] = tmp

		if offset:
			ndm[rs:rs + ifg.num_cells, offset_col + i] = 1  # init offset cols

	return ndm



class OrbitalError(Exception):
	'''Generic class for errors in orbital correction'''
	pass
