'''
Collection of algorithms for PyRate.
Author: Ben Davies, ANUSF
'''

from math import pi
from numpy import sin, cos, unique, histogram

from shared import EpochList, IfgException


# constants
MM_PER_METRE = 1000



def least_squares_covariance(A,b,V):
	"""TODO
	known as lscov() in MATLAB)

	A: design matrix
	b: observations (vector of phase values)
	v: covariances (weights)
	"""


	#%LSCOV	Least squares solution in the presence of known covariance.
	#%	X = LSCOV(A,b,V) returns the vector X that minimizes
	#%	(A*X-b)'*inv(V)*(A*X-b) for the case in which length(b) > length(X).
	#%	This is the over-determined least squares problem with covariance V.
	#%	The solution is found without needing to invert V which is a square
	#%	symmetric matrix with dimensions equal to length(b).
	#%
	#%	The classical linear algebra solution to this problem is:
	#%
	#%	    x = inv(A'*inv(V)*A)*A'*inv(V)*b
	#%
	#%	See also SLASH, NNLS, QR.

	#%	Reference:
	#%	    G. Strang, "Introduction to Applied Mathematics",
	#%	    Wellesley-Cambridge, p. 398, 1986.
	#%	L. Shure 3-31-89
	#%	Copyright (c) 1984-94 by The MathWorks, Inc.

	m, n = A.shape
	if m <= n:
		raise ValueError('Problem must be over-determined')

	raise NotImplementedError

	#% Matlab's 'qr' function is an Orthogonal-triangular Decomposition. Scipy
	#% has a qr funciton in the linalg package

	#[q,r] = qr(A);
	#efg = q'*V*q;
	#e = efg(1:n,1:n);
	#g = efg(n+1:m,n+1:m);
	#cd = q'*b;
	#f = efg(1:n,n+1:m);
	#c = cd(1:n);
	#d = cd(n+1:m);
	#r = r(1:n,1:n);

	#% The 'slash' command in matlab is a matrix division. a\b is roughly equal
	#% to but more efficient than inv(a)*b

	#x = r\(c-f*(g\d));


def wavelength_to_mm(data, wavelength):
	"""Converts ROIPAC phase from radians to millimetres"""
	return data * MM_PER_METRE * (wavelength / (4 * pi))


def los_conversion(phase_data, unit_vec):
	'''
	Converts phase from line-of-sight (LOS) to horizontal/vertical components.
	phase_data - phase band data array (eg. ifg.phase_data)
	unit_vec - 3 component sequence, eg. [EW, NS, vertical]
	'''
	# NB: currently not tested as implementation is too simple
	return phase_data * unit_vec


def unit_vector(incidence, azimuth):
	'''Returns unit vector tuple (east_west, north_south, vertical).'''
	vertical = cos(incidence)
	north_south = sin(incidence) * cos(azimuth)
	east_west = sin(incidence) * sin(azimuth)
	return east_west, north_south, vertical


def ifg_date_lookup(ifgs, date_pair):
	'''
	Returns an Ifg which has a master/slave dates given in 'date_pair'.
	ifgs - list of Ifg objects to search in
	date_pair - a (datetime.date, datetime.date) tuple
	'''
	if len(date_pair) != 2:
		raise IfgException("Need (datetime.date, datetime.date) master/slave pair")

	# check master/slave dates are in order
	try:
		if date_pair[0] > date_pair[1]:
			date_pair = date_pair[1], date_pair[0]
	except:
		raise ValueError("Bad date_pair arg to ifg_date_lookup()")

	for i in ifgs:
		if date_pair == i.DATE12:
			return i

	raise ValueError("Cannot find Ifg with master/slave of %s" % str(date_pair))


def get_epochs(ifgs):
	'''Returns an EpochList derived from all given interferograms.'''
	combined = [i.MASTER for i in ifgs] + [i.SLAVE for i in ifgs]
	dates, n = unique(combined, False, True)
	repeat, _ = histogram(n, bins=len(set(n)))

	# absolute span for each date from the zero/start point
	span = [ (dates[i] - dates[0]).days / 365.25 for i in range(len(dates)) ]
	return EpochList(dates, repeat, span)


def get_all_epochs(ifgs):
	'''Returns sequence of all master and slave dates in given ifgs.'''
	return [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]


def get_epoch_count(ifgs):
	'''Returns count of unique epochs from sequence of ifgs'''
	return len(set(get_all_epochs(ifgs)))


def master_slave_ids(dates):
	'''
	Returns dict of 'date:unique ID' for each date in 'dates'. IDs are ordered
	from oldest to newest, starting at 0. Replaces ifglist.mas|slvnum in Pirate.
	'''
	dset = sorted(set(dates))
	return dict([(date_, i) for i, date_ in enumerate(dset)])
