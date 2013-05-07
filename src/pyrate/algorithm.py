'''
Collection of algorithms for PyRate.
Author: Ben Davies, ANUSF
'''

from math import pi
from numpy import sin, cos, unique, histogram

from shared import EpochList, IfgException


# constants
MM_PER_METRE = 1000


def wavelength_to_mm(data, wavelength):
	"""Converts ROIPAC phase from metres to millimetres"""
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
	north_south = sin(incidence) * sin(azimuth)
	east_west = sin(incidence) * cos(azimuth)
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
	if date_pair[0] > date_pair[1]:
		date_pair = date_pair[1], date_pair[0]

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
