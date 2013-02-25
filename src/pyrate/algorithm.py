'''
Collection of algorithms for PyRate.
Author: Ben Davies, ANUSF
'''

from math import pi
from numpy import sin, cos, unique, histogram

from shared import EpochList


# constants
MM_PER_METRE = 1000


def wavelength_to_mm(data, wavelength):
	"""Converts ROIPAC phase from metres to millimetres"""
	return data * MM_PER_METRE * (wavelength / (4 * pi))


def los_conversion(phase_data, unit_vec_component):
	'''Converts phase from LOS to horizontal/vertical components. Args are
	numpy arrays.'''

	# NB: currently not tested as implementation is too simple
	return phase_data * unit_vec_component


def unit_vector(incidence, azimuth):
	'''Returns unit vector tuple (east_west, north_south, vertical).'''
	vertical = cos(incidence)
	north_south = sin(incidence) * sin(azimuth)
	east_west = sin(incidence) * cos(azimuth)
	return east_west, north_south, vertical


def get_epochs(ifgs):
	'''Returns an EpochList derived from all given interferograms.'''
	masters = [i.MASTER for i in ifgs]
	slaves = [i.SLAVE for i in ifgs]

	combined = masters + slaves
	dates, n = unique(combined, False, True)
	repeat, _ = histogram(n, bins=len(set(n)))

	# absolute span for each date from the zero/start point
	span = [ (dates[i] - dates[0]).days / 365.25 for i in range(len(dates)) ]
	return EpochList(dates, repeat, span)


def master_slave_ids(dates):
	'''Returns dict of 'date:unique ID' for each date in dates. IDs are sorted in
	oldest to newest date order, starting at 0.'''
	dset = sorted(set(dates))
	return dict([(date_, i) for i, date_ in enumerate(dset)])
