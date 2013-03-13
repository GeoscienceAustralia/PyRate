'''
Collection of algorithms for PyRate.
Author: Ben Davies, ANUSF
'''

from math import pi, floor
from numpy import sin, cos, unique, histogram
import pyproj

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


def ifg_date_lookup(ifgs, date_pair):
	'''Returns an Ifg which has a master/slave dates given in 'date_pair'.
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


def master_slave_ids(dates):
	'''Returns dict of 'date:unique ID' for each date in dates.
	IDs are sorted in	oldest to newest date order, starting at 0.
	'''
	dset = sorted(set(dates))
	return dict([(date_, i) for i, date_ in enumerate(dset)])


def utm_zone(longitude):
	'''Returns basic UTM zone for given longitude in degrees. Currently does NOT
	handle the sub-zoning around Scandanavian countries.
	See http://www.dmap.co.uk/utmworld.htm
	'''
	if longitude == 180:
		return 60.0
	return floor((longitude + 180)/6.0) + 1


def cell_size(lat, lon, x_step, y_step):
	'''Converts X|Y_STEP in degrees to X & Y cell length/width in metres.
	lat - latitude in degrees
	lon - longitude in degrees
	x_step - horizontal step size in degrees
	y_step - vertical step size in degrees
	'''
	if lat > 84.0 or lat < -80:
		msg = "No UTM zone for polar region: > 84 degrees N or < 80 degrees S"
		raise ValueError(msg)

	zone = utm_zone(lon)
	p0 = pyproj.Proj(proj='latlong', ellps='WGS84')
	p1 = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')
	assert p0.is_latlong()
	assert not p1.is_latlong()

	x0, y0 = pyproj.transform(p0, p1, lon, lat)
	x1, y1 = pyproj.transform(p0, p1, lon + x_step, lat + y_step)
	return (abs(e) for e in (x1 - x0, y1 - y0))
