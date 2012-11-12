'''
Created
Author: Ben Davies
'''

from math import pi

from numpy import nan, isnan, sum, sin, cos, radians


# constants
MM_PER_METRE = 1000


def wavelength_to_mm(data, wavelength):
	"""Converts ROIPAC phase from metres to millimetres"""
	return data * MM_PER_METRE * (wavelength / (4 * pi))


# TODO: remove this as there is an equivalent func in Ifg class
def nan_fraction(data):
	"""Calculates fraction of cells in data that are NaN"""
	
	# calc number of cells by summing additional dimensions
	shape = data.shape
	ncells = float(shape[0])  # for float division
	if len(shape) > 1:
		for i in data.shape[1:]:
			ncells *= i
	
	nan_mask = isnan(data)
	nan_sum = sum(nan_mask)
	return nan_sum / ncells


def los_conversion(phase_data, unit_vec_component):
	'''Converts phase from LOS to horizontal/vertical components. Args are 
	numpy arrays.'''
	
	# NB: currently not tested as implementation is too simple 
	return phase_data * unit_vec_component 
	

def unit_vector(incidence, azimuth):
	vertical = cos(incidence)
	north_south = sin(incidence) * sin(azimuth)
	east_west = sin(incidence) * cos(azimuth)
	return east_west, north_south, vertical	
