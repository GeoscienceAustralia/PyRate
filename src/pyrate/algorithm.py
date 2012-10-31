'''
Created
Author: Ben Davies
'''

from math import pi

from numpy import nan, isnan, sum

# constants
MM_PER_METRE = 1000


def wavelength_to_mm(data, wavelength):
	"""Converts ROIPAC phase from metres to millimetres"""
	return data * MM_PER_METRE * (wavelength / (4 * pi))


def nan_fraction(data):
	"""Calculates fraction of cells in data that are NaN"""
	nan_mask = isnan(data)
	nan_sum = sum(nan_mask)

	# calc number of cells (sum additional dimensions as required)
	shape = data.shape
	ncells = float(shape[0])  # for float division
	if len(shape) > 1:
		for i in data.shape[1:]:
			ncells *= i

	return nan_sum / ncells
