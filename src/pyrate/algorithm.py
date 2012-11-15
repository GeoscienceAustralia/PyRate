'''
Created
Author: Ben Davies
'''

from math import pi
from itertools import product

from numpy import nan, isnan, sum, sin, cos, radians, unique, histogram
from numpy import array, ndarray
from pygraph.classes.graph import graph
from pygraph.algorithms.minmax import minimal_spanning_tree

from shared import EpochList


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


def get_epochs(ifgs):
	masters = [i.MASTER for i in ifgs]
	slaves = [i.SLAVE for i in ifgs]

	combined = masters + slaves
	dates, n = unique(combined, False, True)
	repeat, _ = histogram(n, bins=len(set(n)))

	# absolute span for each date from the zero/start point
	span = [ (dates[i] - dates[0]).days / 365.25 for i in range(len(dates)) ]
	return EpochList(dates, repeat, span)


def temp_mst(ifgs):
	'''TODO'''

	# TODO: Test speed of cell by cell accesses?
	# TODO: implement rows memory saving option/ row by row access
	# TODO: implement minimum # edges/ifg threshold
	# TODO: does default overall MST need to be implemented for the default value?

	g = graph()
	epochs = get_epochs(ifgs)
	g.add_nodes(epochs.dates) # each acquisition is a node

	# cache all possible edges
	edges = [i.DATE12 for i in ifgs]
	weights = [i.nan_fraction for i in ifgs]
	data_stack = array([i.phase_band.ReadAsArray() for i in ifgs], dtype=object)
	mst_result = ndarray(shape=(i.FILE_LENGTH, i.WIDTH), dtype=object)

	# create MSTs for each pixel in the ifg data stack
	for y, x in product(xrange(i.FILE_LENGTH), xrange(i.WIDTH)):
		values = data_stack[:,y,x] # drill down through all ifgs for a pixel
		if (values == nan).all():
			raise NotImplementedError("All cells are NaN at (y=%s, x=%s)" % (y,x))

		# TODO: implement NaN threshold here
		# TODO: handle case of too few edges to connect the epochs? (ie. too many NaN cells)

		# dynamically adjust graph, removing edges where pixel is NaN
		for value, edge, weight in zip(values, edges, weights):
			if not isnan(value):
				if not g.has_edge(edge):
					g.add_edge(edge, wt=weight)
			else:
				if g.has_edge(edge):
					g.del_edge(edge)

		mst = minimal_spanning_tree(g)

		# discard root node (saves some memory)
		for k in mst.keys():
			if mst[k] is None:
				del mst[k]

		mst_result[y,x] = mst

	return mst_result



