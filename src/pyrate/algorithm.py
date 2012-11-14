'''
Created
Author: Ben Davies
'''

from math import pi

from config import EpochList, get_epochs
from numpy import nan, isnan, sum, sin, cos, radians, unique
from pygraph.classes.graph import graph
from pygraph.algorithms.minmax import minimal_spanning_tree

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


def temp_mst(ifgs):
	'''TODO'''
	
	# TODO: top level loop to handle cell by cell MST
	# TODO: run MST on each pixel
	
	# TODO: cell by cell accesses? of do by rows/numpy arrays? 
	# TODO: class, tuple (or named tuple?) for each epoch with links to the different ifgs it can relate to
	# TODO: what to do in the event of there being not enough edges to connect the epochs? (ie. there are too many NaN cells?) 
	
	epochs = get_epochs(ifgs)
		
	g = graph()
	g.add_nodes(epochs.date)
	for i in ifgs:
		i.open()
		g.add_edge((i.MASTER, i.SLAVE), i.nan_fraction)
		
	# TODO: what happens with py-graph when a node isn't connected
	
	return minimal_spanning_tree(g)
