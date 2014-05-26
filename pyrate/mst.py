'''
Minimum Spanning Tree functionality for PyRate.
Contains functions to calculate MST using interferograms.

Author: Ben Davies, ANUSF
'''

from itertools import product
from numpy import array, nan, isnan, float32, ndarray

from algorithm import get_all_epochs, master_slave_ids
from pygraph.classes.graph import graph
from pygraph.algorithms.minmax import minimal_spanning_tree

# TODO: may need to implement memory saving row-by-row access
# TODO: document weighting by either Nan fraction OR variance


def _remove_root_node(mst):
	"""Discard pygraph's root node from MST dict to conserve memory."""
	for k in mst.keys():
		if mst[k] is None:
			del mst[k]
	return mst

# TODO: this does not appear to be used anywhere
def default_mst(ifgs, noroot=True):
	'''
	Returns default MST dict for the given Ifgs. The MST is calculated using a
	weighting based on the number of incoherent cells in the phase band.
	noroot: True removes the PyGraph default root node from the result.
	'''
	edges = [i.DATE12 for i in ifgs]
	epochs = master_slave_ids(get_all_epochs(ifgs)).keys()
	weights = [i.nan_fraction for i in ifgs]  # NB: other attrs for weights?

	g = _build_graph(epochs, edges, weights)
	mst = minimal_spanning_tree(g)
	if noroot:
		_remove_root_node(mst)
	return mst


def _build_graph(nodes, edges, weights, noroot=True):
	'Convenience graph builder function: returns a new graph obj'
	g = graph()
	g.add_nodes(nodes)
	for edge, weight in zip(edges, weights):
		g.add_edge(edge, wt=weight)
	
	return g


def mst_matrix(ifgs, epochs):
	'''
	Returns array of MST trees from a pixel-by-pixel MST. A MST is calculated
	for each individual pixel, ignoring NODATA values.
	ifgs: sequence of Ifg objs
	epochs: an EpochList object derived from the ifgs
	'''

	# locally cache all edges/weights for on-the-fly graph modification
	# make default MST to optimise result when no Ifg cells in a stack are nans
	edges = [i.DATE12 for i in ifgs]
	weights = [i.nan_fraction for i in ifgs]
	g = _build_graph(epochs.dates, edges, weights)
	dflt_mst = _remove_root_node(minimal_spanning_tree(g)) # default base MST

	# prepare source and dest data arrays
	num_ifgs = len(ifgs)
	data_stack = array([i.phase_data for i in ifgs], dtype=float32)
	mst_result = ndarray(shape=i.phase_data.shape, dtype=object)

	# create MSTs for each pixel in the ifg data stack
	for y, x in product(xrange(i.FILE_LENGTH), xrange(i.WIDTH)):
		values = data_stack[:, y, x] # select stack of all ifg values for a pixel
		nc = sum(isnan(values))

		# optimisations: use precreated results for all nans/no nans
		if nc == 0:
			mst_result[y, x] = dflt_mst
			continue
		elif nc == num_ifgs:
			mst_result[y, x] = nan
			continue

		# otherwise dynamically adjust graph, skipping edges where pixels are NaN
		for value, edge, weight in zip(values, edges, weights):
			if not isnan(value):
				if not g.has_edge(edge):
					g.add_edge(edge, wt=weight)
			else:
				if g.has_edge(edge):
					g.del_edge(edge)

		mst = _remove_root_node(minimal_spanning_tree(g))
		mst_result[y, x] = mst

	return mst_result
