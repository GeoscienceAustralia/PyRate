'''
TODO
'''

from itertools import product
from numpy import array, nan, isnan, float32, ndarray

import algorithm
from pygraph.classes.graph import graph
from pygraph.algorithms.minmax import minimal_spanning_tree



def _remove_root_node(mst):
	"""Discard pygraph's root node from MST dict to conserve memory."""
	for k in mst.keys():
		if mst[k] is None:
			del mst[k]
	return mst


def default_mst(ifgs, noroot=True):
	'''Returns the default MST dict for the given Ifgs. False for noroot prevents
	the root node from being removed from the result.
	'''
	edges = [i.DATE12 for i in ifgs]
	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	epochs = algorithm.master_slave_ids(dates).keys()
	weights = [i.nan_fraction for i in ifgs]  # TODO: user specified attrs for weights?

	g = graph()
	g.add_nodes(epochs) # each acquisition is a node
	for edge, weight in zip(edges, weights):
		g.add_edge(edge, wt=weight)

	mst = minimal_spanning_tree(g)
	if noroot:
		_remove_root_node(mst)
	return mst


def mst_matrix(ifgs, epochs):
	'''Returns array of minimum spanning trees for the Ifgs.'''

	# TODO: implement rows memory saving option/ row by row access?

	# locally cache all edges/weights for on-the-fly graph modification
	edges = [i.DATE12 for i in ifgs]
	weights = [i.nan_fraction for i in ifgs]

	# make default MST to optimise result when no Ifg cells in a stack are nans
	g = graph()
	g.add_nodes(epochs.dates) # each acquisition is a node
	for edge, weight in zip(edges, weights):
		g.add_edge(edge, wt=weight)

	dmst = _remove_root_node(minimal_spanning_tree(g)) # default base MST

	# prepare source and dest data arrays
	# [i.phase_data for i in ifgs]
	num_ifgs = len(ifgs)
	data_stack = array([i.phase_data for i in ifgs], dtype=float32)
	mst_result = ndarray(shape=(i.FILE_LENGTH, i.WIDTH), dtype=object)

	# create MSTs for each pixel in the ifg data stack
	for y, x in product(xrange(i.FILE_LENGTH), xrange(i.WIDTH)):
		values = data_stack[:, y, x] # select stack of all ifg values for a pixel
		nc = sum(isnan(values))

		# optimisations: use precreated results for all nans/no nans
		if nc == 0:
			mst_result[y, x] = dmst
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
