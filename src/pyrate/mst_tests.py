'''
TODO
Author: Ben Davies, ANUSF
'''

import unittest
from numpy import nan
from itertools import product

from mst import mst_matrix, default_mst
import algorithm
from algorithm_tests import sydney_test_setup, MockIfg
from orbital_tests import sydney5_mock_ifgs


class MSTTests(unittest.TestCase):
	'''Basic verification of minimum spanning tree (MST) functionality.'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()
		self.epochs = algorithm.get_epochs(self.ifgs)


	def test_mst_matrix(self):
		# Verifies mst matrix function returns an array with of dicts in each cell
		# assumes pygraph is correct from its unit tests
		res = mst_matrix(self.ifgs, self.epochs)
		ys, xs = res.shape
		for y, x in product(xrange(ys), xrange(xs)):
			r = res[y,x]
			self.assertTrue(hasattr(r, "keys"))
			self.assertTrue(len(r) <= len(self.epochs.dates))


	def test_default_mst(self):
		ifgs = sydney5_mock_ifgs() # no nodes should be lost
		dates = [i.DATE12 for i in ifgs]

		res = default_mst(ifgs)
		num_edges = len(res.keys())
		self.assertEqual(num_edges, len(ifgs))

		# test edges, note node order can be reversed
		for edge in res.iteritems():
			self.assertTrue(edge in dates or (edge[1],edge[0]) in dates)

		# check all nodes exist in this default tree
		mst_dates = set(res.keys() + res.values())
		for i in ifgs:
			for node in i.DATE12:
				self.assertTrue(node in mst_dates)


	def test_partial_nan_pixel_stack(self):
		# Ensure a limited # of cells gives in a smaller node tree

		num_coherent = 3

		def test_result():
			res = mst_matrix(self.mock_ifgs, self.epochs)
			self.assertEqual(len(res[0,0]), num_coherent)

		self.mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
		for m in self.mock_ifgs[num_coherent:]:
			m.phase_data[:] = nan
		test_result()

		# fill in more nans leaving only one ifg
		for m in self.mock_ifgs[1:num_coherent]:
			m.phase_data[:] = nan
		num_coherent = 1
		test_result()


	def test_all_nan_pixel_stack(self):
		self.mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
		for m in self.mock_ifgs:
			m.phase_data[:] = nan

		res = mst_matrix(self.mock_ifgs, self.epochs)
		exp = [nan]

		shape = (self.mock_ifgs[0].FILE_LENGTH, self.mock_ifgs[0].WIDTH)
		self.assertTrue(res.shape == shape)
		self.assertEqual(exp, res)
