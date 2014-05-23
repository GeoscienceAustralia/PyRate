'''
Tests for Minimum Spanning Tree functionality in PyRate
Author: Ben Davies, ANUSF
'''

import unittest
from numpy import array, nan, isnan, where, sum as nsum
from itertools import product

from mst import mst_matrix, default_mst
import algorithm
from tests_common import MockIfg, sydney5_mock_ifgs, sydney_test_setup


class MSTTests(unittest.TestCase):
	'''Basic verification of minimum spanning tree (MST) functionality.'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()
		self.epochs = algorithm.get_epochs(self.ifgs)


	def test_mst_matrix(self):
		# Verifies mst matrix func returns array with dict/trees in each cell
		for i in self.ifgs[3:]:
			i.phase_data[0,1] = 0 # add a large stack of nans to one cell

		for i in self.ifgs:
			i.convert_to_nans()

		nepochs = len(self.epochs.dates)
		res = mst_matrix(self.ifgs, self.epochs)
		ys, xs = res.shape
		for y, x in product(xrange(ys), xrange(xs)):
			r = res[y,x]
			num_nodes = len(r)
			self.assertTrue(num_nodes < len(self.epochs.dates))

			stack = array([i.phase_data[y,x] for i in self.ifgs]) # 17 ifg stack
			self.assertTrue(0 == nsum(stack == 0)) # all 0s should be converted
			nc = nsum(isnan(stack))

			if nc == 0:
				self.assertTrue(num_nodes == (nepochs-1))
			elif nc > 5:
				# rough test: too many nans must reduce the total tree size
				self.assertTrue(num_nodes <= (17-nc) )


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
			res = mst_matrix(mock_ifgs, self.epochs)
			self.assertEqual(len(res[0,0]), num_coherent)

		mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
		for m in mock_ifgs[num_coherent:]:
			m.phase_data[:] = nan
		test_result()

		# fill in more nans leaving only one ifg
		for m in mock_ifgs[1:num_coherent]:
			m.phase_data[:] = nan
		num_coherent = 1
		test_result()


	def test_all_nan_pixel_stack(self):
		mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
		for m in mock_ifgs:
			m.phase_data[:] = nan

		res = mst_matrix(mock_ifgs, self.epochs)
		exp = [nan]

		shape = (mock_ifgs[0].FILE_LENGTH, mock_ifgs[0].WIDTH)
		self.assertTrue(res.shape == shape)
		self.assertEqual(exp, res)
