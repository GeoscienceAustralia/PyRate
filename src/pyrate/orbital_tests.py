'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from numpy import zeros
from numpy.testing import assert_array_almost_equal

from orbital import independent_design_matrix, forward_calculation
from algorithm_tests import MockIfg, sydney_test_setup



class OrbitalTests(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

	def setUp(self):
		xstep, ystep = 0.6, 0.7  # fake cell sizes

		shape = (6,2)
		designm = zeros(shape)
		designm[0] = [0,0]
		designm[1] = [0,xstep]
		designm[2] = [ystep,0]
		designm[3] = [ystep,xstep]
		designm[4] = [2*ystep,0]
		designm[5] = [2*ystep,xstep]

		self.designm = designm
		self.xstep = xstep
		self.ystep = ystep


	def test_design_matrix(self):
		testdir, ifgs = sydney_test_setup()

		xs, ys = 2, 3
		m = MockIfg(ifgs[0], xs, ys)
		m.X_STEP = self.xstep
		m.Y_STEP = self.ystep

		design_mat = independent_design_matrix(m)
		assert_array_almost_equal(design_mat, self.designm)


	def test_linear_model(self):
		# TODO: this is testing orb forward
		model = [2, 4]  # a & b
		res = forward_calculation(self.designm, model)

		# TODO: this should be a model the same shape as the Ifg
		# TODO: should this be a ndarray or some other obj?
		expected = None
		self.assertEqual(expected, res)


	def test_complete(self):
		# create a data matrix D
		# pass D & design mat to some other func, which returns a model
		# do fwd modelling, subtract FM from D, test that result (this will need
		# the linalg/least sqaures in scipy...)
		raise NotImplementedError
