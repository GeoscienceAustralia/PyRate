'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from glob import glob
from numpy import array, reshape, zeros, float32
from numpy.testing import assert_array_equal, assert_array_almost_equal

from shared import Ifg
from orbital import orbital_correction_np
from orbital import independent_design_matrix, forward_calculation
from algorithm_tests import MockIfg, sydney_test_setup



class OrbitalTests(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

	def setUp(self):
		xstep, ystep = 0.6, 0.7  # fake cell sizes

		shape = (6,2)
		designm = zeros(shape) # design matrix
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


	def test_ifg_to_vector(self):
		# test numpy reshaping order
		ifg = zeros((2,3), dtype=float32)
		a = [24, 48, 1000]
		b = [552, 42, 68]
		ifg[0] = a
		ifg[1] = b
		res = reshape(ifg, (len(a+b)))
		exp = array(a + b)
		assert_array_equal(exp, res)


	def test_orbital_correction_np(self):
		paths = sorted(glob("../../tests/sydney_test/obs/geo*.unw"))[:5]
		ifgs = [Ifg(p) for p in paths]
		[i.open() for i in ifgs]
		corrections = orbital_correction_np(ifgs, degree=1, method=1)
		for i,c in zip(ifgs, corrections):
			# are corrections same size as the original array?
			ys, xs = c.shape
			self.assertEqual(i.FILE_LENGTH, ys)
			self.assertEqual(i.WIDTH, xs)
			self.assertTrue(c.ptp() != 0)


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
