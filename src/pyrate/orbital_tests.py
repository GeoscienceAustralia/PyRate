'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from glob import glob
from numpy import nan, isnan, array, reshape, zeros, float32
from numpy.testing import assert_array_equal, assert_array_almost_equal

from shared import Ifg
from orbital import orbital_correction
from orbital import get_design_matrix
from algorithm_tests import MockIfg, sydney_test_setup



class OrbitalTests(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

	def setUp(self):
		xstep, ystep = 0.6, 0.7  # fake cell sizes

		shape = (6, 2)
		designm = zeros(shape) # planar design matrix for 2x3 orig array
		designm[0] = [0, 0]
		designm[1] = [0, xstep]
		designm[2] = [ystep, 0]
		designm[3] = [ystep, xstep]
		designm[4] = [2*ystep, 0]
		designm[5] = [2*ystep, xstep]

		self.designm = designm
		self.xstep = xstep
		self.ystep = ystep


	def test_design_matrix_planar(self):
		_, ifgs = sydney_test_setup()

		xs, ys = 2, 3
		m = MockIfg(ifgs[0], xs, ys)
		m.X_STEP = self.xstep
		m.Y_STEP = self.ystep
		design_mat = get_design_matrix(m, 1)
		assert_array_almost_equal(design_mat, self.designm)


	def test_design_matrix_quadratic(self):
		dm = []
		ys, xs = (3, 5)
		for y in xrange(ys):
			for x in xrange(xs):
				dm.append([(x * self.xstep)**2,
									(y * self.ystep)**2,
									(x * self.xstep) * (y * self.ystep),
									x * self.xstep,
									y * self.ystep])

		exp_dm = array(dm, dtype=float32) # planar design matrix

		# get design matrix from an ifg
		ifg = Ifg("../../tests/sydney_test/obs/geo_060619-061002.unw")
		ifg.open()
		m = MockIfg(ifg, xs, ys)
		m.X_STEP = self.xstep
		m.Y_STEP = self.ystep

		design_mat = get_design_matrix(m, degree=2)
		assert_array_almost_equal(design_mat, exp_dm, decimal=3)


	def test_ifg_to_vector(self):
		# test numpy reshaping order
		ifg = zeros((2, 3), dtype=float32)
		a = [24, 48, 1000]
		b = [552, 42, 68]
		ifg[0] = a
		ifg[1] = b
		res = reshape(ifg, (len(a+b)))
		exp = array(a + b)
		assert_array_equal(exp, res)

		# ensure order is maintained when converting back
		tmp = reshape(res, ifg.shape)
		assert_array_equal(tmp, ifg)


	def test_orbital_correction(self):

		def test_results():
			for i, c in zip(ifgs, corrections):
				# are corrections same size as the original array?
				ys, xs = c.shape
				self.assertEqual(i.FILE_LENGTH, ys)
				self.assertEqual(i.WIDTH, xs)
				self.assertFalse(isnan(i.phase_data).all())
				self.assertFalse(isnan(c).all())
				self.assertTrue(c.ptp() != 0)
				# TODO: do the results need to be checked at all?

		_, ifgs = sydney_test_setup()[:5]
		ifgs[0].phase_data[1, 1:3] = nan # add some NODATA

		corrections = orbital_correction(ifgs, degree=1, method=1)
		test_results()

		# test quadratic model
		corrections = orbital_correction(ifgs, degree=2, method=1)
		test_results()


class OrbitalCorrectionNetwork(unittest.TestCase):

	def test_network_design_matrix(self):
		# verify creation of sparse matrix comprised of smaller design matricies
		raise NotImplementedError
		#exp = None
		#act = get_network_design_matrix(ifgs)
		#assert_array_equal(exp, act)
