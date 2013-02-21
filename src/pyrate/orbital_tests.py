'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from glob import glob
from numpy import nan, isnan, array, reshape, ones, zeros, float32
from numpy.testing import assert_array_equal, assert_array_almost_equal

from shared import Ifg
from orbital import orbital_correction, get_design_matrix, get_network_design_matrix
from orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, QUADRATIC
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

		# no offsets
		design_mat = get_design_matrix(m, 1, offset=False)
		assert_array_almost_equal(design_mat, self.designm)

		# with offset
		design_mat = get_design_matrix(m, 1, True)
		nys, nxs = self.designm.shape
		nxs += 1 # for offset col
		exp = ones((nys,nxs), dtype=float32)
		exp[:,:2] = self.designm
		assert_array_almost_equal(design_mat, exp)


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

		# no offset
		design_mat = get_design_matrix(m, degree=2, offset=False)
		assert_array_almost_equal(design_mat, exp_dm, decimal=3)

		# with offset
		design_mat = get_design_matrix(m, degree=2, offset=True)
		nys, nxs = exp_dm.shape
		nxs += 1 # for offset col
		exp2 = ones((nys,nxs), dtype=float32)
		exp2[:,:5] = exp_dm
		assert_array_almost_equal(design_mat, exp2, decimal=3)


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

		# test both models with no offsets TODO: magic numbers
		corrections = orbital_correction(ifgs, degree=1, method=1, offset=False)
		test_results()
		corrections = orbital_correction(ifgs, degree=2, method=1, offset=False)
		test_results()

		# test both with offsets
		corrections = orbital_correction(ifgs, degree=1, method=1)
		test_results()
		corrections = orbital_correction(ifgs, degree=2, method=1)
		test_results()


class OrbitalCorrectionNetwork(unittest.TestCase):
	'''TODO'''

	# TODO: check DM/correction with NaN rows

	# TODO
	#def test_invalid_ifgs_arg(self):
		# pass in list of ifgs of len=1 (what is the minimum required? 2 ifgs/3 epochs?)

	#def test_invalid_degree_arg(self):
		#pass

	def test_network_design_matrix_planar(self):
		# verify creation of sparse matrix comprised of smaller design matricies

		NUM_IFGS = 6
		_, ifgs = sydney_test_setup()
		ifgs = ifgs[:NUM_IFGS]

		# test planar option
		NUM_PARAMS = 3
		NUM_EPOCHS = 7 # number of ifgs + 1
		head = ifgs[0]
		exp_num_rows = head.FILE_LENGTH * head.WIDTH * NUM_IFGS
		exp_num_cols = NUM_EPOCHS * NUM_PARAMS
		exp_shape = (exp_num_rows, exp_num_cols)

		exp_dm = zeros(exp_shape, dtype=float32) # FIXME
		act_dm = get_network_design_matrix(ifgs, PLANAR, True)
		assert_array_equal(exp_dm, act_dm)


	# TODO:
	#def test_network_design_matrix_quadratic(self):
