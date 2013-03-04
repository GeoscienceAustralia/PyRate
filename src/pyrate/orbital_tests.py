'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from os.path import join
from numpy import nan, isnan, array, reshape, where, float32
from numpy import empty, ones, zeros, meshgrid
from numpy.testing import assert_array_equal, assert_array_almost_equal
from scipy.linalg import pinv

import algorithm
from shared import Ifg
from orbital import OrbitalCorrectionError, orbital_correction
from orbital import get_design_matrix, get_network_design_matrix
from orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, QUADRATIC
from algorithm_tests import MockIfg



class IndependentDesignMatrixTests(unittest.TestCase):
	'''Tests for verifying various forms of design matricies'''

	def setUp(self):
		# fake cell sizes
		self.xsize = 0.6
		self.ysize = 0.7
		self.ifg = Ifg("../../tests/sydney_test/obs/geo_060619-061002.unw")
		self.ifg.open()


	def test_design_matrix_planar(self):
		m = MockIfg(self.ifg, 2, 3)
		m.X_STEP = self.xsize # FIXME: Ifg class needs to use custom X_SIZE
		m.Y_STEP = self.ysize # FIXME: Y_SIZE

		# test with and without offsets option
		exp = unittest_dm(m, self.xsize, self.ysize, PLANAR, True)
		assert_array_almost_equal(exp, get_design_matrix(m, PLANAR, True))
		assert_array_almost_equal(exp[:, :-1], get_design_matrix(m, PLANAR, False))


	def test_design_matrix_quadratic(self):
		m = MockIfg(self.ifg, 3, 5)
		m.X_STEP = self.xsize
		m.Y_STEP = self.ysize

		# test with and without offsets option
		exp_dm = unittest_dm(m, self.xsize, self.ysize, QUADRATIC, True)
		design_mat = get_design_matrix(m, QUADRATIC, False) # no offset
		assert_array_almost_equal(exp_dm[:, :-1], design_mat)
		assert_array_almost_equal(exp_dm, get_design_matrix(m, QUADRATIC, True))


class OrbitalCorrection(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

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


	def test_independent_correction(self):
		# verify top level orbital correction function

		def test_results():
			for i, c in zip(ifgs, corrections):
				# are corrections same size as the original array?
				ys, xs = c.shape
				self.assertEqual(i.FILE_LENGTH, ys)
				self.assertEqual(i.WIDTH, xs)
				self.assertFalse(isnan(i.phase_data).all())
				self.assertFalse(isnan(c).all())
				self.assertTrue(c.ptp() != 0) # ensure range of values in grid
				# TODO: do the results need to be checked at all?

		base = "../../tests/sydney_test/obs"
		ifgs = [Ifg(join(base, p)) for p in IFMS5.split()]
		for ifg in ifgs:
			ifg.open()

		ifgs[0].phase_data[1, 1:3] = nan # add some NODATA

		# test both models with no offsets
		for m in [PLANAR, QUADRATIC]:
			corrections = orbital_correction(ifgs, m, INDEPENDENT_METHOD, False)
			test_results()

		# test both with offsets
		for m in [PLANAR, QUADRATIC]:
			corrections = orbital_correction(ifgs, m, INDEPENDENT_METHOD)
			test_results()



class NetworkDesignMatrixTests(unittest.TestCase):
	'''Tests for the networked correction method'''
	# TODO: check correction with NaN rows
	# FIXME: check location of constant column - its at the end in HW's version

	def test_invalid_ifgs_arg(self):
		oex = OrbitalCorrectionError
		args = (PLANAR, True) # some default args
		for invalid in ([], [None]):
			self.assertRaises(oex, get_network_design_matrix, invalid, *args)
		# TODO: what is the minimum required? 2 ifgs/3 epochs?)


	def test_invalid_degree_arg(self):
		oex = OrbitalCorrectionError
		ifgs = sydney5_ifgs()

		for deg in range(-5, 1):
			self.assertRaises(oex, get_network_design_matrix, ifgs, deg, True)
		for deg in range(3, 6):
			self.assertRaises(oex, get_network_design_matrix, ifgs, deg, True)


	def test_design_matrix_shape(self):
		# verify shape of design matrix is correct
		ifgs = sydney5_ifgs()
		num_ifgs = len(ifgs)
		num_epochs = 6
		exp_num_rows = ifgs[0].num_cells * num_ifgs

		# with offsets
		for num_params, offset in zip((2,3), (False, True)):
			exp_num_cols = num_epochs * num_params
			act_dm = get_network_design_matrix(ifgs, PLANAR, offset)
			self.assertEqual((exp_num_rows, exp_num_cols), act_dm.shape)

		# quadratic method
		for num_params, offset in zip((5,6), (False, True)):
			exp_num_cols = num_epochs * num_params
			act_dm = get_network_design_matrix(ifgs, QUADRATIC, offset)
			self.assertEqual((exp_num_rows, exp_num_cols), act_dm.shape)


	def test_network_design_matrix(self):
		# verify creation of sparse matrix comprised of smaller design matricies
		# TODO: do master/slave indices need to be sorted? Sorted = less ambiguity
		ifgs = sydney5_ifgs()
		date_ids = get_date_ids(ifgs)

		# try all combinations of methods with and without the offsets
		for ncoef in [2, 3, 5, 6]:
			offset = ncoef in (3, 6)
			mth = PLANAR if ncoef < 5 else QUADRATIC
			act_dm = get_network_design_matrix(ifgs, mth, offset)
			self.assertNotEqual(act_dm.ptp(), 0)

			for i, ifg in enumerate(ifgs):
				# FIXME: sizes instead of steps
				exp_dm = unittest_dm(ifg, ifg.X_STEP, ifg.Y_STEP, mth, offset)

				# use slightly refactored version of Hua's code to test
				ib1 = i * ifg.num_cells # start row for subsetting the sparse matrix
				ib2 = (i+1) * ifg.num_cells # last row of subset of sparse matrix
				jbm = date_ids[ifg.MASTER] * ncoef # starting row index for master
				jbs = date_ids[ifg.SLAVE] * ncoef # row start for slave
				assert_array_almost_equal(-exp_dm, act_dm[ib1:ib2, jbm:jbm+ncoef])
				assert_array_almost_equal(exp_dm, act_dm[ib1:ib2, jbs:jbs+ncoef])


	def test_network_correction_planar(self):
		ifgs = sydney5_mock_ifgs()
		for i in ifgs: i.open()

		ERR = 3
		ifgs[0].phase_data[0,0:ERR] = nan # add NODATA as rasters are complete

		# reshape phase data to vectors
		num_ifgs = len(ifgs)
		data = empty(num_ifgs * ifgs[0].num_cells, dtype=float32)
		for i, ifg in enumerate(ifgs):
			st = i * ifg.num_cells
			data[st:st + ifg.num_cells] = ifg.phase_data.reshape(ifg.num_cells)

		dm = get_network_design_matrix(ifgs, PLANAR, offset=False)
		dm = dm[~isnan(data)]
		fd = data[~isnan(data)]
		self.assertTrue(len(dm) == len(fd) == (num_ifgs * ifg.num_cells) - ERR)

		# inversion
		params = pinv(dm, 1e-6) * fd

		# TODO: replace with a more general function call
		from orbital import _get_net_correction
		act = _get_net_correction(ifgs, PLANAR, False)
		assert_array_almost_equal(act, params)

		# TODO: fwd correction


# FIXME: add derived field to ifgs to convert X|Y_STEP degrees to metres
def unittest_dm(ifg, xs, ys, degree, offset=False):
	'''Convenience function to create design matrices'''
	ncoef = 2 if degree == PLANAR else 5
	if offset is True:
		ncoef += 1

	out = ones((ifg.num_cells, ncoef), dtype=float32)
	x, y = meshgrid(range(ifg.WIDTH), range(ifg.FILE_LENGTH))
	x = x.reshape(ifg.num_cells) * xs
	y = y.reshape(ifg.num_cells) * ys

	if degree == PLANAR:
		out[:, 0] = x
		out[:, 1] = y
	elif degree == QUADRATIC:
		out[:, 0] = x**2
		out[:, 1] = y**2
		out[:, 2] = x * y
		out[:, 3] = x
		out[:, 4] = y
	else:
		raise Exception("Degree is invalid")

	return out


# small dummy ifg list to limit overall # of ifgs
IFMS5 = """geo_060828-061211.unw
geo_061106-061211.unw
geo_061106-070115.unw
geo_061106-070326.unw
geo_070326-070917.unw
"""

def sydney5_ifgs():
	'''Convenience func to return a subset of 5 linked Ifgs from the testdata'''
	base = "../../tests/sydney_test/obs"
	return [Ifg(join(base, p)) for p in IFMS5.split()]


def sydney5_mock_ifgs(xs=3, ys=4):
	'''Returns smaller mocked version of sydney Ifgs for testing'''
	ifgs = sydney5_ifgs()
	for i in ifgs: i.open()
	mocks = [MockIfg(i, xs, ys) for i in ifgs]
	for i,m in zip(ifgs, mocks):
		m.phase_data = i.phase_data[:ys,:xs]
		del m.nan_fraction

	return mocks


def get_date_ids(ifgs):
	'''Returns master/slave date IDs from the given Ifgs'''
	dates = []
	for ifg in ifgs:
		dates += list(ifg.DATE12)
	return algorithm.master_slave_ids(dates)
