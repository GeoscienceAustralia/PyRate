'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from os.path import join

from numpy.linalg import pinv
from numpy import nan, isnan, array, reshape, float32
from numpy import empty, ones, zeros, meshgrid
from numpy.testing import assert_array_equal, assert_array_almost_equal

import algorithm
from shared import Ifg
from orbital import OrbitalCorrectionError, orbital_correction
from orbital import get_design_matrix, get_network_design_matrix
from orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, QUADRATIC
from tests_common import sydney5_ifgs, MockIfg, IFMS5



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
			'''Helper method for result verification'''
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



class NetworkDesignMatrixErrorTests(unittest.TestCase):
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


class NetworkDesignMatrixTests(unittest.TestCase):

	def setUp(self):
		self.ifgs = sydney5_ifgs()


	def test_design_matrix_shape(self):
		# verify shape of design matrix is correct
		num_ifgs = len(self.ifgs)
		num_epochs = 6
		exp_num_rows = self.ifgs[0].num_cells * num_ifgs

		# with offsets
		for num_params, offset in zip((2, 3), (False, True)):
			exp_num_cols = num_epochs * num_params
			act_dm = get_network_design_matrix(self.ifgs, PLANAR, offset)
			self.assertEqual((exp_num_rows, exp_num_cols), act_dm.shape)

		# quadratic method
		for num_params, offset in zip((5, 6), (False, True)):
			exp_num_cols = num_epochs * num_params
			act_dm = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
			self.assertEqual((exp_num_rows, exp_num_cols), act_dm.shape)


	def test_network_design_matrix(self):
		# verify creation of sparse matrix comprised of smaller design matricies
		# TODO: do master/slave indices need to be sorted? Sorted = less ambiguity
		date_ids = get_date_ids(self.ifgs)

		# try all combinations of methods with and without the offsets
		for ncoef in [2, 3, 5, 6]:
			offset = ncoef in (3, 6)
			mth = PLANAR if ncoef < 5 else QUADRATIC
			act_dm = get_network_design_matrix(self.ifgs, mth, offset)
			self.assertNotEqual(act_dm.ptp(), 0)

			for i, ifg in enumerate(self.ifgs):
				# FIXME: sizes instead of steps
				exp_dm = unittest_dm(ifg, ifg.X_STEP, ifg.Y_STEP, mth, offset)

				# use slightly refactored version of Hua's code to test
				ib1 = i * ifg.num_cells # start row for subsetting the sparse matrix
				ib2 = (i+1) * ifg.num_cells # last row of subset of sparse matrix
				jbm = date_ids[ifg.MASTER] * ncoef # starting row index for master
				jbs = date_ids[ifg.SLAVE] * ncoef # row start for slave
				assert_array_almost_equal(-exp_dm, act_dm[ib1:ib2, jbm:jbm+ncoef])
				assert_array_almost_equal(exp_dm, act_dm[ib1:ib2, jbs:jbs+ncoef])


	def test_network_correct_planar(self):
		'''Verifies planar form of network method of correction'''
		for i in self.ifgs:
			i.open()

		err = 3
		self.ifgs[0].phase_data[0, 0:err] = nan # add NODATA as rasters are complete

		# reshape phase data to vectors
		num_ifgs = len(self.ifgs)
		data = empty(num_ifgs * self.ifgs[0].num_cells, dtype=float32)
		for i, ifg in enumerate(self.ifgs):
			st = i * ifg.num_cells
			data[st:st + ifg.num_cells] = ifg.phase_data.reshape(ifg.num_cells)

		dm = get_network_design_matrix(self.ifgs, PLANAR, False)[~isnan(data)]
		fd = data[~isnan(data)]
		self.assertEqual(len(dm), len(fd))
		self.assertEqual(len(dm), (num_ifgs * self.ifgs[0].num_cells) - err)

		params = pinv(dm, 1e-6) * fd
		act = orbital_correction(self.ifgs, PLANAR, NETWORK_METHOD, False)  # TODO: replace with a more internal function call?
		assert_array_almost_equal(act, params)
		# TODO: fwd correction
		# FIXME: with offsets


	def test_network_correct_quadratic(self):
		'''Verifies quadratic form of network method of correction'''
		for i in self.ifgs:
			i.open()

		self.ifgs[0].phase_data[3, 0:7] = nan # add NODATA as rasters are complete

		# reshape phase data to vectors
		data = empty(len(self.ifgs) * self.ifgs[0].num_cells, dtype=float32)
		for i, ifg in enumerate(self.ifgs):
			st = i * ifg.num_cells
			data[st:st + ifg.num_cells] = ifg.phase_data.reshape(ifg.num_cells)

		dm = get_network_design_matrix(self.ifgs, QUADRATIC, False)[~isnan(data)]
		params = pinv(dm, 1e-6) * data[~isnan(data)]

		act = orbital_correction(self.ifgs, QUADRATIC, NETWORK_METHOD, False)  # TODO: replace with a more internal function call?
		assert_array_almost_equal(act, params)
		# TODO: fwd correction
		# FIXME: with offsets


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


def get_date_ids(ifgs):
	'''Returns master/slave date IDs from the given Ifgs'''
	dates = []
	for ifg in ifgs:
		dates += list(ifg.DATE12)
	return algorithm.master_slave_ids(dates)


if __name__ == "__main__":
	unittest.main()
