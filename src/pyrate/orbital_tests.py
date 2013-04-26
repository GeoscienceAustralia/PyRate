'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from itertools import product

from numpy.linalg import pinv
from numpy import nan, isnan, array, reshape, float32
from numpy import empty, zeros, dot
from numpy.testing import assert_array_equal, assert_array_almost_equal

import algorithm
from shared import Ifg
from orbital import OrbitalError, orbital_correction
from orbital import get_design_matrix, get_network_design_matrix
from orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, QUADRATIC
from tests_common import sydney5_mock_ifgs, MockIfg


# FIXME: check the offset cols as 1s in network tests

# TODO: test lstsq() here against Hua's manual method

# ORDER of work:
# fix up offsets options for QUAD DM tests
# multilooking (do ML at start when multilooking the originals)
# high level functions/workflow


class SingleDesignMatrixTests(unittest.TestCase):
	'''
	Unittests verify correctness of basic planar & quadratic design matricies	or
	DMs. This class serves two purposes, ensuring the independent method DMs are
	produced correctly. Secondly, these indivdual DMs are subsets of the larger
	DM 'grid' required for the networked orbital correction method.
	'''

	def setUp(self):
		# fake cell sizes
		self.xs = 0.75
		self.ys = 0.8
		self.ifg = Ifg("../../tests/sydney_test/obs/geo_060619-061002.unw")
		self.ifg.open()

		self.m = MockIfg(self.ifg, 3, 4)
		self.m.X_SIZE = self.xs
		self.m.Y_SIZE = self.ys


	def test_create_planar_dm(self):
		offset = False
		act = get_design_matrix(self.m, PLANAR, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 2))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset)
		assert_array_equal(act, exp)


	def test_create_planar_dm_offsets(self):
		offset = True
		act = get_design_matrix(self.m, PLANAR, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 3))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset)
		assert_array_almost_equal(act, exp)


	# tests for quadratic mode

	def test_create_quadratic_dm(self):
		offset = False
		act = get_design_matrix(self.m, QUADRATIC, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 5))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset)
		assert_array_equal(act, exp)


	def test_create_quadratic_dm_offsets(self):
		offset = True
		act = get_design_matrix(self.m, QUADRATIC, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 6))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset)
		assert_array_equal(act, exp)


	# tests for unittest_dm() assuming network method

	def test_create_planar_dm_network(self):
		# networked method planar version should not have offsets col
		ncol_exp = 2
		exp = unittest_dm(self.m, NETWORK_METHOD, PLANAR, False)
		self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
		exp2 = unittest_dm(self.m, NETWORK_METHOD, PLANAR, True)
		self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
		assert_array_equal(exp, exp2)


	def test_create_quadratic_dm_network(self):
		# quadratic version with networked method does not have offsets col
		ncol_exp = 5
		exp = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, False)
		self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
		exp2 = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, True)
		self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
		assert_array_equal(exp, exp2)



class OrbitalCorrection(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

	# FIXME: review these tests

	def test_ifg_to_vector_reshaping(self):
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
		# TODO: check correction with NaN rows

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

		ifgs = sydney5_mock_ifgs()
		for ifg in ifgs:
			ifg.X_SIZE = 90.0
			ifg.Y_SIZE = 89.5
			ifg.open()

		ifgs[0].phase_data[1, 1:3] = nan # add some NODATA

		# test both models with no offsets
		for m in [PLANAR, QUADRATIC]:
			corrections = orbital_correction(ifgs, m, INDEPENDENT_METHOD, None, False)
			test_results()

		# test both with offsets
		for m in [PLANAR, QUADRATIC]:
			corrections = orbital_correction(ifgs, m, INDEPENDENT_METHOD)
			test_results()



class ErrorTests(unittest.TestCase):
	'''Tests for the networked correction method'''

	# FIXME: review these tests

	def test_invalid_ifgs_arg(self):
		# min requirement is 1 ifg, can still subtract one epoch from the other
		self.assertRaises(OrbitalError, get_network_design_matrix, [], PLANAR, True)


	def test_invalid_degree_arg(self):
		ifgs = sydney5_mock_ifgs()
		for d in range(-5, 1):
			self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)
		for d in range(3, 6):
			self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)


	def test_invalid_method(self):
		ifgs = sydney5_mock_ifgs()
		for m in [None, 5, -1, -3, 45.8]:
			self.assertRaises(OrbitalError, orbital_correction, ifgs, PLANAR, m, None, True)


	def test_multilooked_ifgs_arg(self):
		# check a variety of bad args for network method multilooked ifgs
		ifgs = sydney5_mock_ifgs()
		args = [[None, None, None, None, None], ["X"] * 5]
		for a in args:
			args = (ifgs, PLANAR, NETWORK_METHOD, a)
			self.assertRaises(OrbitalError, orbital_correction, *args)

		# ensure failure if uneven ifgs lengths
		args = (ifgs, PLANAR, NETWORK_METHOD, ifgs[:4])
		self.assertRaises(OrbitalError, orbital_correction, *args)



class NetworkDesignMatrixTests(unittest.TestCase):
	'''Contains tests verifying creation of sparse network design matrix.'''
	# TODO: add nodata/nans to several layers for realism

	def setUp(self):
		self.ifgs = sydney5_mock_ifgs()
		self.nifgs = len(self.ifgs)
		self.nc = self.ifgs[0].num_cells
		self.date_ids = get_date_ids(self.ifgs)
		self.nepochs = len(self.date_ids)
		assert self.nepochs == 6

		for ifg in self.ifgs:
			ifg.X_SIZE = 90.0
			ifg.Y_SIZE = 89.5


	def test_planar_network_dm(self):
		ncoef = 2
		offset = False
		act = get_network_design_matrix(self.ifgs, PLANAR, offset)
		self.assertEqual(act.shape, (self.nc * self.nifgs, ncoef * self.nepochs))
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_planar_network_dm_offset(self):
		ncoef = 2 # NB: doesn't include offset col
		offset = True
		act = get_network_design_matrix(self.ifgs, PLANAR, offset)
		self.assertEqual(act.shape[0], self.nc * self.nifgs)
		self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)

		# TODO: check offsets cols thoroughly
		self.assertTrue(act[-1, -1] == 1)
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_network_design_matrix_quad(self):
		ncoef = 5
		offset = False
		act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
		self.assertEqual(act.shape, (self.nc * self.nifgs, ncoef * self.nepochs))
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_network_design_matrix_quad_offset(self):
		ncoef = 5
		offset = True
		act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
		exp = (self.nc * self.nifgs, (self.nepochs * ncoef) + self.nifgs )
		self.assertEqual(act.shape, exp)

		# TODO: check offsets cols thoroughly
		self.assertTrue(act[-1, -1] == 1)
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def check_equality(self, ncoef, dm, ifgs, offset):
		'''
		Internal test function to check subsets against network design matrix
		ncoef - base number of coefficients, without extra col for offsets
		dm - network design matrix to check the results
		ifgs - sequence of Ifg objs
		offset - boolean to include extra parameters for model offsets
		'''
		self.assertTrue(ncoef in [2,5])
		deg = PLANAR if ncoef < 5 else QUADRATIC
		np = ncoef * self.nepochs # base offset for offset col testing

		for i, ifg in enumerate(ifgs):
			exp = unittest_dm(ifg, NETWORK_METHOD, deg, offset)
			self.assertEqual(exp.shape, (ifg.num_cells, ncoef)) # subset DMs shouldn't have an offsets col

			# NB: this code is Hua Wang's slightly modified for Py
			ib1 = i * self.nc # start row for subsetting the sparse matrix
			ib2 = (i+1) * self.nc # last row of subset of sparse matrix
			jbm = self.date_ids[ifg.MASTER] * ncoef # starting col index for master
			jbs = self.date_ids[ifg.SLAVE] * ncoef # col start for slave
			assert_array_almost_equal(-exp, dm[ib1:ib2, jbm:jbm+ncoef])
			assert_array_almost_equal( exp, dm[ib1:ib2, jbs:jbs+ncoef])

			# ensure rest of row is zero
			assert_array_equal(0, dm[ib1:ib2, :jbm])
			assert_array_equal(0, dm[ib1:ib2, jbm+ncoef:jbs])

			# check offset cols
			if offset is True:
				self.assertTrue(bool((dm[ib1:ib2, i+np] == 1).all()))
				assert_array_equal(0, dm[ib1:ib2, jbs+ncoef:i+np])
				assert_array_equal(0, dm[ib1:ib2, i+np+1:])
			else:
				assert_array_equal(0, dm[ib1:ib2, jbs+ncoef:])


	def test_network_correct_planar(self):
		'''Verifies planar form of network method of correction'''

		# add some nans to the data
		self.ifgs[0].phase_data[0, :] = nan # 3 error cells
		self.ifgs[1].phase_data[2, 1:3] = nan # 2 error cells
		self.ifgs[2].phase_data[3, 2:3] = nan # 1 err
		self.ifgs[3].phase_data[1, 2] = nan # 1 err
		self.ifgs[4].phase_data[1, 1:3] = nan # 2 err
		err = sum([i.nan_count for i in self.ifgs])

		# reshape phase data of all ifgs to single vector
		data = empty(self.nifgs * self.nc, dtype=float32)
		for i, ifg in enumerate(self.ifgs):
			st = i * self.nc
			data[st:st + self.nc] = ifg.phase_data.reshape(self.nc)

		for off in [False, True]:
			dm = get_network_design_matrix(self.ifgs, PLANAR, off)[~isnan(data)]
			new_sh = ( (self.nc * self.nifgs) - err, 1)
			fd = data[~isnan(data)].reshape(new_sh)
			self.assertEqual(len(dm), len(fd))

			ncoef = 3 if off else 2
			nepochs = len(self.ifgs) + 1
			params = dot(pinv(dm, 1e-6), fd)
			self.assertEqual(params.shape, (nepochs * ncoef, 1) ) # needs to be matrix multiplied/reduced

			# calculate forward correction
			sdm = unittest_dm(ifg, NETWORK_METHOD, PLANAR) # base DM to multiply out to real data
			raise NotImplementedError("TODO: determine how many ncoefs this should be tested for in offsets=TRUE")
			self.assertEqual(sdm.shape, (12,ncoef))

			corrections = []
			for i, ifg in enumerate(self.ifgs):
				jbm = self.date_ids[ifg.MASTER] * ncoef # starting row index for master
				jbs = self.date_ids[ifg.SLAVE] * ncoef # row start for slave
				par = params[jbs:jbs + ncoef] - params[jbm:jbm + ncoef]
				corrections.append(dot(sdm, par).reshape(ifg.phase_data.shape))

			act = orbital_correction(self.ifgs, PLANAR, NETWORK_METHOD, None, off)
			assert_array_almost_equal(act, corrections, decimal=4) # TODO: fails intermittently on default decimal


	def test_network_correct_quadratic(self):
		'''Verifies quadratic form of network method of correction'''

		self.ifgs[0].phase_data[3, 0:7] = nan # add NODATA as rasters are complete

		# reshape phase data to vectors
		data = empty(self.nifgs * self.nc, dtype=float32)
		for i, ifg in enumerate(self.ifgs):
			st = i * self.nc
			data[st:st + self.nc] = ifg.phase_data.reshape(self.nc)

		dm = get_network_design_matrix(self.ifgs, QUADRATIC, False)[~isnan(data)]

		raise NotImplementedError("use dot() for matrix mult of pinv() results")
		params = pinv(dm, 1e-6) * data[~isnan(data)] # FIXME: use dot() instead

		# TODO: add fwd correction calc & testing here

		raise NotImplementedError("TODO: add forward correction implementation")

		act = orbital_correction(self.ifgs, QUADRATIC, NETWORK_METHOD, None, False)
		assert_array_almost_equal(act, params, decimal=5) # TODO: fails occasionally on default decimal=6

		# FIXME: expand test to include offsets


def unittest_dm(ifg, method, degree, offset=False):
	'''Helper/test func to create design matrix segments. Includes handling for
	making quadratic DM segments for use in network method.
	ifg - source interferogram to model design matrix on
	method - INDEPENDENT_METHOD or NETWORK_METHOD
	degree - PLANAR or QUADRATIC
	offset - True/False to include additional cols for offsets
	'''
	assert method in [INDEPENDENT_METHOD, NETWORK_METHOD]
	assert degree in [PLANAR, QUADRATIC]

	NX = ncoef = 2 if degree == PLANAR else 5
	if offset is True:
		if method == INDEPENDENT_METHOD:
			ncoef += 1
		else:
			offset = False # prevent offsets in DM sections for network method

	data = empty((ifg.num_cells, ncoef), dtype=float32)
	rows = iter(data)
	yr = xrange(ifg.FILE_LENGTH)
	xr = xrange(ifg.WIDTH)

	if degree == PLANAR:
		for y,x in product(yr, xr):
			row = rows.next()
			row[:NX] = [x * ifg.X_SIZE, y * ifg.Y_SIZE]
	else:
		for y,x in product(yr, xr):
			ys = y * ifg.Y_SIZE
			xs = x * ifg.X_SIZE
			row = rows.next()
			row[:NX] = [xs**2, ys**2, xs*ys, xs, ys]

	if offset:
		data[:, -1] = 1

	return data


def get_date_ids(ifgs):
	'''Returns unique master/slave date IDs from the given Ifgs'''
	dates = []
	for ifg in ifgs:
		dates += list(ifg.DATE12)
	return algorithm.master_slave_ids(dates)


if __name__ == "__main__":
	unittest.main()
