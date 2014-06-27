'''
Collection of tests for orbital correction.

Created on 31/3/13
@author: Ben Davies
'''


import unittest
from os.path import join
from itertools import product

from numpy.linalg import pinv, inv
from numpy import nan, isnan, array, reshape, median
from numpy import empty, dot, concatenate, float32
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyrate import algorithm
from pyrate.shared import Ifg
from pyrate.orbital import OrbitalError, orbital_correction
from pyrate.orbital import get_design_matrix, get_network_design_matrix, get_num_params
from pyrate.orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, QUADRATIC, PART_CUBIC
from common import sydney5_mock_ifgs, MockIfg, SYD_TEST_TIF
from scipy.linalg import lstsq

DEG_LOOKUP = { 2 : PLANAR,
				5 : QUADRATIC,
				6 : PART_CUBIC, }

NUM_COEF_LOOKUP = { PLANAR : 2,
					QUADRATIC : 5,
					PART_CUBIC : 6, }


class SingleDesignMatrixTests(unittest.TestCase):
	'''
	Tests to verify correctness of basic planar & quadratic design matricies or
	DMs. This class serves two purposes, ensuring the independent method DMs are
	produced correctly. Secondly, these indivdual DMs are subsets of the larger
	DM 'grid' required for the networked orbital correction method.
	'''

	def setUp(self):
		# faked cell sizes
		self.xs = 0.75
		self.ys = 0.8
		self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_060619-061002.tif'))
		self.ifg.open()

		self.m = MockIfg(self.ifg, 3, 4)
		self.m.x_size = self.xs
		self.m.y_size = self.ys


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

	# tests for partial cubic model

	def test_create_partcubic_dm(self):
		offset = False
		act = get_design_matrix(self.m, PART_CUBIC, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 6))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset)
		assert_array_equal(act, exp)


	def test_create_partcubic_dm_offsets(self):
		offset = True
		act = get_design_matrix(self.m, PART_CUBIC, offset)
		self.assertEqual(act.shape, (self.m.num_cells, 7))
		exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset)
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

	def test_create_partcubic_dm_network(self):
		# partial cubic version with networked method does not have offsets col
		ncol_exp = 6
		exp = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, False)
		self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
		exp2 = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, True)
		self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
		assert_array_equal(exp, exp2)


class IndependentCorrectionTests(unittest.TestCase):
	'''Test cases for the orbital correction component of PyRate.'''

	def setUp(self):
		self.ifgs = sydney5_mock_ifgs()
		_add_nodata(self.ifgs)

		for ifg in self.ifgs:
			ifg.x_size = 90.0
			ifg.y_size = 89.5
			ifg.open()


	def alt_orbital_correction(self, ifg, deg, offset):
		# almost 1:1 copy of MATLAB version
		data = ifg.phase_data.reshape(ifg.num_cells)
		dm = get_design_matrix(ifg, deg, offset)[~isnan(data)]
		fd = data[~isnan(data)].reshape((dm.shape[0], 1))

		dmt = dm.T
		invNbb = inv(dmt.dot(dm))
		params = invNbb.dot(dmt.dot(fd))
		alt_params = lstsq(dm, fd)[0]
		assert_array_almost_equal(params, alt_params, decimal=2) # FIXME: precision

		dm2 = get_design_matrix(ifg, deg, offset)
		fwd_correction = reshape(dot(dm2, params), ifg.phase_data.shape)
		return ifg.phase_data - fwd_correction


	def check_correction(self, degree, method, offset):
		orig = array([c.phase_data.copy() for c in self.ifgs])
		exp = [self.alt_orbital_correction(i, degree, offset) for i in self.ifgs]
		orbital_correction(self.ifgs, degree, method, None, offset)
		corrected = array([c.phase_data for c in self.ifgs])

		self.assertFalse((orig == corrected).all())
		self.check_results(self.ifgs, orig) # test shape, data is non zero

		# FIXME: is decimal=2 close enough?
		for i, (e, a) in enumerate(zip(exp, corrected)):
			assert_array_almost_equal(e, a, decimal=2)


	def check_results(self, ifgs, corrections):
		'''Helper method for result verification'''
		for i, c in zip(ifgs, corrections):
			ys, xs = c.shape
			self.assertEqual(i.nrows, ys)
			self.assertEqual(i.ncols, xs)

			# ensure there is real data
			self.assertFalse(isnan(i.phase_data).all())
			self.assertFalse(isnan(c).all())
			self.assertTrue(c.ptp() != 0) # ensure range of values in grid


	def test_independent_correction_planar(self):
		self.check_correction(PLANAR, INDEPENDENT_METHOD, False)


	def test_independent_correction_planar_offsets(self):
		self.check_correction(PLANAR, INDEPENDENT_METHOD, True)


	def test_independent_correction_quadratic(self):
		self.check_correction(QUADRATIC, INDEPENDENT_METHOD, False)


	def test_independent_correction_quadratic_offsets(self):
		self.check_correction(QUADRATIC, INDEPENDENT_METHOD, True)


	def test_independent_correction_partcubic(self):
		self.check_correction(PART_CUBIC, INDEPENDENT_METHOD, False)


	def test_independent_correction_partcubic_offsets(self):
		self.check_correction(PART_CUBIC, INDEPENDENT_METHOD, True)

class ErrorTests(unittest.TestCase):
	'''Tests for the networked correction method'''

	def test_invalid_ifgs_arg(self):
		# min requirement is 1 ifg, can still subtract one epoch from the other
		self.assertRaises(OrbitalError, get_network_design_matrix, [], PLANAR, True)


	def test_invalid_degree_arg(self):
		# test failure of a few different args for 'degree'
		ifgs = sydney5_mock_ifgs()
		for d in range(-5, 1):
			self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)
		for d in range(4, 7):
			self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)


	def test_invalid_method(self):
		# test failure of a few different args for 'method'
		ifgs = sydney5_mock_ifgs()
		for m in [None, 5, -1, -3, 45.8]:
			self.assertRaises(OrbitalError, orbital_correction, ifgs, PLANAR, m, None)


	def test_multilooked_ifgs_arg(self):
		# check some bad args for network method with multilooked ifgs
		ifgs = sydney5_mock_ifgs()
		args = [[None, None, None, None, None], ["X"] * 5]
		for a in args:
			args = (ifgs, PLANAR, NETWORK_METHOD, a)
			self.assertRaises(OrbitalError, orbital_correction, *args)

		# ensure failure if # ifgs doesn't match # mlooked ifgs
		args = (ifgs, PLANAR, NETWORK_METHOD, ifgs[:4])
		self.assertRaises(OrbitalError, orbital_correction, *args)



class NetworkDesignMatrixTests(unittest.TestCase):
	'''Contains tests verifying creation of sparse network design matrix.'''

	def setUp(self):
		self.ifgs = sydney5_mock_ifgs()
		_add_nodata(self.ifgs)
		self.nifgs = len(self.ifgs)
		self.ncells = self.ifgs[0].num_cells
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
		self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_planar_network_dm_offset(self):
		ncoef = 2 # NB: doesn't include offset col
		offset = True
		act = get_network_design_matrix(self.ifgs, PLANAR, offset)
		self.assertEqual(act.shape[0], self.ncells * self.nifgs)
		self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_quadratic_network_dm(self):
		ncoef = 5
		offset = False
		act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
		self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_quadratic_network_dm_offset(self):
		ncoef = 5
		offset = True
		act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
		self.assertEqual(act.shape[0], self.ncells * self.nifgs)
		self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_partcubic_network_dm(self):
		ncoef = 6
		offset = False
		act = get_network_design_matrix(self.ifgs, PART_CUBIC, offset)
		self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
		self.assertNotEqual(act.ptp(), 0)
		self.check_equality(ncoef, act, self.ifgs, offset)


	def test_partcubic_network_dm_offset(self):
		ncoef = 6
		offset = True
		act = get_network_design_matrix(self.ifgs, PART_CUBIC, offset)
		self.assertEqual(act.shape[0], self.ncells * self.nifgs)
		self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
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
		deg = DEG_LOOKUP[ncoef]
		np = ncoef * self.nepochs # index of 1st offset col

		for i, ifg in enumerate(ifgs):
			exp = unittest_dm(ifg, NETWORK_METHOD, deg, offset)
			self.assertEqual(exp.shape, (ifg.num_cells, ncoef))

			# NB: this is Hua Wang's MATLAB code slightly modified for Py
			ib1, ib2 = [x * self.ncells for x in (i, i+1)] # row start/end
			jbm = ncoef * self.date_ids[ifg.master] # starting col index for master
			jbs = ncoef * self.date_ids[ifg.slave] # col start for slave
			assert_array_almost_equal(-exp, dm[ib1:ib2, jbm:jbm+ncoef])
			assert_array_almost_equal( exp, dm[ib1:ib2, jbs:jbs+ncoef])

			# ensure remaining rows/cols are zero for this ifg NOT inc offsets
			assert_array_equal(0, dm[ib1:ib2, :jbm]) # all cols leading up to master
			assert_array_equal(0, dm[ib1:ib2, jbm + ncoef:jbs]) # cols btwn mas/slv
			assert_array_equal(0, dm[ib1:ib2, jbs + ncoef:np]) # to end of non offsets

			# check offset cols for 1s and 0s
			if offset is True:
				ip1 = i + np # offset column index
				assert_array_equal(1, dm[ib1:ib2, ip1])
				assert_array_equal(0, dm[ib1:ib2, np:ip1]) # cols before offset col
				assert_array_equal(0, dm[ib1:ib2, ip1 + 1:]) # cols after offset col



class NetworkCorrectionTests(unittest.TestCase):
	'''Tests to verify Ifg correction using network method'''

	def setUp(self):
		# fake some real ifg data by adding nans
		self.ifgs = sydney5_mock_ifgs()
		_add_nodata(self.ifgs)
		self.err = sum([i.nan_count for i in self.ifgs])

		for ifg in self.ifgs:
			# use different sizes to differentiate axes results
			ifg.X_SIZE = 90.0
			ifg.Y_SIZE = 89.5

		# precalc other useful vars
		self.tol = 1e-5
		self.nc_tol = 1e-6
		self.nifgs = len(self.ifgs)
		self.ncells = self.ifgs[0].num_cells
		self.date_ids = get_date_ids(self.ifgs)
		self.nepochs = len(self.date_ids)
		assert self.nepochs == 6


	def test_offset_inversion(self):
		'Ensure pinv(DM)*obs gives equal results given constant change to fd'

		def get_orbital_params():
			'''Helper func: returns pseudo-inverse of the DM'''
			data = concatenate([i.phase_data.reshape(self.ncells) for i in self.ifgs])
			dm = get_network_design_matrix(self.ifgs, PLANAR, True)[~isnan(data)]
			fd = data[~isnan(data)].reshape((dm.shape[0], 1))
			return dot(pinv(dm, self.tol), fd)

		params0 = get_orbital_params()

		# include a constant change to the observed values (fd)
		for value in [5.2, -23.5, 0]:
			for i in self.ifgs: # change ifgs in place
				i.phase_data += value
				self.assertTrue(isnan(i.phase_data).any())

			params = get_orbital_params()
			diff = params - params0
			self.assertTrue( (diff[:-self.nifgs] < self.tol).all() )
			assert_array_almost_equal(diff[-self.nifgs:], value, decimal=5)

			# reset back to orig data
			for i in self.ifgs:
				i.phase_data -= value


	# These 4 functions test full size data for orbital correction. The options
	# are separated as the ifg.phase_data arrays are modified in place, allowing
	# setUp() renew the phase data between tests.

	def test_network_correction_planar(self):
		self.network_correction(PLANAR, False)

	def test_network_correction_planar_offset(self):
		self.network_correction(PLANAR, True)

	def test_network_correction_quadratic(self):
		self.network_correction(QUADRATIC, False)

	def test_network_correction_quadratic_offset(self):
		self.network_correction(QUADRATIC, True)

	def test_network_correction_partcubic(self):
		self.network_correction(PART_CUBIC, False)

	def test_network_correction_partcubic_offset(self):
		self.network_correction(PART_CUBIC, True)


	def network_correction(self, deg, off):
		'''
		Compares results of orbital_correction() to alternate implementation.
		deg - PLANAR, QUADRATIC or PART_CUBIC
		off - True/False to calculate correction with offsets
		'''
		data = concatenate([i.phase_data.reshape(self.ncells) for i in self.ifgs])

		dm = get_network_design_matrix(self.ifgs, deg, off)[~isnan(data)]
		fd = data[~isnan(data)].reshape((dm.shape[0], 1))
		params = dot(pinv(dm, self.nc_tol), fd)
		self.assertEqual(len(dm), len(fd))
		self.assertEqual(params.shape, (dm.shape[1], 1) )

		# calculate forward correction
		sdm = unittest_dm(self.ifgs[0], NETWORK_METHOD, deg)
		ncoef = get_num_params(deg, offset=False) # NB: ignore offsets for network method
		self.assertEqual(sdm.shape, (self.ncells, ncoef) )
		orbs = self._get_corrections(self.ifgs, sdm, params, ncoef, off)

		# tricky: get expected before orbital_correction() modifies ifg phase
		exp = [i.phase_data - orb for i, orb in zip(self.ifgs, orbs)]
		orbital_correction(self.ifgs, deg, NETWORK_METHOD, None, off)
		act = [i.phase_data for i in self.ifgs]
		assert_array_almost_equal(act, exp, decimal=5)


	def _get_corrections(self, ifgs, dm, params, ncoef, offsets):
		'''
		Convenience func returns model converted to data points.
		dm: design matrix (do not filter/remove nan cells)
		params: model parameters array from pinv() * dm
		ncoef: number of model coefficients (2 planar, 5 quadratic)
		offsets: True/False to calculate correction with offsets
		'''
		corrections = []
		for ifg in ifgs:
			jbm = self.date_ids[ifg.master] * ncoef # starting row index for master
			jbs = self.date_ids[ifg.slave] * ncoef # row start for slave
			par = params[jbs:jbs + ncoef] - params[jbm:jbm + ncoef]

			# estimate orbital correction effects
			# corresponds to "fullorb = B*parm + offset" in orbfwd.m
			cor = dm.dot(par).reshape(ifg.phase_data.shape)

			if offsets:
				cor += self._estimate_offset(ifg, cor)

			corrections.append(cor)
		return corrections

	def _estimate_offset(self, ifg, mod):
		'''
		Returns estimated median offset
		ifgs - interferograms for use in estimation
		mods - basic orbital correction arrays (surfaces)
		'''
		off = (ifg.phase_data - mod).reshape(ifg.num_cells)
		return median(off[~isnan(off)])


	# These 4 functions test multilooked data for orbital correction. The options
	# are separated as the ifg.phase_data arrays are modified in place, allowing
	# setUp() renew the phase data between tests.

	def test_mlooked_network_correction_planar(self):
		self.multilooked_network_correction(PLANAR, False)


	def test_mlooked_network_correction_planar_offset(self):
		self.multilooked_network_correction(PLANAR, True)


	def test_mlooked_network_correction_quadratic(self):
		self.multilooked_network_correction(QUADRATIC, False)


	def test_mlooked_network_correction_quadratic_offset(self):
		self.multilooked_network_correction(QUADRATIC, True)

	# TODO: refactor into network_correction() as they are very similar
	def multilooked_network_correction(self, deg, off):
		'''
		Compares results of orbital_correction() to alternate implementation.
		deg - PLANAR or QUADRATIC
		off - True/False to calculate correction with offsets
		'''
		# treat default sydney mock data as multilooked
		xs, ys = 6, 8
		full = sydney5_mock_ifgs(xs, ys) # 2x as much data
		mldata = concatenate([i.phase_data.reshape(self.ncells) for i in self.ifgs])

		# calculate params from dummy mlooked
		dm = get_network_design_matrix(self.ifgs, deg, off)[~isnan(mldata)]
		fd = mldata[~isnan(mldata)].reshape((dm.shape[0], 1))
		params = dot(pinv(dm, self.nc_tol), fd)

		# scale params to full size orbital error matrix
		sdm = unittest_dm(full[0], NETWORK_METHOD, deg)
		ncoef = 2 if deg == PLANAR else 5
		self.assertEqual(sdm.shape, (full[0].num_cells, ncoef))
		orbs = self._get_corrections(full, sdm, params, ncoef, off)

		# tricky: get expected result before correction (which modifies phase_data)
		exp = [i.phase_data - orb for i, orb in zip(full, orbs)]
		orbital_correction(full, deg, NETWORK_METHOD, self.ifgs, off)
		act = [i.phase_data for i in full]
		assert_array_almost_equal(act, exp, decimal=5)


def unittest_dm(ifg, method, degree, offset=False, scale=100.0):
	'''Helper/test func to create design matrix segments. Includes handling for
	making quadratic DM segments for use in network method.
	ifg - source interferogram to model design matrix on
	method - INDEPENDENT_METHOD or NETWORK_METHOD
	degree - PLANAR, QUADRATIC or PART_CUBIC
	offset - True/False to include additional cols for offsets
	'''
	assert method in [INDEPENDENT_METHOD, NETWORK_METHOD]

	xlen = ncoef = NUM_COEF_LOOKUP[degree]
	if offset and method == INDEPENDENT_METHOD:
		ncoef += 1
	else:
		offset = False # prevent offsets in DM sections for network method

	# NB: do NOT use meshgrid as it copies the production implementation
	data = empty((ifg.num_cells, ncoef), dtype=float32)
	rows = iter(data)
	yr = xrange(1, ifg.nrows+1) # simulate meshgrid starting from 1
	xr = xrange(1, ifg.ncols+1)

	xsz, ysz = [i/scale for i in [ifg.x_size, ifg.y_size]]

	if degree == PLANAR:
		for y,x in product(yr, xr):
			row = rows.next()
			row[:xlen] = [x * xsz, y * ysz]
	elif degree ==  QUADRATIC:
		for y,x in product(yr, xr):
			ys = y * ysz
			xs = x * xsz
			row = rows.next()
			row[:xlen] = [xs**2, ys**2, xs*ys, xs, ys]
	else:
		for y,x in product(yr, xr):
			ys = y * ysz
			xs = x * xsz
			row = rows.next()
			row[:xlen] = [xs*ys**2, xs**2, ys**2, xs*ys, xs, ys]

	if offset:
		data[:, -1] = 1

	return data


def get_date_ids(ifgs):
	'''Returns unique master/slave date IDs from the given Ifgs'''
	dates = []
	for ifg in ifgs:
		dates += [ifg.master, ifg.slave]
	return algorithm.master_slave_ids(dates)


def _add_nodata(ifgs):
	'''Adds some NODATA/nan cells to the sydney mock ifgs'''
	ifgs[0].phase_data[0, :] = nan # 3 error cells
	ifgs[1].phase_data[2, 1:3] = nan # 2 error cells
	ifgs[2].phase_data[3, 2:3] = nan # 1 err
	ifgs[3].phase_data[1, 2] = nan # 1 err
	ifgs[4].phase_data[1, 1:3] = nan # 2 err


if __name__ == "__main__":
	unittest.main()
