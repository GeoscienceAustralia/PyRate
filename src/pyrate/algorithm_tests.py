
from os.path import join
import glob, unittest, datetime
from math import pi, cos, sin, radians
from itertools import product

import numpy
from numpy import array, reshape, squeeze, ndarray
from numpy import isnan, nan, float32, std, mean
from numpy.testing import assert_array_almost_equal

import algorithm
from algorithm import RefPixelError
from shared import Ifg
from config import ConfigException
from config import REFX, REFNX, REFY, REFNY, REF_CHIP_SIZE, REF_MIN_FRAC


class AlgorithmTests(unittest.TestCase):

	def test_wavelength_conversion(self):
		# ROIPAC is in metres, verify conversion to mm
		xs, ys = 5,7
		data = (numpy.arange(xs * ys) - 1.7) * 0.1 # fake a range of values
		data = numpy.where(data == 0, numpy.nan, data)
		wavelen = 0.0562356424
		exp = (data * wavelen * 1000) / (4 * pi)
		act = algorithm.wavelength_to_mm(data, wavelen)
		numpy.testing.assert_array_almost_equal(exp, act)


	def test_nan_fraction(self):
		data = numpy.ones((4,3))
		data[1] = nan
		exp = 3 / float(12)
		act = algorithm.nan_fraction(data)
		self.assertEqual(act, exp)


	def test_unit_vector(self):
		incidence = [radians(x) for x in (34.3, 39.3, 29.3, 22.8) ]
		azimuth = [radians(x) for x in (77.8, 77.9, 80.0, 80.3)  ]

		vert, ns, ew = [], [], []
		for i, a in zip(incidence, azimuth):
			vert.append(cos(i))
			ns.append(sin(i) * sin(a))
			ew.append(sin(i) * cos(a))

		sh = (2,2)
		unitv = [array(ew), array(ns), array(vert)]
		unitv = [a.reshape(sh) for a in unitv]

		# TODO: assumes rad input for now
		act = algorithm.unit_vector(reshape(incidence, sh), reshape(azimuth, sh))
		for a,e in zip(act, unitv):
			assert_array_almost_equal(squeeze(a), e)

# TODO: InitialModelTests
#class InitialModelTests(unittest.TestCase):

#	def test_initial_model(self):
		# TODO: fake an RSC file with coords
		# TODO: fake a ones(shape)  # could also make a ramp etc
		# data is single band of DISPLACEMENT
		#raise NotImplementedError


def default_params():
	return { REFNX : 5, REFNY : 7,
					REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : 3 }


class ReferencePixelInputTests(unittest.TestCase):
	'''Verifies error checking capabilities of the reference pixel function'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()


	def test_missing_chipsize(self):
		params = default_params()
		del params[REF_CHIP_SIZE]
		self.assertRaises(ConfigException, algorithm.ref_pixel, params, self.ifgs)


	def test_chipsize_valid(self):
		params = {}
		for illegal in [0, -1, -15, 1, 2, self.ifgs[0].WIDTH+1, 4, 6, 10, 20]:
			params[REF_CHIP_SIZE] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_minimum_fraction_missing(self):
		params = default_params()
		del params[REF_MIN_FRAC]
		self.assertRaises(ConfigException, algorithm.ref_pixel, params, self.ifgs)


	def test_minimum_fraction_threshold(self):
		params = default_params()
		for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
			params[REF_MIN_FRAC] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_predefined_reference_pixel(self):
		# return reference pixel coords if already set in config
		exp_coord = 3, 7
		Y, X = exp_coord
		params = { REFX : X, REFY : Y }
		act = algorithm.ref_pixel(params, self.ifgs)


	def test_invalid_reference_pixel(self):
		# ensure refx & refy are within the grid (if not 0)
		params = default_params()

		for illegal in [-5, -1, self.ifgs[0].WIDTH+1]:
			params[REFX] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)

		params[REFX] = 5 # valid coord to ensure testing of REFY
		for illegal in [-5, -1, self.ifgs[0].FILE_LENGTH+1]:
			params[REFY] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	# TODO: determine action for step of 1? can't be done in corner
	def test_search_windows(self):
		params = default_params()
		for illegal in [-5, -1, 0, 46, 50, 100]: # 45 is max # cells a width 3 sliding window can iterate over
			params[REFNX] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)

		params[REFNX] = 3
		for illegal in [-5, -1, 0, 71, 85, 100]: # 40 is max # cells a width 3 sliding window can iterate over
			params[REFNX] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_missing_search_windows(self):
		params = default_params()
		del params[REFNX]
		self.assertRaises(ConfigException, algorithm.ref_pixel, params, self.ifgs)

		params[REFNX] = 3 # reset
		del params[REFNY]
		self.assertRaises(ConfigException, algorithm.ref_pixel, params, self.ifgs)



class ReferencePixelTests(unittest.TestCase):
	'''Tests results of the reference pixel search'''

	# TODO: test case: 5x5 view over a 5x5 ifg with 1 window/ref pix search
	# TODO: test result where one window is < thresh

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()


	def test_all_below_threshold_exception(self):
		# test failure when no valid stacks in dataset
		params = { REFNX : 2, REFNY : 2, REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : 3 }

		# rig mock data to be below threshold
		self.mock_ifgs = [MockIfg(i, 6, 7) for i in self.ifgs]
		for m in self.mock_ifgs:
			m.phase_data[:1] = nan
			m.phase_data[1:5] = 0.1
			m.phase_data[5:] = nan

		self.assertRaises(RefPixelError, algorithm.ref_pixel, params, self.mock_ifgs)


	def test_step(self):
		# test different search windows to verify x/y step calculation

		# convenience testing function
		def test_equal(actual, expected):
			for a, e in zip(actual, expected):
				self.assertEqual(a, e)

		# start with simple corner only test
		width = 47
		radius = 2
		refnx = 2
		exp = [2, 44]
		act = algorithm._step(width, refnx, radius)
		test_equal(act, exp)

		# test with 3 windows
		refnx = 3
		exp = [2, 23, 44]
		act = algorithm._step(width, refnx, radius)
		test_equal(act, exp)

		# test 4 search windows
		refnx = 4
		exp = [2, 16, 30, 44]
		act = algorithm._step(width, refnx, radius)
		test_equal(act, exp)


	# TODO: try this data but NaN out the first corner
	def test_ref_pixel(self):
		params = default_params()
		params[REFNX] = 2 # just use the corners
		params[REFNY] = 2
		params[REF_CHIP_SIZE] = 5

		exp_refpx = (2,2) # calculated manually from _expected_ref_pixel()
		act_refpx = algorithm.ref_pixel(params, self.ifgs)
		self.assertNotEqual(act_refpx, (0,0))
		self.assertEqual(act_refpx, exp_refpx)


def _expected_ref_pixel(ifgs, cs):
	'''Helper function for finding reference pixel when refnx/y=2'''

	# calculate expected data
	data = [i.phase_data for i in ifgs] # len 17 list of ndarrays
	ul = [ i[:cs,:cs] for i in data] # upper left corner stack
	ur = [ i[:cs,-cs:] for i in data]
	ll = [ i[-cs:,:cs] for i in data]
	lr = [ i[-cs:,-cs:] for i in data]

	ulm = mean([std(i) for i in ul]) # mean std of all the layers
	urm = mean([std(i) for i in ur])
	llm = mean([std(i) for i in ll])
	lrm = mean([std(i) for i in lr])
	assert isnan([ulm, urm, llm, lrm]).any() == False

	# coords of the smallest mean is the result
	mn = [ulm, urm, llm, lrm]
	print mn, min(mn)



class EpochListTests(unittest.TestCase):

	def test_get_epochs(self):
		def str2date(s):
			segs = s[:4], s[4:6], s[6:] # year, month, day
			return datetime.date(*[int(s) for s in segs])

		raw_date = ['20060619', '20060828', '20061002', '20061106', '20061211',
							'20070115', '20070219', '20070326', '20070430', '20070604',
							'20070709', '20070813', '20070917']
		exp_dates = [str2date(d) for d in raw_date]
		exp_repeat = [1,1,3,3,4,3,3,3,3,3,3,2,2]

		exp_spans = [0, 0.1916, 0.2875, 0.3833, 0.4791, 0.5749, 0.6708, 0.7666,
							0.8624, 0.9582, 1.0541, 1.1499, 1.2457]

		# test against Hua's results
		base = "../../tests/sydney_test/obs"
		paths = join(base, "ifms_17")
		with open(paths) as f:
			ifgs = [Ifg(join(base, path)) for path in f.readlines()]
			epochs = algorithm.get_epochs(ifgs)

		self.assertTrue((exp_dates == epochs.dates).all())
		self.assertTrue((exp_repeat == epochs.repeat).all())
		assert_array_almost_equal(exp_spans, epochs.spans, decimal=4)



class MSTTests(unittest.TestCase):
	'''Basic verification of minimum spanning tree (MST) functionality.'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()
		self.epochs = algorithm.get_epochs(self.ifgs)


	def test_mst_matrix(self):
		# Verifies mst matrix function returns an array with of dicts in each cell
		# assumes pygraph is correct from its unit tests
		res = algorithm.mst_matrix(self.ifgs, self.epochs)
		ys, xs = res.shape
		for y, x in product(xrange(ys), xrange(xs)):
			r = res[y,x]
			self.assertTrue(hasattr(r, "keys"))
			self.assertTrue(len(r) <= len(self.epochs.dates))


	def test_partial_nan_pixel_stack(self):
		# Ensure a limited # of cells gives in a smaller node tree

		num_coherent = 3

		def test_result():
			res = algorithm.mst_matrix(self.mock_ifgs, self.epochs)
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

		res = algorithm.mst_matrix(self.mock_ifgs, self.epochs)
		exp = [nan]

		shape = (self.mock_ifgs[0].FILE_LENGTH, self.mock_ifgs[0].WIDTH)
		self.assertTrue(res.shape == shape)
		self.assertEqual(exp, res)


def sydney_test_setup():
	testdir = "../../tests/sydney_test/obs"
	datafiles = glob.glob( join(testdir, "*.unw") )
	ifgs = [Ifg(i) for i in datafiles]
	for i in ifgs:
		i.open()

	return testdir, ifgs



class MockIfg(object):
	'''Mock Ifg for detailed testing'''

	def __init__(self, src_ifg, xsize=None, ysize=None):
		self.MASTER = src_ifg.MASTER
		self.SLAVE = src_ifg.SLAVE
		self.DATE12 = src_ifg.DATE12

		self.FILE_LENGTH = ysize
		self.WIDTH = xsize
		self.phase_data = ndarray((self.FILE_LENGTH, self.WIDTH), dtype=float32)
		self.nan_fraction = src_ifg.nan_fraction # use existing overall nan fraction
