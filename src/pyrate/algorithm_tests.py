
from os.path import join
import glob, unittest, datetime
from math import pi, cos, sin, radians
from itertools import product

import numpy
from numpy import array, nan, reshape, squeeze, ndarray, float32
from numpy.testing import assert_array_almost_equal

import algorithm
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


class InitialModelTests(unittest.TestCase):

	def test_TODO(self):
		# TODO: fake an RSC file with coords
		# TODO: fake a ones(shape)  # could also make a ramp etc
		# data is single band of DISPLACEMENT
		raise NotImplementedError



class ReferencePixelTests(unittest.TestCase):
	'''TODO'''
	# TODO: warnings vs exceptions?
	# TODO: test missing items in params

	def setUp(self):
		self.testdir, self.ifgs = sydney_test_setup()


	def default_params(self):
		return { REFNX : 50,
							REFNY : 50,
							REF_MIN_FRAC : 0.7,
							REF_CHIP_SIZE : 3,  }


	def test_missing_chipsize(self):
		params = self.default_params()
		del params[REF_CHIP_SIZE]
		self.assertRaises(ConfigException, algorithm.ref_pixel, params, self.ifgs)


	def test_chipsize_valid(self):
		params = {}
		for illegal in [0, -1, -15, 1, 2, self.ifgs[0].WIDTH+1, 4, 6, 10, 20]:
			params[REF_CHIP_SIZE] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_ref_window_size(self):
		params = { REFNX : None, REFNY : None }
		illegal_values = [0, -1, -15]

		for i in illegal_values:
			params[REF_CHIP_SIZE] = i
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)

		params[REFNX] = 21 # dummy but valid, ensure Y axis tests happen
		for i in illegal_values:
			params[REFNY] = i
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_minimum_fraction_threshold(self):
		params = self.default_params()
		for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
			params[REF_MIN_FRAC] = illegal
			self.assertRaises(ValueError, algorithm.ref_pixel, params, self.ifgs)


	def test_predefined_coords(self):
		# TODO: refx, refy are within the grid (if not 0)
		raise NotImplementedError



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

		self.mock_ifgs = [MockIfg(i) for i in self.ifgs]
		for m in self.mock_ifgs[num_coherent:]:
			m.phase_data[:] = nan
		test_result()

		# fill in more nans leaving only one ifg
		for m in self.mock_ifgs[1:num_coherent]:
			m.phase_data[:] = nan
		num_coherent = 1
		test_result()


	def test_all_nan_pixel_stack(self):
		self.mock_ifgs = [MockIfg(i) for i in self.ifgs]
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
	'''Mock 1x1 Ifg for fine grain testing of values'''

	def __init__(self, src_ifg):
		self.MASTER = src_ifg.MASTER
		self.SLAVE = src_ifg.SLAVE
		self.DATE12 = src_ifg.DATE12

		self.FILE_LENGTH = 1
		self.WIDTH = 1
		self.phase_data = ndarray((self.FILE_LENGTH, self.WIDTH), dtype=float32)
		self.nan_fraction = src_ifg.nan_fraction # use existing overall nan fraction
