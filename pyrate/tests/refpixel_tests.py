'''
Collection of tests for validating PyRate's reference pixel code.
Author: Ben Davies
'''

import unittest
from numpy import nan, mean, std, isnan

from common import sydney_data_setup, MockIfg
from pyrate.refpixel import ref_pixel, RefPixelError, _step

from pyrate.config import ConfigException
from pyrate.config import REFX, REFNX, REFY, REFNY, REF_CHIP_SIZE, REF_MIN_FRAC


def default_params():
	return { REFNX : 5, REFNY : 7, REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : 3 }


class ReferencePixelInputTests(unittest.TestCase):
	'''Verifies error checking capabilities of the reference pixel function'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_data_setup()


	def test_missing_chipsize(self):
		params = default_params()
		del params[REF_CHIP_SIZE]
		self.assertRaises(ConfigException, ref_pixel, params, self.ifgs)


	def test_chipsize_valid(self):
		params = {}
		for illegal in [0, -1, -15, 1, 2, self.ifgs[0].WIDTH+1, 4, 6, 10, 20]:
			params[REF_CHIP_SIZE] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)


	def test_minimum_fraction_missing(self):
		params = default_params()
		del params[REF_MIN_FRAC]
		self.assertRaises(ConfigException, ref_pixel, params, self.ifgs)


	def test_minimum_fraction_threshold(self):
		params = default_params()
		for illegal in [-0.1, 1.1, 1.000001, -0.0000001]:
			params[REF_MIN_FRAC] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)


	def test_predefined_reference_pixel(self):
		# return reference pixel coords if already set in config
		for exp_coord in [(3, 7), (0, 0)]:
			Y, X = exp_coord
			params = { REFX : X, REFY : Y }
			act = ref_pixel(params, self.ifgs)
			self.assertEqual(exp_coord, act)


	def test_invalid_reference_pixel(self):
		# ensure refx & refy are within the grid (if not < 0)
		params = default_params()

		base = self.ifgs[0].WIDTH
		for illegal in [base + 1, base + 7]:
			params[REFX] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)

		params[REFX] = 5 # valid coord to ensure testing of REFY
		base = self.ifgs[0].FILE_LENGTH
		
		for illegal in [base + 1, base + 9]:
			params[REFY] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)


	def test_search_windows(self):
		params = default_params()

		# 45 is max # cells a width 3 sliding window can iterate over
		for illegal in [-5, -1, 0, 46, 50, 100]:
			params[REFNX] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)

		params[REFNX] = 3
		# 40 is max # cells a width 3 sliding window can iterate over
		for illegal in [-5, -1, 0, 71, 85, 100]:
			params[REFNX] = illegal
			self.assertRaises(ValueError, ref_pixel, params, self.ifgs)


	def test_missing_search_windows(self):
		params = default_params()
		del params[REFNX]
		self.assertRaises(ConfigException, ref_pixel, params, self.ifgs)

		params[REFNX] = 3 # reset
		del params[REFNY]
		self.assertRaises(ConfigException, ref_pixel, params, self.ifgs)



class ReferencePixelTests(unittest.TestCase):
	'''Tests results of the reference pixel search'''

	def setUp(self):
		self.testdir, self.ifgs = sydney_data_setup()
		self.mock_ifgs = None


	def test_empty_refxy(self):
		# no REFX|Y params in config should cause refpixel search to occur
		params = default_params()
		self.assertTrue(REFX not in params)
		self.assertTrue(REFY not in params)
		
		refpx = ref_pixel(params, self.ifgs)
		self.assertEqual(len(refpx), 2)
		self.assertNotEqual(refpx, (0,0))
		self.assertNotEqual(refpx, (-1,-1))


	def test_subzero_refxy(self):
		# subzero REFX|Y params in config should also result in refpixel search
		params = default_params()
		
		for v in [-1, -12]:		
			params[REFX] = params[REFY] = v
			refpx = ref_pixel(params, self.ifgs)
			self.assertEqual(len(refpx), 2)
			self.assertNotEqual(refpx, (0,0))
			self.assertNotEqual(refpx, (-1,-1))


	def test_all_below_threshold_exception(self):
		# test failure when no valid stacks in dataset
		params = { REFNX : 2, REFNY : 2, REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : 3 }

		# rig mock data to be below threshold
		self.mock_ifgs = [MockIfg(i, 6, 7) for i in self.ifgs]
		for m in self.mock_ifgs:
			m.phase_data[:1] = nan
			m.phase_data[1:5] = 0.1
			m.phase_data[5:] = nan

		self.assertRaises(RefPixelError, ref_pixel, params, self.mock_ifgs)


	def test_refnxy_1(self):
		# test step of 1 for refnx|y gets the reference pixel for axis centre
		params = { REFNX : 1, REFNY : 1, REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : 3 }

		self.mock_ifgs = [MockIfg(i, 47, 72) for i in self.ifgs]
		for m in self.mock_ifgs:
			m.phase_data[:1] = 0.2
			m.phase_data[1:5] = 0.1
			m.phase_data[5:] = 0.3

		exp_refpx = (36,23)
		act_refpx = ref_pixel(params, self.mock_ifgs)
		self.assertEqual(exp_refpx, act_refpx)


	def test_large_window(self):
		# 5x5 view over a 5x5 ifg with 1 window/ref pix search
		chps = 5
		params = { REFNX : 1, REFNY : 1, REF_MIN_FRAC : 0.7, REF_CHIP_SIZE : chps }
		self.mock_ifgs = [MockIfg(i, chps, chps) for i in self.ifgs]
		act_refpx = ref_pixel(params, self.mock_ifgs)
		self.assertEqual((2,2), act_refpx)


	def test_step(self):
		# test different search windows to verify x/y step calculation

		# convenience testing function
		def assert_equal(actual, expected):
			for a, e in zip(actual, expected):
				self.assertEqual(a, e)

		# start with simple corner only test
		width = 47
		radius = 2
		refnx = 2
		exp = [2, 44]
		act = _step(width, refnx, radius)
		assert_equal(act, exp)

		# test with 3 windows
		refnx = 3
		exp = [2, 23, 44]
		act = _step(width, refnx, radius)
		assert_equal(act, exp)

		# test 4 search windows
		refnx = 4
		exp = [2, 16, 30, 44]
		act = _step(width, refnx, radius)
		assert_equal(act, exp)


	def test_ref_pixel(self):
		params = default_params()
		params[REFNX] = 2 # just use the corners
		params[REFNY] = 2
		params[REF_CHIP_SIZE] = 5

		exp_refpx = (2,2) # calculated manually from _expected_ref_pixel()
		act_refpx = ref_pixel(params, self.ifgs)
		self.assertNotEqual(act_refpx, (0,0))
		self.assertEqual(act_refpx, exp_refpx)

		# Invalidate first stack and test result again
		for i in self.ifgs:
			i.phase_data[:3,:5] = nan

		exp_refpx = (2,44) # calculated manually from _expected_ref_pixel()
		act_refpx = ref_pixel(params, self.ifgs)
		self.assertEqual(act_refpx, exp_refpx)


def _expected_ref_pixel(ifgs, cs):
	'''Helper function for finding reference pixel when refnx/y=2'''

	# calculate expected data
	data = [i.phase_data for i in ifgs] # len 17 list of ndarrays
	ul = [ i[:cs,:cs] for i in data] # upper left corner stack
	ur = [ i[:cs,-cs:] for i in data]
	ll = [ i[-cs:,:cs] for i in data]
	lr = [ i[-cs:,-cs:] for i in data]

	ulm = mean([std(i[~isnan(i)]) for i in ul]) # mean std of all the layers
	urm = mean([std(i[~isnan(i)]) for i in ur])
	llm = mean([std(i[~isnan(i)]) for i in ll])
	lrm = mean([std(i[~isnan(i)]) for i in lr])
	assert isnan([ulm, urm, llm, lrm]).any() == False

	# coords of the smallest mean is the result
	mn = [ulm, urm, llm, lrm]
	print mn, min(mn), mn.index(min(mn))

