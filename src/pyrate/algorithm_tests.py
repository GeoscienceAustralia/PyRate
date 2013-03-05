'''
Collection of algorithms used in PyRate.
Author: Ben Davies
'''

from os.path import join
from datetime import date
from unittest import TestCase
from math import pi, cos, sin, radians

import numpy
from numpy import array, reshape, squeeze
from numpy.testing import assert_array_almost_equal, assert_allclose

import algorithm
from shared import Ifg
from tests_common import sydney5_mock_ifgs



class AlgorithmTests(TestCase):
	'''Misc unittests for functions in the algorithm module.'''

	def test_wavelength_conversion(self):
		# ROIPAC is in metres, verify conversion to mm
		xs, ys = 5, 7
		data = (numpy.arange(xs * ys) - 1.7) * 0.1 # fake a range of values
		data = numpy.where(data == 0, numpy.nan, data)
		wavelen = 0.0562356424
		exp = (data * wavelen * 1000) / (4 * pi)
		act = algorithm.wavelength_to_mm(data, wavelen)
		assert_allclose(exp, act)


	def test_unit_vector(self):
		incidence = [radians(x) for x in (34.3, 39.3, 29.3, 22.8)]
		azimuth = [radians(x) for x in (77.8, 77.9, 80.0, 80.3)]

		vert, ns, ew = [], [], []
		for i, a in zip(incidence, azimuth):
			vert.append(cos(i))
			ns.append(sin(i) * sin(a))
			ew.append(sin(i) * cos(a))

		sh = (2, 2)
		unitv = [array(ew), array(ns), array(vert)]
		unitv = [a.reshape(sh) for a in unitv]

		# NB: assumes radian inputs
		act = algorithm.unit_vector(reshape(incidence, sh), reshape(azimuth, sh))
		for a, e in zip(act, unitv):
			assert_array_almost_equal(squeeze(a), e)


	def test_master_slave_ids(self):
		d0 = date(2006, 6, 19)
		d1 = date(2006, 8, 28)
		d2 = date(2006, 10, 02)
		d3 = date(2006, 11, 06)
		exp = { d0: 0, d1: 1, d2: 2, d3: 3}

		# test unordered and with duplicates
		self.assertEqual(exp, algorithm.master_slave_ids([d3, d0, d2, d1]))
		self.assertEqual(exp, algorithm.master_slave_ids([d3, d0, d2, d1, d3, d0]))


	def test_ifg_date_lookup(self):
		# check reverse lookup of ifg given a master and slave date tuple
		ifgs = sydney5_mock_ifgs()
		date_pair = (date(2006, 8, 28), date(2006, 12, 11))
		i = algorithm.ifg_date_lookup(ifgs, date_pair)
		self.assertEqual(ifgs[0], i)

		# test with reversed date tuple, should reorder it according to age
		date_pair = (date(2006, 12, 11), date(2006, 11, 6))
		i = algorithm.ifg_date_lookup(ifgs, date_pair)
		self.assertEqual(ifgs[1], i)


	def test_ifg_date_lookup_failure(self):
		# error when lookup cannot find an ifg given a date pair
		ifgs = sydney5_mock_ifgs()
		date_pair = (date(2006, 12, 11), date(2007, 3, 26))
		self.assertRaises(ValueError, algorithm.ifg_date_lookup, ifgs, date_pair)

	# TODO: test bad inputs for date_pair?



# TODO: InitialModelTests
#class InitialModelTests(unittest.TestCase):

#	def test_initial_model(self):
		# TODO: fake an RSC file with coords
		# TODO: fake a ones(shape)  # could also make a ramp etc
		# data is single band of DISPLACEMENT
		#raise NotImplementedError



class EpochListTests(TestCase):
	'''Unittests for the EpochList class.'''

	def test_get_epochs(self):
		def str2date(s):
			segs = s[:4], s[4:6], s[6:] # year, month, day
			return date(*[int(s) for s in segs])

		raw_date = ['20060619', '20060828', '20061002', '20061106', '20061211',
							'20070115', '20070219', '20070326', '20070430', '20070604',
							'20070709', '20070813', '20070917']
		exp_dates = [str2date(d) for d in raw_date]
		exp_repeat = [1, 1, 3, 3, 4, 3, 3, 3, 3, 3, 3, 2, 2]

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
