
from os.path import join 
import glob, unittest
from math import pi, cos, sin, radians

import numpy
from numpy import array, nan, reshape, squeeze
from numpy.testing import assert_array_almost_equal

import algorithm
from shared import Ifg


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


class MSTTests(unittest.TestCase):
	
	def setUp(self):
		self.testdir = "../../tests/sydney_test/obs"
		self.datafiles = glob.glob( join(self.testdir, "*.unw") )
		
	
	def test_temp_mst(self):
		ifgs = [Ifg(i) for i in self.datafiles]
		res = algorithm.temp_mst(ifgs)
		self.assertTrue(res is not None)
		
		raise NotImplementedError("Test results of MST")
		