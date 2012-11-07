
import unittest
from math import pi, cos, sin, radians

import numpy
from numpy import array, nan

import algorithm


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
		azimuth = radians(77.8)
		incidence = radians(34.3)
		vert = cos(incidence)
		ns = sin(incidence) * sin(azimuth)
		ew = sin(incidence) * cos(azimuth)
		unitv = [ew, ns, vert]
		
		# TODO: assumes rad input for now
		act = algorithm.unit_vector(array([incidence]), array([azimuth]))
		for a,e in zip(list(act), unitv):
			self.assertEqual(a,e)
