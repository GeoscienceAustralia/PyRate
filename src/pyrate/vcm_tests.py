'''
Tests for the PyRate's Variance/Covariance matrix functionality.

Created on 14/03/2013
@author: Ben Davies
'''

import unittest
from numpy import array

from shared import Ifg
from vcm import cvd
from tests_common import sydney5_mock_ifgs


class VCMTests(unittest.TestCase):

	def test_basic(self):
		ifgs = sydney5_mock_ifgs(5,9)

		for i in ifgs:
			i.X_SIZE = 92.0
			i.Y_SIZE = 95.0
			i.X_CENTRE = (i.WIDTH / 2) # spectrally, the grid is not doubled in size
			i.Y_CENTRE = (i.FILE_LENGTH / 2)

		maxvar, alpha = cvd(ifgs[0])
		self.assertTrue(maxvar is not None)
		self.assertTrue(alpha is not None)
		print "maxvar, alpha:", maxvar, alpha


	def test_larger(self):
		# test with larger interferogram
		ifg = Ifg("../../tests/sydney_test/obs/geo_060619-061002.unw")
		ifg.open()

		if bool((ifg.phase_data == 0).all()) is True:
			raise Exception("All zero - aaaieeee")

		maxvar, alpha = cvd(ifg)
		self.assertTrue(maxvar is not None)
		self.assertTrue(alpha is not None)
		print "maxvar, alpha:", maxvar, alpha
