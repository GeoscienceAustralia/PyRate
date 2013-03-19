'''
Tests for the PyRate's Variance/Covariance matrix functionality.

Created on 14/03/2013
@author: Ben Davies
'''


import unittest

from numpy import array

from vcm import cvd
from tests_common import sydney5_mock_ifgs


class VCMTests(unittest.TestCase):

	def test_basic(self):
		ifgs = sydney5_mock_ifgs(5,9)

		#print
		#for i in ifgs:
		#	print i.phase_data
		#	print "/" * 35


		ifgs[0].X_SIZE = 92.0
		ifgs[0].Y_SIZE = 95.0

		ifgs[0].X_CENTRE = (ifgs[0].WIDTH / 2) # spectrally, the grid is not doubled in size
		ifgs[0].Y_CENTRE = (ifgs[0].FILE_LENGTH / 2)

		print cvd(ifgs[0])
