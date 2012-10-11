'''
Created on 12/09/2012
@author: bpd900
'''

import unittest, os

from gdal import Dataset, UseExceptions
UseExceptions()

import shared


class SharedTests(unittest.TestCase):

	def test_create_ifg(self):
		ifg = shared.Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')
		self.assertTrue(os.path.exists(ifg.hdr_path)) # validate header path

	def test_headers_as_attr(self):
		attrs = ['WIDTH', 'FILE_LENGTH', 'X_FIRST', 'X_STEP',
			'Y_FIRST', 'Y_STEP', 'WAVELENGTH', 'MASTER', 'SLAVE']

		ifg = shared.Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')
		for a in attrs:
			self.assertTrue(getattr(ifg, a))

	def test_open(self):
		ifg = shared.Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')
		self.assertTrue(ifg.dataset is None)
		ifg.open()
		self.assertTrue(ifg.dataset is not None)
		self.assertTrue(isinstance(ifg.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(shared.IfgException, ifg.open)
		os.remove(ifg.ehdr_path)


if __name__ == "__main__":
	unittest.main()
