'''
Created on 12/09/2012
@author: bpd900
'''

import unittest

import shared


class SharedTests(unittest.TestCase):

	def test_create_ifg(self):
		ifg = shared.Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')

	def test_headers_as_attr(self):
		attrs = ['WIDTH', 'FILE_LENGTH', 'X_FIRST', 'X_STEP',
			'Y_FIRST', 'Y_STEP', 'WAVELENGTH', 'MASTER', 'SLAVE']

		ifg = shared.Ifg('../../tests/sydney_test/obs/geo_060619-061002.unw')
		for a in attrs:
			self.assertTrue(getattr(ifg, a))


if __name__ == "__main__":
	unittest.main()
