'''
Unittests for PyRate shared utility objects.
Created on 12/09/2012
@author: Ben Davies, NCI
'''

import os
import shutil
import unittest
from os.path import join, basename
from stat import S_IRGRP, S_IWGRP, S_IWOTH, S_IROTH, S_IRUSR, S_IWUSR

from itertools import product
from numpy import isnan, where, nan
from numpy.testing import assert_array_equal

from gdal import Open, Dataset, UseExceptions
UseExceptions()

from pyrate.shared import Ifg, DEM, RasterException
from common import SYD_TEST_TIF, SYD_TEST_DEM_TIF


class IfgTests(unittest.TestCase):
	'''Unit tests for the Ifg/interferogram class.'''

	def setUp(self):
		self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_060619-061002.tif'))
		self.ifg.open()

	def test_headers_as_attr(self):
		for a in ['ncols', 'nrows', 'x_first', 'x_step',
				'y_first', 'y_step', 'wavelength', 'master', 'slave']:
			self.assertTrue(getattr(self.ifg, a) is not None)

	def test_convert_to_nans(self):
		self.ifg.convert_to_nans(0)
		self.assertTrue(self.ifg.nan_converted)

	def test_xylast(self):
		# ensure the X|Y_LAST header element has been created
		self.assertAlmostEqual(self.ifg.x_last, 150.9491667)
		self.assertAlmostEqual(self.ifg.y_last, -34.23)

	def test_num_cells(self):
		# test cell size from header elements
		data = self.ifg.phase_band.ReadAsArray()
		ys, xs = data.shape
		exp_ncells = ys * xs
		self.assertEqual(exp_ncells, self.ifg.num_cells)

	def test_shape(self):
		self.assertEqual(self.ifg.shape, self.ifg.phase_data.shape)

	def test_nan_count(self):
		num_nan = 0
		for row in self.ifg.phase_data:
			for v in row:
				if isnan(v):
					num_nan += 1

		self.assertEqual(num_nan, self.ifg.nan_count)

	def test_phase_band(self):
		data = self.ifg.phase_band.ReadAsArray()
		self.assertEqual(data.shape, (72, 47) )

	def test_nan_fraction(self):
		# NB: source data lacks 0 -> NaN conversion
		data = self.ifg.phase_data
		data = where(data == 0, nan, data) # fake 0 -> nan for the count below

		# manually count # nan cells
		nans = 0
		ys, xs = data.shape
		for y, x in product(xrange(ys), xrange(xs)):
			if isnan(data[y, x]):
				nans += 1
		del data

		num_cells = float(ys * xs)
		self.assertTrue(nans > 0)
		self.assertTrue(nans <= num_cells)
		self.assertEqual(nans / num_cells, self.ifg.nan_fraction)

	def test_xy_size(self):
		self.assertFalse(self.ifg.ncols is None)
		self.assertFalse(self.ifg.nrows is None)

		# test with tolerance from base 90m cell
		self.assertTrue(self.ifg.y_size > 88.0) # within 2% of cells over sydney?
		self.assertTrue(self.ifg.y_size < 92.0, 'Got %s' % self.ifg.y_size)

		syd_width = 76.9 # from nearby pirate coords
		self.assertTrue(self.ifg.x_size > 0.97 * syd_width) # ~3% tolerance
		self.assertTrue(self.ifg.x_size < 1.03 * syd_width)

	def test_centre_latlong(self):
		lat_exp = self.ifg.y_first + ((self.ifg.nrows / 2) * self.ifg.y_step)
		long_exp = self.ifg.x_first + ((self.ifg.ncols / 2) * self.ifg.x_step)
		self.assertEqual(lat_exp, self.ifg.lat_centre)
		self.assertEqual(long_exp, self.ifg.long_centre)

	def test_centre_cell(self):
		self.assertEqual(self.ifg.x_centre, 23)
		self.assertEqual(self.ifg.y_centre, 36)

	def test_time_span(self):
		self.assertAlmostEqual(self.ifg.time_span, 0.287474332649)

	def test_wavelength(self):
		self.assertEqual(self.ifg.wavelength, 0.0562356424)


class IfgIOTests(unittest.TestCase):

	def setUp(self):
		self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_060619-061002.tif'))

	def test_open(self):
		self.assertTrue(self.ifg.dataset is None)
		self.assertTrue(self.ifg.is_open is False)
		self.ifg.open(readonly=True)
		self.assertTrue(self.ifg.dataset is not None)
		self.assertTrue(self.ifg.is_open is True)
		self.assertTrue(isinstance(self.ifg.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(RasterException, self.ifg.open, True)

	def test_write(self):
		base = "/tmp"
		src = self.ifg.data_path
		dest = join(base, basename(self.ifg.data_path))

		# shutil.copy needs to copy writeable permission from src
		os.chmod(src, S_IRGRP | S_IWGRP | S_IWOTH | S_IROTH | S_IRUSR | S_IWUSR)
		shutil.copy(src, dest)
		os.chmod(src, S_IRGRP | S_IROTH | S_IRUSR) # revert

		i = Ifg(dest)
		i.open()
		i.phase_data[0,1:] = nan
		i.write_phase()
		del i

		# reopen to ensure data/nans can be read back out
		i = Ifg(dest)
		i.open(readonly=True)
		assert_array_equal(True, isnan(i.phase_data[0,1:]) )
		os.remove(dest)


	def test_readonly_permission_failure(self):
		# ensure failure if opening R/O permission file as writeable/GA_Update
		self.assertRaises(IOError, self.ifg.open, False)


	def test_write_fails_on_readonly(self):
		# check readonly status is same before and after open() for readonly file
		self.assertTrue(self.ifg.is_read_only)
		self.ifg.open(readonly=True)
		self.assertTrue(self.ifg.is_read_only)
		self.assertRaises(IOError, self.ifg.write_phase)

	def test_phase_band_unopened_ifg(self):
		try:
			_ = self.ifg.phase_band
			self.fail("Should not be able to access band without open dataset")
		except RasterException:
			pass

	def test_nan_fraction_unopened(self):
		try:
			# NB: self.assertRaises doesn't work here (as it is a property?)
			_ = self.ifg.nan_fraction
			self.fail("Shouldn't be able to call nan_fraction() with unopened Ifg")
		except RasterException:
			pass

	def test_phase_data_properties(self):
		# Use raw GDAL to isolate raster reading from Ifg functionality
		ds = Open(self.ifg.data_path)
		data = ds.GetRasterBand(1).ReadAsArray()
		del ds

		self.ifg.open()

		# test full array and row by row access
		assert_array_equal(data, self.ifg.phase_data)
		for y, row in enumerate(self.ifg.phase_rows):
			assert_array_equal(data[y], row)

		# test the data is cached if changed
		crd = (5,4)
		orig = self.ifg.phase_data[crd]
		self.ifg.phase_data[crd] *= 2
		nv = self.ifg.phase_data[crd] # pull new value out again
		self.assertEqual(nv, 2 * orig)


# FIXME: 
# class IncidenceFileTests(unittest.TestCase):
# 	'Unit tests to verify operations on GeoTIFF format Incidence rasters'
# 	
# 	def setUp(self):
# 		raise NotImplementedError
# 		self.inc = Incidence(join(INCID_TEST_DIR, '128x2.tif'))
# 		self.inc.open()
# 
# 
# 	def test_incidence_data(self):
# 		# check incidences rises while traversing the scene
# 		data = self.inc.incidence_data
# 		diff = data.ptp()
# 		self.assertTrue(diff > 0.5, "Got ptp() diff of %s" % diff)
# 
# 		# ascending pass, values should increase from W->E across scene
# 		for i in range(2):
# 			d = data[i]
# 			self.assertFalse((d == 0).any()) # ensure no NODATA
# 			self.assertFalse((isnan(d)).any())
# 
# 			diff = array([d[i+1] - d[i] for i in range(len(d)-1)])
# 			res = abs(diff[diff < 0])
# 			self.assertTrue((res < 1e-4).all()) # TODO: check if this is normal
# 
# 
# 	def test_azimuth_data(self):
# 		# ensure azimuth is fairly constant
# 
# 		az = self.inc.azimuth_data
# 		self.assertFalse((az == 0).all())
# 		az = az[az != 0] # filter NODATA cells
# 
# 		# azimuth should be relatively constant
# 		ptp = az.ptp()
# 		self.assertTrue(ptp < 0.1, msg="min -> max diff is %s" % ptp)


class DEMTests(unittest.TestCase):
	'Unit tests to verify operations on GeoTIFF format DEMs'

	def setUp(self):
		self.ras = DEM(SYD_TEST_DEM_TIF)


	def test_create_raster(self):
		self.assertTrue(os.path.exists(self.ras.data_path)) # validate header path


	def test_headers_as_attr(self):
		self.ras.open()
		attrs = ['ncols', 'nrows', 'x_first', 'x_step', 'y_first', 'y_step' ]

		# TODO: are 'projection' and 'datum' attrs needed?
		for a in attrs:
			self.assertTrue(getattr(self.ras, a) is not None)


	def test_is_dem(self):
		self.ras = DEM(join(SYD_TEST_TIF, 'geo_060619-061002.tif'))
		self.assertFalse(hasattr(self.ras, 'datum'))


	def test_open(self):
		self.assertTrue(self.ras.dataset is None)
		self.ras.open()
		self.assertTrue(self.ras.dataset is not None)
		self.assertTrue(isinstance(self.ras.dataset, Dataset))

		# ensure open cannot be called twice
		self.failUnlessRaises(RasterException, self.ras.open)


	def test_band(self):
		# test accessing bands with open and unopened datasets
		try:
			_ = self.ras.height_band
			self.fail("Should not be able to access band without open dataset")
		except RasterException:
			pass

		self.ras.open()
		data = self.ras.height_band.ReadAsArray()
		self.assertEqual(data.shape, (72, 47) )


if __name__ == "__main__":
	unittest.main()
