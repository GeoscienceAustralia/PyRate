'''
Contains objects common to multiple parts of PyRate

Created on 12/09/2012
@author: bpd900
'''

import os
import numpy

try:
	from osgeo import gdal
	from gdalconst import GA_Update
except:
	import gdal, gdalconst

gdal.UseExceptions()

import roipac
from ifgconstants import DATUM


# Constants
AMPLITUDE_BAND = 1
PHASE_BAND = 2

# TODO: add phase_data and amplitude_data properties?
#     Problem: property links to FULL dataset, which may be slower than row by row access
#         row by row access would be efficient, but needes wavelength converted layers
#         NB: can just by row regardless, through the numpy interface: for row in data: ...


class RasterBase(object):
	'''Base class for ROIPAC format raster datasets.'''

	def __init__(self, path, hdr_path=None):
		'''Handles common task of bundling various forms of ROIPAC header files with
		a binary data layer.'''

		if hdr_path:
			# handle non default header (eg. for look files with different naming)
			self.data_path, self.hdr_path = path, hdr_path
		else:
			# default the header path
			self.data_path, self.hdr_path = roipac.filename_pair(path)

		# dynamically include header items as class attrs
		header = roipac.parse_header(self.hdr_path)
		self.__dict__.update(header)

		self.ehdr_path = None # path to EHdr format header
		self.dataset = None # for GDAL dataset obj
		self._readonly = None
		self.num_cells = None


	def __str__(self):
		name = self.__class__.__name__
		return "%s('%s')" % (name, self.data_path)


	def __repr__(self):
		name = self.__class__.__name__
		return "%s('%s', '%s')" % (name, self.data_path, self.hdr_path)


	def open(self, readonly=True):
		'''Opens generic raster dataset. Creates ESRI/EHdr format header in the data
		dir, creating a recogniseable header file for GDAL (as per ROIPAC doco).'''
		if self.ehdr_path is None:
			self.ehdr_path = roipac.to_ehdr_header(self.hdr_path)
			args = (self.data_path,) if readonly else (self.data_path, GA_Update)
			self.dataset = gdal.Open(*args)

			if self.dataset is None:
				raise RasterException("Error opening %s" % self.data_path)

		else:
			if self.dataset is not None:
				msg = "open() already called for %s" % self
				raise RasterException(msg)

		self._readonly = readonly
		self.num_cells = self.dataset.RasterYSize * self.dataset.RasterXSize


	@property
	def is_open(self):
		return self.dataset is not None


	def _get_band(self, band):
		'''Wrapper (with error checking) for GDAL's Band.GetRasterBand() method.'''
		if self.dataset is not None:
			return self.dataset.GetRasterBand(band)
		else:
			raise RasterException("Raster %s has not been opened" % self.data_path)



class Ifg(RasterBase):
	"""Interferogram class, representing the difference between two acquisitions.
	Ifg objects double as a container for related data."""

	def __init__(self, path, hdr_path=None):
		'''Interferogram constructor, for 2 band ROIPAC Ifg raster datasets.'''
		RasterBase.__init__(self, path, hdr_path)
		self._amp_band = None
		self._phase_band = None

		# creating code needs to set this flag after 0 -> NaN replacement
		self.nan_converted = False

		# TODO: what are these for?
		self.max_variance = None
		self.alpha = None
		self._nan_fraction = None


	@property
	def amp_band(self):
		'''Returns a GDAL Band object for the amplitude band'''
		if self._amp_band is None:
			self._amp_band = self._get_band(AMPLITUDE_BAND)
		return self._amp_band


	@property
	def phase_band(self):
		'''Returns a GDAL Band object for the phase band'''
		if self._phase_band is None:
			self._phase_band = self._get_band(PHASE_BAND)
		return self._phase_band


	@property
	def phase_data(self):
		'''Returns entire phase band as an array'''
		return self.phase_band.ReadAsArray()


	@property
	def phase_rows(self):
		'''TODO'''
		# TODO: is a pre-created buffer more efficient?
		for y in xrange(self.FILE_LENGTH):
			row = self.phase_band.ReadAsArray(yoff=y, win_xsize=self.WIDTH, win_ysize=1)
			yield row[0] # squeezes row from (1, WIDTH) to 1D array


	@property
	def nan_fraction(self):
		'''Returns 0-1 (float) proportion of NaN cells for the phase band'''

		# TODO: cache nan_count for readonly datasets? Perf benefit vs temp changes to data?
		data = self.phase_band.ReadAsArray()
		nan_count = numpy.sum(numpy.isnan(data))

		# handle datasets with no 0 -> NaN replacement
		if self.nan_converted is False and nan_count == 0:
			nan_count = numpy.sum(data == 0)

		return nan_count / float(self.num_cells)



class DEM(RasterBase):
	"""Generic raster class for ROIPAC single band DEM files"""

	def __init__(self, path, hdr_path=None):
		'''DEM constructor.'''
		RasterBase.__init__(self, path, hdr_path)
		self._band = None


	@property
	def height_band(self):
		if self._band is None:
			self._band = self._get_band(1)
		return self._band



class IfgException(Exception):
	'''Generic exception class for interferogram errors'''
	pass

class RasterException(Exception):
	'''Generic exception for raster errors'''
	pass

class PyRateException(Exception):
	'''Generic exception class for PyRate S/W errors'''
	pass


class EpochList(object):
	'''TODO'''

	def __init__(self, dates=None, repeat=None, spans=None):
		self.dates = dates # list of unique dates from all the ifgs
		self.repeat = repeat
		self.spans = spans # time span from earliest ifg
