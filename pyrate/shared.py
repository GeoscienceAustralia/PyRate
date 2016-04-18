"""
Contains objects common to multiple parts of PyRate

Created on 12/09/2012

.. codeauthor:: Ben Davies, Sudipta Basak
"""

import os
import math
from datetime import date
import logging
from numpy import where, nan, isnan, sum as nsum
import numpy as np
import random
import string
from functools import wraps
import time
import logging

import pyrate.ifgconstants as ifc

try:
    from osgeo import gdal
    from osgeo.gdalconst import GA_Update, GA_ReadOnly, GF_Write
except ImportError:
    import gdal

gdal.UseExceptions()

from pyrate.geodesy import cell_size

# Constants
PHASE_BAND = 1
META_UNITS = 'PHASE_UNITS'
MILLIMETRES = 'MILLIMETRES'
MM_PER_METRE = 1000

# GDAL projection list
GDAL_X_CELLSIZE = 1
GDAL_Y_CELLSIZE = 5
GDAL_X_FIRST = 0
GDAL_Y_FIRST = 3


class RasterBase(object):
    """
    Base class for PyRate GeoTIFF based raster datasets.
    """

    def __init__(self, path):
        self.data_path = path
        self.dataset = None # for GDAL dataset obj
        self._readonly = not os.access(path, os.R_OK | os.W_OK)

        if self._readonly is None:
            raise NotImplementedError # os.access() has failed?

    def __str__(self):
        name = self.__class__.__name__
        return "%s('%s')" % (name, self.data_path)

    def __repr__(self):
        name = self.__class__.__name__
        return "%s('%s')" % (name, self.data_path)

    def open(self, readonly=None):
        """
        Opens generic raster dataset.
        """
        if self.dataset is not None:
            msg = "open() already called for %s" % self
            raise RasterException(msg)

        # unless read only, by default open files as writeable
        if readonly not in [True, False, None]:
            raise ValueError("readonly must be True, False or None")

        if readonly is False and self._readonly is True:
            raise IOError("Cannot open write protected file for writing")

        flag = GA_ReadOnly if self._readonly else GA_Update

        if not os.path.exists(self.data_path):
            raise IOError('The file {path} does not exist.'
                    ' Consider running prepifg'.format(path=self.data_path))

        self.dataset = gdal.Open(self.data_path, flag)

        if self.dataset is None:
            raise RasterException("Error opening %s" % self.data_path)

        # add some geographic data
        self.x_centre = self.ncols / 2
        self.y_centre = self.nrows / 2
        self.lat_centre = self.y_first + (self.y_step * self.y_centre)
        self.long_centre = self.x_first + (self.x_step * self.x_centre)

        # use cell size from centre of scene
        self.x_size, self.y_size = cell_size(self.lat_centre, self.long_centre,
                                            self.x_step, self.y_step)

    @property
    def ncols(self):
        return self.dataset.RasterXSize

    @property
    def nrows(self):
        return self.dataset.RasterYSize

    @property
    def x_step(self):        
        return float(self.dataset.GetGeoTransform()[GDAL_X_CELLSIZE])

    @property
    def y_step(self):
        return float(self.dataset.GetGeoTransform()[GDAL_Y_CELLSIZE])

    @property
    def x_first(self):
        return float(self.dataset.GetGeoTransform()[GDAL_X_FIRST])

    @property
    def x_last(self):
        return self.x_first + (self.x_step * self.ncols)

    @property
    def y_first(self):
        return float(self.dataset.GetGeoTransform()[GDAL_Y_FIRST])

    @property
    def y_last(self):
        return self.y_first + (self.y_step * self.nrows)

    @property
    def shape(self):
        """
        Returns tuple of (Y,X) shape of the raster (as per numpy.shape).
        """

        return self.dataset.RasterYSize, self.dataset.RasterXSize

    @property
    def num_cells(self):
        if self.is_open:
            return self.dataset.RasterXSize * self.dataset.RasterYSize
        else:
            raise RasterException('Dataset not open')

    @property
    def is_open(self):
        """
        Returns True if the underlying dataset has been opened by GDAL.
        """

        return self.dataset is not None

    def close(self):
        """
        Explicitly closes file opened by gdal.Open()
        This is required in windows, otherwise opened files can not be removed,
        because windows locks opened files.
        """
        if self.is_open:
            self.dataset = None

    @property
    def is_read_only(self):
        return self._readonly

    def _get_band(self, band):
        """
        Wrapper (with error checking) for GDAL's Band.GetRasterBand() method.

        :param band: number of band, starting at 1
        """

        if self.dataset is not None:
            return self.dataset.GetRasterBand(band)
        else:
            raise RasterException("Raster %s has not been opened"
                                  % self.data_path)


class Ifg(RasterBase):
    """
    Interferogram class, represents the difference between two acquisitions.
    Ifg objects double as a container for related data.
    """

    def __init__(self, path):
        """
        Interferogram constructor, for 2 band ROIPAC Ifg raster datasets.
        """

        RasterBase.__init__(self, path)
        self._phase_band = None
        self._phase_data = None
        self.master = None
        self.slave = None
        self.nan_converted = False
        self.mm_converted = False
        self.meta_data = None
        self.wavelength = None

    def open(self, readonly=None):
        """
        :param bool readonly: True/False, or None to open as underlying file setting
        """

        RasterBase.open(self, readonly)
        self._init_dates()

        md = self.dataset.GetMetadata()
        self.wavelength = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        self.meta_data = md

        # creating code needs to set this flag after 0 -> NaN replacement
        self.nan_converted = False

    def _init_dates(self):
        def _to_date(datestr):
            year, month, day = [int(i) for i in datestr.split('-')]
            return date(year, month, day)

        md = self.dataset.GetMetadata()
        datestrs = [md[k] for k in [ifc.PYRATE_DATE, ifc.PYRATE_DATE2]]

        if all(datestrs):
            self.master, self.slave = [_to_date(s) for s in datestrs]
            self.time_span = (self.slave - self.master).days / ifc.DAYS_PER_YEAR
        else:
            msg = 'Missing master and/or slave date in %s' % self.data_path
            raise IfgException(msg)

    def convert_to_nans(self, val=0.0):
        """
        Converts given values in phase data to NaNs
        :param val: value to convert, default is 0
        """
        self.phase_data = where(np.isclose(self.phase_data, val, atol=1e-6),
                                nan, self.phase_data)
        self.nan_converted = True

    @property
    def phase_band(self):
        """
        Returns a GDAL Band object for the phase band.
        """

        if self._phase_band is None:
            self._phase_band = self._get_band(PHASE_BAND)
        return self._phase_band

    @property
    def phase_data(self):
        """
        Returns entire phase band as an array.
        """

        if self._phase_data is None:
            self._phase_data = self.phase_band.ReadAsArray()
        return self._phase_data

    def convert_to_mm(self):
        """
        :param ifg: ifg file
        :return: convert wavelength from radians to mm
        """
        self.mm_converted = True
        if self.dataset.GetMetadataItem(META_UNITS) == MILLIMETRES:
            msg = '%s: ignored as previous wavelength conversion detected'
            logging.debug(msg % self.data_path)
            self.phase_data = self.phase_data
            return

        self.phase_data = wavelength_radians_to_mm(self.phase_data,
                                                      self.wavelength)
        self.dataset.SetMetadataItem(META_UNITS, MILLIMETRES)
        msg = '%s: converted wavelength to millimetres'
        logging.debug(msg % self.data_path)

    @phase_data.setter
    def phase_data(self, data):
        self._phase_data = data

    @property
    def phase_rows(self):
        """
        Generator returning each row of the phase data.
        """

        for y in xrange(self.nrows):
            r = self.phase_band.ReadAsArray(yoff=y,
                                            win_xsize=self.ncols, win_ysize=1)
            yield r[0] # squeezes row from (1, WIDTH) to 1D array

    @property
    def nan_count(self):
        """
        Returns number of NaN cells in the phase data.
        """
        return nsum(isnan(self.phase_data))

    @property
    def nan_fraction(self):
        """
        Returns 0-1 (float) proportion of NaN cells for the phase band.
        """
        # don't cache nan_count as client code may modify phase data
        nan_count = self.nan_count

        # handle datasets with no 0 -> NaN replacement
        if self.nan_converted is False and nan_count == 0:
            nan_count = nsum(np.isclose(self.phase_data, 0.0, atol=1e-6))
        return nan_count / float(self.num_cells)

    def write_modified_phase(self):
        """
        Writes phase data to disk.
        For this to work, a copy of the original file
        """

        if self.is_read_only:
            raise IOError("Cannot write to read only Ifg")
        """
        # keep this block
        if new_data_path is None:
            self.dataset = gdal.Open(self.data_path, GA_Update)
        else:
            self.dataset = gdal.Open(new_data_path, GA_Update)
        self._phase_band = None
        """
        self.phase_band.WriteArray(self.phase_data)


class Incidence(RasterBase):

    def __init__(self, path):
        """
        Incidence obj constructor.
        """

        RasterBase.__init__(self, path)
        self._incidence_band = None
        self._azimuth_band = None
        self._incidence_data = None
        self._azimuth_data = None


    @property
    def incidence_band(self):
        '''
        Returns the GDALBand for the incidence angle layer.
        '''

        if self._incidence_band is None:
            self._incidence_band = self._get_band(1)
        return self._incidence_band


    @property
    def incidence_data(self):
        """
        Returns the entire incidence band as an array.
        """

        if self._incidence_data is None:
            self._incidence_data = self.incidence_band.ReadAsArray()
        return self._incidence_data


    @property
    def azimuth_band(self):
        """
        Returns the GDALBand for the azimuth layer.
        """

        if self._azimuth_band is None:
            self._azimuth_band = self._get_band(2)
        return self._azimuth_band


    @property
    def azimuth_data(self):
        """
        Returns the entire incidence band as an array.
        """

        if self._azimuth_data is None:
            self._azimuth_data = self.azimuth_band.ReadAsArray()
        return self._azimuth_data



class DEM(RasterBase):
    """
    Generic raster class for ROIPAC single band DEM files.
    """

    def __init__(self, path):
        '''DEM constructor.'''
        RasterBase.__init__(self, path)
        self._band = None


    @property
    def height_band(self):
        """
        Returns the GDALBand for the elevation layer.
        """

        if self._band is None:
            self._band = self._get_band(1)
        return self._band



class IfgException(Exception):
    """
    Generic exception class for interferogram errors.
    """

    pass

class RasterException(Exception):
    """
    Generic exception for raster errors.
    """

    pass

class PyRateException(Exception):
    """
    Generic exception class for PyRate S/W errors.
    """

    pass


class EpochList(object):
    """
    Metadata container for epoch related information.
    """

    def __init__(self, dates=None, repeat=None, spans=None):
        self.dates = dates # list of unique dates from all the ifgs
        self.repeat = repeat
        self.spans = spans # time span from earliest ifg

    def __str__(self):
        return "EpochList: %s" % str(self.dates)

    def __repr__(self):
        return "EpochList: %s" % repr(self.dates)


def wavelength_radians_to_mm(data, wavelength):
    """
    Translates phase from radians to millimetres
    data: ifg phase data
    wavelength: normally included with SAR instrument pass data
    """

    # '4' is 2*2, the 1st 2 is that the raw signal is 'there and back', to get
    # the vector length between satellite and ground, half the signal is needed
    # second 2*pi is because one wavelength is equal to 2 radians
    return data * MM_PER_METRE * (wavelength / (4 * math.pi))


def generate_random_string(N=10):
    return ''.join(random.SystemRandom().choice(
        string.ascii_letters + string.digits)
                   for _ in range(N))


def timer(f):
    """
    A simple timing decorator for the entire process.

    """
    @wraps(f)
    def wrap(*args, **kwargs):
        t1 = time.time()
        res = f(*args, **kwargs)

        tottime = time.time() - t1
        msg = "%02d:%02d:%02d " % \
              reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
                     [(tottime,), 60, 60])

        logging.info("Time for {0}: {1}".format(f.func_name, msg))
        return res

    return wrap