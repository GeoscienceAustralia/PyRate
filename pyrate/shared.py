"""
Contains objects common to multiple parts of PyRate

Created on 12/09/2012

.. codeauthor:: Ben Davies, Sudipta Basak, Matt Garthwaite
"""
import errno

import os, struct
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
import pkg_resources
from pyrate import roipac, gamma

import pyrate.ifgconstants as ifc

try:
    from osgeo import osr, gdal
    from osgeo.gdalconst import GA_Update, GA_ReadOnly, GF_Write
except ImportError:
    import gdal

gdal.UseExceptions()

from pyrate.geodesy import cell_size

# Constants
PHASE_BAND = 1
RADIANS = 'RADIANS'
MILLIMETRES = 'MILLIMETRES'
MM_PER_METRE = 1000
GAMMA = 'GAMMA'
ROIPAC = 'ROIPAC'

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
        if isinstance(path, gdal.Dataset):
            self.dataset = path  # path will be Dataset in this case
            self.data_path = self.dataset  # data_path dummy
            self.add_geographic_data()
        else:
            self.data_path = path
            self.dataset = None  # for GDAL dataset obj
            self._readonly = not os.access(path, os.R_OK | os.W_OK)

            if self._readonly is None:
                raise NotImplementedError  # os.access() has failed?

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

        if not os.path.exists(self.data_path):
            raise IOError('The file {path} does not exist.'
                    ' Consider running prepifg'.format(path=self.data_path))

        # unless read only, by default open files as writeable
        if readonly not in [True, False, None]:
            raise ValueError("readonly must be True, False or None")

        if readonly is False and self._readonly is True:
            raise IOError("Cannot open write protected file for writing")

        flag = GA_ReadOnly if self._readonly else GA_Update

        self.dataset = gdal.Open(self.data_path, flag)

        if self.dataset is None:
            raise RasterException("Error opening %s" % self.data_path)

        self.add_geographic_data()

    def add_geographic_data(self):
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
        self._nodata_value = None

    def open(self, readonly=None):
        """
        :param bool readonly: True/False, or None to open as underlying file setting
        """
        RasterBase.open(self, readonly)
        self.initialize()

    def initialize(self):
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

    def convert_to_nans(self):
        """
        Converts given values in phase data to NaNs
        :param val: value to convert, default is 0
        """
        if (self._nodata_value is None) or (self.dataset is None):
            msg = 'nodata value needs to be set for nan conversion.' \
                  'Use ifg.nodata_value = NoDataValue to set nodata_value'
            raise RasterException(msg)

        self.phase_data = where(np.isclose(self.phase_data,
                                           self._nodata_value, atol=1e-6),
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
    def nodata_value(self):
        """
        Returns a GDAL Band object for the phase band.
        """
        return self._nodata_value

    @nodata_value.setter
    def nodata_value(self, val):
        self._nodata_value = val

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
        if self.dataset.GetMetadataItem(ifc.PYRATE_PHASE_UNITS) == MILLIMETRES:
            msg = '%s: ignored as previous phase unit conversion already applied'
            logging.debug(msg % self.data_path)
            self.phase_data = self.phase_data
            return
        #elif self.dataset.GetMetadataItem(ifc.PYRATE_PHASE_UNITS) == RADIANS:
        self.phase_data = convert_radians_to_mm(self.phase_data,
                                                      self.wavelength)
        self.dataset.SetMetadataItem(ifc.PYRATE_PHASE_UNITS, MILLIMETRES)
        msg = '%s: converted phase units to millimetres'
        logging.debug(msg % self.data_path)
        # TODO: implement test for when units neither mm or radians
        #else:
        #    msg = 'Phase units are not millimetres or radians'
        #    raise IfgException(msg)


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
        if (self._nodata_value is None) or (self.dataset is None):
            msg = 'nodata_value needs to be set for nan fraction calc.' \
                  'Use ifg.nondata = NoDataValue to set nodata'
            raise RasterException(msg)
        # don't cache nan_count as client code may modify phase data
        nan_count = self.nan_count
        # handle datasets with no 0 -> NaN replacement
        if not self.nan_converted and (nan_count == 0):
            nan_count = nsum(np.isclose(self.phase_data,
                                        self._nodata_value, atol=1e-6))
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


class IfgPart(object):
    """
    slice of Ifg data object
    """

    def __init__(self, data_path, r_start, r_end, c_start, c_end):

        """
        :param ifg: original ifg
        :param r_start: starting tow of the original ifg
        :param r_end: ending row of the original ifg
        :return:
        """
        self.data_path = data_path
        ifg = Ifg(data_path)
        ifg.open(readonly=True)
        ifg.nodata_value = 0
        self._phase_data = ifg.phase_data
        self.nan_fraction = ifg.nan_fraction
        self.master = ifg.master
        self.slave = ifg.slave
        ifg.close()
        self.r_start = r_start
        self.r_end = r_end
        self._phase_data_part = None
        self.c_start = c_start
        self.c_end = c_end

    @property
    def nrows(self):
        return self.r_end - self.r_start

    @property
    def ncols(self):
        return self.c_end - self.c_start

    @property
    def phase_data(self):
        if self._phase_data_part is None:
            return self._phase_data[self.r_start:self.r_end,
                   self.c_start:self.c_end]


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


def convert_radians_to_mm(data, wavelength):
    """
    Translates phase from radians to millimetres
    data: interferogram phase data
    wavelength: radar wavelength; normally included with SAR instrument metadata
    """
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


def nanmedian(x):
    """
    :param x:
    :return:
    """
    version = [int(i) for i in
               pkg_resources.get_distribution("numpy").version.split('.')]
    if version[0] == 1 and version[1] > 9:
        return np.nanmedian(x)
    else:
        return np.median(x[~np.isnan(x)])


def write_geotiff(header, data_path, dest, nodata):
    """
    Writes image data to GeoTIFF image with PyRate metadata
    """
    is_ifg = ifc.PYRATE_WAVELENGTH_METRES in header
    ifg_proc = header[ifc.PYRATE_INSAR_PROCESSOR]
    ncols = header[ifc.PYRATE_NCOLS]
    nrows = header[ifc.PYRATE_NROWS]

    # need to have gamma/roipac functionality here?
    if ifg_proc == ROIPAC:
        roipac._check_raw_data(is_ifg, data_path, ncols, nrows)
        roipac._check_step_mismatch(header)
    else: #GAMMA
        gamma._check_raw_data(data_path, ncols, nrows)
        gamma._check_step_mismatch(header)

    driver = gdal.GetDriverByName("GTiff")
    dtype = gdal.GDT_Float32 if is_ifg else gdal.GDT_Int16
    ds = driver.Create(dest, ncols, nrows, 1, dtype)

    # write pyrate parameters to headers
    if is_ifg:
        for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                  ifc.PYRATE_INSAR_PROCESSOR,
                  ifc.PYRATE_DATE, ifc.PYRATE_DATE2, ifc.PYRATE_PHASE_UNITS]:
            ds.SetMetadataItem(k, str(header[k]))

    # position and projection data
    ds.SetGeoTransform([header[ifc.PYRATE_LONG], header[ifc.PYRATE_X_STEP], 0,
                        header[ifc.PYRATE_LAT], 0, header[ifc.PYRATE_Y_STEP]])

    srs = osr.SpatialReference()
    res = srs.SetWellKnownGeogCS(header[ifc.PYRATE_DATUM])

    if res:
        msg = 'Unrecognised projection: %s' % header[ifc.PYRATE_DATUM]
        raise GeotiffException(msg)

    ds.SetProjection(srs.ExportToWkt())

    # copy data from the binary file
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    if ifg_proc == GAMMA:
        fmtstr = '!' + ('f' * ncols)  # data format is big endian float32s
        bytes_per_col = 4
    elif ifg_proc == ROIPAC:
        if is_ifg:
            fmtstr = '<' + ('f' * ncols)  # roipac ifgs are little endian float32s
            bytes_per_col = 4
        else:
            fmtstr = '<' + ('h' * ncols)  # roipac DEM is little endian signed int16
            bytes_per_col = 2
    else:
        msg = 'Unrecognised InSAR Processor: %s' % ifg_proc
        raise GeotiffException(msg)

    row_bytes = ncols * bytes_per_col

    with open(data_path, 'rb') as f:
        for y in range(nrows):
            if ifg_proc == ROIPAC:
                if is_ifg:
                    f.seek(row_bytes, 1)  # skip interleaved band 1

            data = struct.unpack(fmtstr, f.read(row_bytes))
            #else: # GAMMA
            #    data = struct.unpack(fmtstr, f.read(ncols * 4))

            band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)

    # Needed? Only in ROIPAC code
    ds = None  # manual close
    del ds


def write_output_geotiff(md, gt, wkt, data, dest, nodata):
    """
    Writes PyRate output data to a GeoTIFF file.
    md is a dictionary containing PyRate metadata
    gt is the GDAL geotransform for the data
    wkt is the GDAL projection information for the data
    """

    driver = gdal.GetDriverByName("GTiff")
    nrows, ncols = data.shape
    ds = driver.Create(dest, ncols, nrows, 1, gdal.GDT_Float32)
    # set spatial reference for geotiff
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    ds.SetMetadataItem(ifc.PYRATE_DATE, str(md[ifc.PYRATE_DATE]))

    # write data to geotiff
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)


class GeotiffException(Exception):
    pass

if __name__ == '__main__':
    print nanmedian([1, 2, nan, 2, nan, 4])


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise