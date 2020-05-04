#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# coding: utf-8
"""
This Python module contains utilities and classes shared by
all other PyRate modules
"""
# pylint: disable=too-many-lines
import re
from typing import List, Union

import errno
import math
from math import floor
import os
from os.path import basename, join
import struct
from datetime import date
from itertools import product
import numpy as np
from numpy import where, nan, isnan, sum as nsum, isclose
import pyproj
import pkg_resources
try:
    from osgeo import osr, gdal
    from osgeo.gdalconst import GA_Update, GA_ReadOnly
except ImportError:
    import gdal

from pyrate.core import ifgconstants as ifc, mpiops, config as cf
from pyrate.core.logger import pyratelogger as log


gdal.UseExceptions()

VERBOSE = True

# Constants
PHASE_BAND = 1
RADIANS = 'RADIANS'
MILLIMETRES = 'MILLIMETRES'
GAMMA = 'GAMMA'
ROIPAC = 'ROIPAC'

# GDAL projection list
GDAL_X_CELLSIZE = 1
GDAL_Y_CELLSIZE = 5
GDAL_X_FIRST = 0
GDAL_Y_FIRST = 3


def joblib_log_level(level: str) -> int:
    """
    Convert python log level to joblib int verbosity.
    """ 
    if level == 'INFO':
        return 0
    else:
        return 60

def mkdir_p(path):
    """
    Make new directory and create parent directories as necessary.
    Copied from stackoverflow:
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python

    :param str path: Path name for new directory
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


class RasterBase(object):
    """
    Base class for PyRate GeoTIFF based raster datasets.
    """
    # pylint: disable=missing-docstring
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path: Union[gdal.Dataset, str]):
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
            raise IOError('The file {path} does not exist. Consider first running prepifg'.format(path=self.data_path))

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
        """
        Determine and add geographic data to object
        """
        # add some geographic data
        self.x_centre = int(self.ncols / 2)
        self.y_centre = int(self.nrows / 2)
        self.lat_centre = self.y_first + (self.y_step * self.y_centre)
        self.long_centre = self.x_first + (self.x_step * self.x_centre)
        # use cell size from centre of scene
        self.x_size, self.y_size = cell_size(self.lat_centre, self.long_centre, self.x_step, self.y_step)

    @property
    def ncols(self):
        """
        Number of raster columns
        """
        return self.dataset.RasterXSize

    @property
    def nrows(self):
        """
        Number of raster rows
        """
        return self.dataset.RasterYSize

    @property
    def x_step(self):
        """
        Raster pixel size in X (easting) dimension
        """
        return float(self.dataset.GetGeoTransform()[GDAL_X_CELLSIZE])

    @property
    def y_step(self):
        """
        Raster pixel size in Y (northing) dimension
        """
        return float(self.dataset.GetGeoTransform()[GDAL_Y_CELLSIZE])

    @property
    def x_first(self):
        """
        Raster western bounding coordinate
        """
        return float(self.dataset.GetGeoTransform()[GDAL_X_FIRST])

    @property
    def x_last(self):
        """
        Raster eastern bounding coordinate
        """
        return self.x_first + (self.x_step * self.ncols)

    @property
    def y_first(self):
        """
        Raster northern bounding coordinate
        """
        return float(self.dataset.GetGeoTransform()[GDAL_Y_FIRST])

    @property
    def y_last(self):
        """
        Raster southern bounding coordinate
        """
        return self.y_first + (self.y_step * self.nrows)

    @property
    def shape(self):
        """
        Returns tuple of (Y,X) shape of the raster (as per numpy.shape).
        """
        return self.dataset.RasterYSize, self.dataset.RasterXSize

    @property
    def num_cells(self):
        """
        Total number of pixels in raster dataset
        """
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
        because windows locks open files.
        """
        if self.is_open:
            self.dataset = None

    @property
    def is_read_only(self):
        """
        Determines file permissions
        """
        return self._readonly

    def _get_band(self, band):
        """
        Wrapper (with error checking) for GDAL's Band.GetRasterBand() method.

        :param int band: number of band, starting at 1
        """
        if self.dataset is not None:
            return self.dataset.GetRasterBand(band)
        else:
            raise RasterException("Raster %s has not been opened" % self.data_path)


class Ifg(RasterBase):
    """
    Interferogram (Ifg) class objects; double as a container for
    interferometric phase raster band data and related data.
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path):
        """
        Interferogram constructor, for 2-band Ifg raster datasets.

        :param str path: Path to interferogram file
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
        self.time_span = None

    def open(self, readonly=None):
        """
        Open raster file.

        :param bool readonly: True/False, or None to open as underlying file setting
        """
        RasterBase.open(self, readonly)
        self.initialize()

    def initialize(self):
        """
        Read basic interferogram properties upon opening interferogram.
        """
        self._init_dates()
        md = self.dataset.GetMetadata()
        self.wavelength = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        self.meta_data = md
        self.nan_converted = False # This flag set True after NaN conversion

    def _init_dates(self):
        """
        Determine master and slave dates, and interferogram timespan
        """
        def _to_date(datestr):
            year, month, day = [int(i) for i in datestr.split('-')]
            return date(year, month, day)

        md = self.dataset.GetMetadata()
        datestrs = [md[k] for k in [ifc.MASTER_DATE, ifc.SLAVE_DATE]]

        if all(datestrs):
            self.master, self.slave = [_to_date(s) for s in datestrs]
            self.time_span = (self.slave - self.master).days/ifc.DAYS_PER_YEAR
        else:
            msg = 'Missing master and/or slave date in %s' % self.data_path
            raise IfgException(msg)

    def convert_to_nans(self):
        """
        Convert phase data of given value to NaN
        """
        if (self._nodata_value is None) \
                or (self.dataset is None):  # pragma: no cover
            msg = 'nodata value needs to be set for nan conversion.' \
                  'Use ifg.nodata_value = NoDataValue to set nodata_value'
            log.warning(msg)
            raise RasterException(msg)
        if ((self.dataset.GetMetadataItem(ifc.NAN_STATUS) == ifc.NAN_CONVERTED)
                or self.nan_converted):
            self.phase_data = self.phase_data
            self.nan_converted = True
            msg = '{}: ignored as previous nan ' \
                  'conversion detected'.format(self.data_path)
            log.debug(msg)
            return
        else:
            self.phase_data = where(
                isclose(self.phase_data, self._nodata_value, atol=1e-6),
                nan,
                self.phase_data)
            self.meta_data[ifc.NAN_STATUS] = ifc.NAN_CONVERTED
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
        Determine the no-data value in phase band.
        """
        return self._nodata_value

    @nodata_value.setter
    def nodata_value(self, val):
        """
        Set the no-data value for phase band.
        """
        self._nodata_value = val

    @property
    def phase_data(self):
        """
        Returns phase band as an array.
        """
        if self._phase_data is None:
            self._phase_data = self.phase_band.ReadAsArray()
        return self._phase_data

    def convert_to_mm(self):
        """
        Convert phase data units from radians to millimetres.
        """
        self.mm_converted = True
        if self.dataset.GetMetadataItem(ifc.DATA_UNITS) == MILLIMETRES:
            msg = '{}: ignored as previous phase unit conversion ' \
                  'already applied'.format(self.data_path)
            log.debug(msg)
            self.phase_data = self.phase_data
            return
        elif self.dataset.GetMetadataItem(ifc.DATA_UNITS) == RADIANS:
            self.phase_data = convert_radians_to_mm(self.phase_data,
                                                    self.wavelength)
            self.meta_data[ifc.DATA_UNITS] = MILLIMETRES
            # self.write_modified_phase()
            # otherwise NaN's don't write to bytecode properly
            # and numpy complains
            # self.dataset.FlushCache()
            msg = '{}: converted phase units to millimetres'.format(self.data_path)
            log.debug(msg)
        else:  # pragma: no cover
            msg = 'Phase units are not millimetres or radians'
            raise IfgException(msg)

    @phase_data.setter
    def phase_data(self, data):
        """
        Set phase data value
        """
        self._phase_data = data

    @property
    def phase_rows(self):
        """
        Generator returning each row of the phase data.
        """
        for y in range(self.nrows):
            r = self.phase_band.ReadAsArray(yoff=y,
                                            win_xsize=self.ncols, win_ysize=1)
            yield r[0] # squeezes row from (1, WIDTH) to 1D array

    @property
    def nan_count(self):
        """
        Returns total number of NaN cells in the phase data.
        """
        return nsum(isnan(self.phase_data))

    @property
    def nan_fraction(self):
        """
        Returns decimal fraction of NaN cells in the phase band.
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

    def write_modified_phase(self, data=None):
        """
        Write updated phase data to file on disk.
        """
        if self.is_read_only:
            raise IOError("Cannot write to read only Ifg")

        # keep this block
        # if new_data_path is None:
        #     self.dataset = gdal.Open(self.data_path, GA_Update)
        # else:
        #     self.dataset = gdal.Open(new_data_path, GA_Update)
        # self._phase_band = None

        if data is not None:
            assert isinstance(data, np.ndarray)
            data_r, data_c = data.shape
            assert data_r == self.nrows and data_c == self.ncols
            self.phase_data = data
        self.phase_band.WriteArray(self.phase_data)
        for k, v in self.meta_data.items():
            self.dataset.SetMetadataItem(k, v)
        self.dataset.FlushCache()


class IfgPart(object):
    """
    Create a tile (subset) of an Ifg data object
    """
    def __init__(self, ifg_or_path, tile, ifg_dict=None, params=None):
        """
        Interferogram tile constructor.

        :param str path: Path to interferogram file
        """
        self.tile = tile
        self.r_start = self.tile.top_left_y
        self.r_end = self.tile.bottom_right_y
        self.c_start = self.tile.top_left_x
        self.c_end = self.tile.bottom_right_x
        if ifg_dict is not None:  # should be used with MPI
            ifg = ifg_dict[ifg_or_path]
            self.nan_fraction = ifg.nan_fraction
            self.master = ifg.master
            self.slave = ifg.slave
            self.time_span = ifg.time_span
            phase_file = 'phase_data_{}_{}.npy'.format(basename(ifg_or_path).split('.')[0], tile.index)
            self.phase_data = np.load(join(params[cf.TMPDIR], phase_file))
        else:
            # check if Ifg was sent.
            if isinstance(ifg_or_path, Ifg):
                ifg = ifg_or_path
            else:
                self.data_path = ifg_or_path
                ifg = Ifg(ifg_or_path)
            self.phase_data = None
            self.nan_fraction = None
            self.master = None
            self.slave = None
            self.time_span = None
        if isinstance(ifg, Ifg):
            self.read_required(ifg)

    def read_required(self, ifg):
        """
        Read interferogram file if not already open.
        """
        if not ifg.is_open:
            ifg.open(readonly=True)
        ifg.nodata_value = 0
        self.phase_data = ifg.phase_data[self.r_start:self.r_end,
                                         self.c_start:self.c_end]
        self.nan_fraction = ifg.nan_fraction
        self.master = ifg.master
        self.slave = ifg.slave
        self.time_span = ifg.time_span
        ifg.phase_data = None
        ifg.close()  # close base ifg

    @property
    def nrows(self):
        """
        Determine number of rows in tile.
        """
        return self.r_end - self.r_start

    @property
    def ncols(self):
        """
        Determine number of columns in tile.
        """
        return self.c_end - self.c_start


class Incidence(RasterBase):   # pragma: no cover
    """
    Class for storing viewing geometry data.
    e.g. incidence and azimuth raster values
    """

    def __init__(self, path):
        """
        Incidence object constructor.
        """
        RasterBase.__init__(self, path)
        self._incidence_band = None
        self._azimuth_band = None
        self._incidence_data = None
        self._azimuth_data = None

    @property
    def incidence_band(self):
        """
        Returns the GDALBand for the incidence angle layer.
        """

        if self._incidence_band is None:
            self._incidence_band = self._get_band(1)
        return self._incidence_band

    @property
    def incidence_data(self):
        """
        Returns the incidence band as an array.
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
        Returns the azimuth band as an array.
        """
        if self._azimuth_data is None:
            self._azimuth_data = self.azimuth_band.ReadAsArray()
        return self._azimuth_data


class DEM(RasterBase):
    """
    Generic raster class for single band DEM files.
    """

    def __init__(self, path):
        """
        DEM constructor.
        """
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


class RasterException(Exception):
    """
    Generic exception for raster errors.
    """


class EpochList(object):
    """
    Metadata container for epoch related information.
    """

    def __init__(self, dates=None, repeat=None, spans=None):
        """
        Construct epochlist object
        """
        self.dates = dates # list of unique dates from all the ifgs
        self.repeat = repeat
        self.spans = spans  # time span from earliest ifg

    def __str__(self):
        return "EpochList: %s" % str(self.dates)

    def __repr__(self):
        return "EpochList: %s" % repr(self.dates)


def convert_radians_to_mm(data, wavelength):
    """
    Function to translates phase in units of radians to units in millimetres.

    :param ndarray data: Interferogram phase data array
    :param float wavelength: Radar wavelength in metres

    :return: data: converted phase data
    :rtype: ndarray
    """
    return data * ifc.MM_PER_METRE * (wavelength / (4 * math.pi))


def nanmedian(x):
    """
    Determine the median of values excluding nan values.
    Use different numpy algorithm dependent on numpy version.

    :param ndarray x: array of numeric data.

    :return: y: median value
    :rtype: float
    """
    # pylint: disable=no-member
    version = [int(i) for i in pkg_resources.get_distribution("numpy").version.split('.')[:2]]
    if version[0] == 1 and version[1] > 9:
        return np.nanmedian(x)
    else:   # pragma: no cover
        return np.median(x[~np.isnan(x)])


def _is_interferogram(hdr):
    """
    Convenience function to determine if file is interferogram
    """
    return ifc.PYRATE_WAVELENGTH_METRES in hdr


def _is_incidence(hdr):
    """
    Convenience function to determine if incidence file
    """
    return 'FILE_TYPE' in hdr


def write_fullres_geotiff(header, data_path, dest, nodata):
    # pylint: disable=too-many-statements
    """
    Creates a copy of input image data (interferograms, DEM, incidence maps
    etc) in GeoTIFF format with PyRate metadata.

    :param dict header: Interferogram metadata dictionary
    :param str data_path: Input file
    :param str dest: Output destination file
    :param float nodata: No-data value

    :return: None, file saved to disk
    """
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-locals
    ifg_proc = header[ifc.PYRATE_INSAR_PROCESSOR]
    ncols = header[ifc.PYRATE_NCOLS]
    nrows = header[ifc.PYRATE_NROWS]
    bytes_per_col, fmtstr = _data_format(ifg_proc, _is_interferogram(header), ncols)
    if _is_interferogram(header) and ifg_proc == ROIPAC:
        # roipac ifg has 2 bands
        _check_raw_data(bytes_per_col*2, data_path, ncols, nrows)
    else:
        _check_raw_data(bytes_per_col, data_path, ncols, nrows)

    _check_pixel_res_mismatch(header)

    # position and projection data
    gt = [header[ifc.PYRATE_LONG], header[ifc.PYRATE_X_STEP], 0, header[ifc.PYRATE_LAT], 0, header[ifc.PYRATE_Y_STEP]]
    srs = osr.SpatialReference()
    res = srs.SetWellKnownGeogCS(header[ifc.PYRATE_DATUM])
    if res:
        msg = 'Unrecognised projection: %s' % header[ifc.PYRATE_DATUM]
        raise GeotiffException(msg)

    wkt = srs.ExportToWkt()
    dtype = 'float32' if (_is_interferogram(header) or _is_incidence(header)) else 'int16'

    # get subset of metadata relevant to PyRate
    md = collate_metadata(header)

    # create GDAL object
    ds = gdal_dataset(
        dest, ncols, nrows, driver="GTiff", bands=1, dtype=dtype, metadata=md, crs=wkt, geotransform=gt,
        creation_opts=["compress=packbits"]
    )

    # copy data from the binary file
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    row_bytes = ncols * bytes_per_col

    with open(data_path, 'rb') as f:
        for y in range(nrows):
            if ifg_proc == ROIPAC:
                if _is_interferogram(header):
                    f.seek(row_bytes, 1)  # skip interleaved band 1

            data = struct.unpack(fmtstr, f.read(row_bytes))
            band.WriteArray(np.array(data).reshape(1, ncols), yoff=y)

    ds = None  # manual close
    del ds


def gdal_dataset(out_fname, columns, rows, driver="GTiff", bands=1,
                 dtype='float32', metadata=None, crs=None,
                 geotransform=None, creation_opts=None):
    """
    Initialises a py-GDAL dataset object for writing image data.
    """
    if dtype == 'float32':
        gdal_dtype = gdal.GDT_Float32
    elif dtype == 'int16':
        gdal_dtype = gdal.GDT_Int16
    else:
        # assume gdal.GDT val is passed to function
        gdal_dtype = dtype

    # create output dataset
    driver = gdal.GetDriverByName(driver)
    outds = driver.Create(out_fname, columns, rows, bands, gdal_dtype, options=creation_opts)

    # geospatial info
    outds.SetGeoTransform(geotransform)
    outds.SetProjection(crs)

    # add metadata
    if metadata is not None:
        for k, v in metadata.items():
            outds.SetMetadataItem(k, str(v))

    return outds


def collate_metadata(header):
    """
    Grab metadata relevant to PyRate from input metadata

    :param dict header: Input file metadata dictionary

    :return: dict of relevant metadata for PyRate
    """
    md = dict()
    if _is_interferogram(header):
        for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                  ifc.PYRATE_INSAR_PROCESSOR,
                  ifc.MASTER_DATE, ifc.SLAVE_DATE,
                  ifc.DATA_UNITS, ifc.DATA_TYPE]:
            md.update({k: str(header[k])})
        if header[ifc.PYRATE_INSAR_PROCESSOR] == GAMMA:
            for k in [ifc.MASTER_TIME, ifc.SLAVE_TIME, ifc.PYRATE_INCIDENCE_DEGREES]:
                md.update({k: str(header[k])})
    elif _is_incidence(header):
        md.update({ifc.DATA_TYPE:ifc.INCIDENCE})
    else: # must be dem
        md.update({ifc.DATA_TYPE:ifc.DEM})

    return md


def _data_format(ifg_proc, is_ifg, ncols):
    """
    Convenience function to determine the bytesize and format of input files
    """
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
    else:  # pragma: no cover
        msg = 'Unrecognised InSAR Processor: %s' % ifg_proc
        raise GeotiffException(msg)
    return bytes_per_col, fmtstr


def _check_raw_data(bytes_per_col, data_path, ncols, nrows):
    """
    Convenience function to check the file size is as expected
    """
    size = ncols * nrows * bytes_per_col
    act_size = os.stat(data_path).st_size
    if act_size != size:
        msg = '%s should have size %s, not %s. Is the correct file being used?'
        raise GeotiffException(msg % (data_path, size, act_size))


def _check_pixel_res_mismatch(header):
    """
    Convenience function to check equality of pixel resolution in X and Y dimensions
    """
    # pylint: disable=invalid-name
    xs, ys = [abs(i) for i in [header[ifc.PYRATE_X_STEP], header[ifc.PYRATE_Y_STEP]]]

    if xs != ys:
        msg = 'X and Y cell sizes do not match: %s & %s'
        raise GeotiffException(msg % (xs, ys))


def write_unw_from_data_or_geotiff(geotif_or_data, dest_unw, ifg_proc):
    """
    Function to write numpy array data or a geotiff to a GAMMA-format
    big-endian float32 .unw file

    :param str geotif_or_data: path name of geotiff file to convert
        OR
    :param ndarray geotif_or_data: numpy array of data to convert
    :param str dest_unw: destination unw file
    :param int ifg_proc: processor type, GAMMA=1, ROIPAC=0

    :return: None, file saved to disk
    """
    if ifg_proc != 1:
        raise NotImplementedError('only supports GAMMA format for now')
    if isinstance(geotif_or_data, str):
        assert os.path.exists(geotif_or_data), 'make sure geotif exists'
        ds = gdal.Open(geotif_or_data)
        data = ds.ReadAsArray()
        ds = None
    else:
        data = geotif_or_data

    nrows, ncols = data.shape
    fmtstr = '!' + ('f' * ncols)  # data format is big endian float32s

    with open(dest_unw, 'wb') as f:
        for y in range(nrows):
            col_data = struct.pack(fmtstr, *data[y, :])
            f.write(col_data)


# This function may be able to be deprecated
def write_output_geotiff(md, gt, wkt, data, dest, nodata):
    # pylint: disable=too-many-arguments
    """
    Writes PyRate output data to a GeoTIFF file.

    :param dict md: Dictionary containing PyRate metadata
    :param list gt: GDAL geotransform for the data
    :param list wkt: GDAL projection information for the data
    :param ndarray data: Output data array to save
    :param str dest: Destination file name
    :param float nodata: No data value of data

    :return None, file saved to disk
    """

    driver = gdal.GetDriverByName("GTiff")
    nrows, ncols = data.shape
    ds = driver.Create(dest, ncols, nrows, 1, gdal.GDT_Float32, options=['compress=packbits'])
    # set spatial reference for geotiff
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    ds.SetMetadataItem(ifc.EPOCH_DATE, str(md[ifc.EPOCH_DATE]))

    # set other metadata
    ds.SetMetadataItem('DATA_TYPE', str(md['DATA_TYPE']))
    # sequence position for time series products
    if "SEQUENCE_POSITION" in md:
        ds.SetMetadataItem("SEQUENCE_POSITION", str(md["SEQUENCE_POSITION"]))
    if "PYRATE_REFPIX_LAT" in md:
        ds.SetMetadataItem("PYRATE_REFPIX_LAT", str(md["PYRATE_REFPIX_LAT"]))
    if "PYRATE_REFPIX_LON" in md:
        ds.SetMetadataItem("PYRATE_REFPIX_LON", str(md["PYRATE_REFPIX_LON"]))
    if "PYRATE_REFPIX_X" in md:
        ds.SetMetadataItem("PYRATE_REFPIX_X", str(md["PYRATE_REFPIX_X"]))
    if "PYRATE_REFPIX_Y" in md:
        ds.SetMetadataItem("PYRATE_REFPIX_Y", str(md["PYRATE_REFPIX_Y"]))
    if "PYRATE_MEAN_REF_AREA" in md:
        ds.SetMetadataItem("PYRATE_MEAN_REF_AREA", str(md["PYRATE_MEAN_REF_AREA"]))
    if "STANDARD_DEVIATION_REF_AREA" in md:
        ds.SetMetadataItem("PYRATE_STANDARD_DEVIATION_REF_AREA", str(md["PYRATE_STANDARD_DEVIATION_REF_AREA"]))

    # write data to geotiff
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)


def write_geotiff(data, outds, nodata):
    # pylint: disable=too-many-arguments
    """
    A generic routine for writing a NumPy array to a geotiff.

    :param ndarray data: Output data array to save
    :param obj outds: GDAL destination object
    :param float nodata: No data value of data

    :return None, file saved to disk
    """
    # only support "2 <= dims <= 3"
    if data.ndim == 3:
        count, height, width = data.shape
    elif data.ndim == 2:
        height, width = data.shape
    else:
        msg = "Only support dimensions of '2 <= dims <= 3'."
        raise GeotiffException(msg)

    # write data to geotiff
    band = outds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)

    outds = None
    band = None
    del outds
    del band


class GeotiffException(Exception):
    """
    Geotiff exception class
    """

def create_tiles(shape, nrows=2, ncols=2):
    """
    Return a list of tiles containing nrows x ncols with each tile preserving
    the physical layout of original array. The number of rows can be changed
    (increased) such that the resulting tiles with float32's do not exceed
    500MB in memory. When the array shape (rows, columns) are not divisible
    by (nrows, ncols) then some of the array dimensions can change according
    to numpy.array_split.

    :param tuple shape: Shape tuple (2-element) of interferogram.
    :param int nrows: Number of rows of tiles
    :param int ncols: Number of columns of tiles

    :return: List of Tile class instances.
    :rtype: list
    """

    if len(shape) != 2:
        raise ValueError('shape must be a length 2 tuple')

    no_y, no_x = shape

    if ncols > no_x or nrows > no_y:
        raise ValueError('nrows/cols must be greater than ifg dimensions')
    col_arr = np.array_split(range(no_x), ncols)
    row_arr = np.array_split(range(no_y), nrows)
    return [Tile(i, (r[0], c[0]), (r[-1]+1, c[-1]+1)) for i, (r, c) in enumerate(product(row_arr, col_arr))]


def get_tiles(ifg_path, rows, cols):
    """
    Break up the interferograms into smaller tiles based on user supplied
    rows and columns.

    :param list ifg_path: List of destination geotiff file names
    :param int rows: Number of rows to break each interferogram into
    :param int cols: Number of columns to break each interferogram into

    :return: tiles: List of shared.Tile instances
    :rtype: list
    """
    ifg = Ifg(ifg_path)
    ifg.open(readonly=True)
    tiles = create_tiles(ifg.shape, nrows=rows, ncols=cols)
    ifg.close()
    return tiles


class Tile():
    """
    Tile class for containing a sub-part of an interferogram
    """
    def __init__(self, index, top_left, bottom_right):
        """
        Parameters
        ----------
        index: int
            identifying index of a tile
        top_left: tuple
            ifg index of top left of tile
        bottom_right: tuple
            ifg index of bottom right of tile
        """

        self.index = index
        self.top_left = top_left
        self.bottom_right = bottom_right
        self.top_left_y, self.top_left_x = top_left
        self.bottom_right_y, self.bottom_right_x = bottom_right

    def __str__(self):
        return "Convenience Tile class containing tile co-ordinates"


def nan_and_mm_convert(ifg, params):
    """
    Perform millimetre and nan conversion on interferogram data

    :param Ifg instance ifg: Interferogram class instance
    :param dict params: Dictionary of parameters

    :return: None, data modified internally
    """
    nan_conversion = params[cf.NAN_CONVERSION]
    if nan_conversion:  # nan conversion happens here in networkx mst
        # if not ifg.nan_converted:
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        ifg.convert_to_nans()
    if not ifg.mm_converted:
        ifg.convert_to_mm()


def cell_size(lat, lon, x_step, y_step):
    # pylint: disable=invalid-name
    """
    Converts X|Y_STEP in degrees to X & Y cell size in metres.
    This function depends on PyProj/PROJ4 to implement the function

    :param float lat: Latitude in degrees
    :param float lon: Longitude in degrees
    :param float x_step: Horizontal step size in degrees
    :param float y_step: Vertical step size in degrees

    :return: tuple of X and Y cell size floats
    :rtype: tuple
    """
    if lat > 84.0 or lat < -80:
        msg = "No UTM zone for polar region: > 84 degrees N or < 80 degrees S. " \
              "Provided values are lat: "+str(lat) +" long: " +str(lon)
        raise ValueError(msg)

    zone = _utm_zone(lon)
    p0 = pyproj.Proj(proj='latlong', ellps='WGS84')
    p1 = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')

    x0, y0 = pyproj.transform(p0, p1, lon, lat, errcheck=True)
    x1, y1 = pyproj.transform(p0, p1, lon + x_step, lat + y_step, errcheck=True)

    return tuple(abs(e) for e in (x1 - x0, y1 - y0))


def _utm_zone(longitude):
    """
    Returns basic UTM zone for given longitude in degrees.
    Currently does NOT handle the sub-zoning around Scandanavian countries.
    See http://www.dmap.co.uk/utmworld.htm
    """
    if longitude == 180:
        return 60.0
    return floor((longitude + 180) / 6.0) + 1


class PrereadIfg():
    """
    Convenience class for handling pre-calculated ifg params
    """
    # pylint: disable=too-many-arguments
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path, nan_fraction, master, slave, time_span,
                 nrows, ncols, metadata):
        self.path = path
        self.nan_fraction = nan_fraction
        self.master = master
        self.slave = slave
        self.time_span = time_span
        self.nrows = nrows
        self.ncols = ncols
        self.shape = (nrows, ncols)
        self.metadata = metadata


def _prep_ifg(ifg_path, params):
    """
    Wrapper for reading an interferogram file and creating an Ifg object

    :param str ifg_path: Interferogram file path
    :param dict params: Dictionary of configuration parameters

    :return: ifg: Interferogram class instance
    :rtype: xxxx (eg flaot)
    """
    # Only used in pyrate.scripts.run_pyrate?
    ifg = Ifg(ifg_path)
    ifg.open()
    nan_and_mm_convert(ifg, params)
    return ifg


def save_numpy_phase(ifg_paths, tiles, params):
    """
    Save interferogram phase data as numpy array file on disk.

    :param list ifg_paths: List of strings for interferogram paths
    :param list tiles: List of pyrate.shared.Tile instances
    :param dict params: Dictionary of configuration parameters

    :return: None, file saved to disk
    """
    process_ifgs = mpiops.array_split(ifg_paths)
    outdir = params[cf.TMPDIR]
    if not os.path.exists(outdir):
        mkdir_p(outdir)
    for ifg_path in process_ifgs:
        ifg = Ifg(ifg_path)
        ifg.open()
        phase_data = ifg.phase_data
        bname = basename(ifg_path).split('.')[0]
        for t in tiles:
            p_data = phase_data[t.top_left_y:t.bottom_right_y,
                                t.top_left_x:t.bottom_right_x]
            phase_file = 'phase_data_{}_{}.npy'.format(bname, t.index)
            np.save(file=join(outdir, phase_file),
                    arr=p_data)
        ifg.close()
    mpiops.comm.barrier()


def get_geotiff_header_info(ifg_path):
    """
    Return information from a geotiff interferogram header using GDAL methods.

    :param str ifg_path: path to interferogram geotiff file

    :return: md: PyRate metadata
    :rtype: list
    :return: gt: GDAL geotransform for the data
    :rtype: list
    :return: wkt: GDAL projection information for the data
    :rtype: list
    """
    ds = gdal.Open(ifg_path)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    gt = ds.GetGeoTransform()  # get geographical bounds of data
    wkt = ds.GetProjection()  # get projection of data
    ds = None  # close dataset
    return gt, md, wkt


def warp_required(xlooks, ylooks, crop):
    """
    Check if a crop or multi-look operation is required.

    :param int xlooks: Resampling/multi-looking in x dir
    :param int ylooks: Resampling/multilooking in y dir
    :param int crop: Interferogram crop option

    :return: True if params show rasters need to be cropped and/or resized
    :rtype: bool
    """
    if xlooks > 1 or ylooks > 1:
        return True
    if crop is None:
        return False
    return True


def output_tiff_filename(inpath, outpath):
    """
    Output geotiff filename for a given input filename.

    :param str inpath: path of input file location
    :param str outpath: path of output file location

    :return: Geotiff filename for the given file.
    :rtype: str
    """
    fname, ext = os.path.basename(inpath).split('.')
    outpath = os.path.dirname(inpath) if outpath is None else outpath
    if ext == 'tif':
        name = os.path.join(outpath, fname + '.tif')
    else:
        name = os.path.join(outpath, fname + '_' + ext + '.tif')
    return name


def check_correction_status(ifgs, meta):  # pragma: no cover
    """
    Generic function for checking if a correction has already been performed
    in a previous run by interrogating PyRate meta data entries

    :param preread_ifgs: Dictionary of pre-read interferogram information
    :param str meta: Meta data flag to check for

    :return: True if correction has been performed, otherwise False
    :rtype: bool
    """
    def close_all(ifgs):    
        for ifg in ifgs:
            ifg.close()
    
    if not isinstance(ifgs[0], Ifg):
        ifgs = [Ifg(ifg_path) for ifg_path in ifgs]

    for ifg in ifgs:
        if not ifg.is_open:
            ifg.open()

    flags = [meta in ifg.meta_data for ifg in ifgs]
    if all(flags):
        log.info('Skipped: interferograms already corrected')
        return True
    elif not all(flags) and any(flags):
        log.debug('Detected mix of corrected and uncorrected interferograms')
        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '{}: correction detected'.format(i.data_path)
            else:
                msg = '{}: correction NOT detected'.format(i.data_path)
            log.debug(msg)
            close_all(ifgs)
            raise CorrectionStatusError(msg)
    else:
        log.debug('Calculating corrections')
        close_all(ifgs)
        return False


class CorrectionStatusError(Exception):
    """
    Generic class for correction status errors.
    """


def extract_epochs_from_filename(filename_with_epochs: str) -> List[str]:
    src_epochs = re.findall(r"(\d{8})", str(filename_with_epochs))
    if not len(src_epochs) > 0:
        src_epochs = re.findall(r"(\d{6})", str(filename_with_epochs))
    return src_epochs


def mpi_vs_multiprocess_logging(step, params):
    if mpiops.size > 1:  # Over-ride input options if this is an MPI job
        log.info(f"Running {step} step using mpi processing. Disabling parallel processing.")
        params[cf.PARALLEL] = 0
    else:
        if params[cf.PARALLEL] == 1:
            log.info(f"Running {step} using {params[cf.PROCESSES]} processes")
        else:
            log.info(f"Running {step} serially")


def dem_or_ifg(data_path):
    """
    Returns an Ifg or DEM class object from input geotiff file.

    :param str data_path: file path name

    :return: Interferogram or DEM object from input file
    :rtype: Ifg or DEM class object
    """
    ds = gdal.Open(data_path)
    md = ds.GetMetadata()
    if ifc.MASTER_DATE in md:  # ifg
        return Ifg(data_path)
    else:
        return DEM(data_path)