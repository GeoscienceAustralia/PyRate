#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
"""
This Python module contains utilities and classes shared by
all other PyRate modules
"""
# pylint: disable=too-many-lines
from __future__ import print_function
import errno
import logging
import math
from math import floor
import os
from os.path import basename, dirname, join
import shutil
import stat
import struct
from datetime import date
from itertools import product
import numpy as np
from numpy import where, nan, isnan, sum as nsum, isclose
import pyproj
import pkg_resources

from pyrate import ifgconstants as ifc, mpiops
from pyrate import roipac, gamma, config as cf

VERBOSE = True
log = logging.getLogger(__name__)

try:
    from osgeo import osr, gdal
    from osgeo.gdalconst import GA_Update, GA_ReadOnly
except ImportError:
    import gdal

gdal.UseExceptions()

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


def mkdir_p(path):
    """ copied from stackoverflow"""
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
            raise IOError('The file {path} does not exist. Consider running '
                          'prepifg'.format(path=self.data_path))

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
        self.x_centre = int(self.ncols / 2)
        self.y_centre = int(self.nrows / 2)
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
    Interferrogram class, represents the difference between two acquisitions.
    Ifg objects double as a container for related data.
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path):
        """
        Interferogram constructor, for 2 band ROIPAC Ifg raster datasets.
        Parameters
        ----------
        path: str
            path to ifg
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
        :param bool readonly: True/False, or None to open as underlying file setting
        """
        RasterBase.open(self, readonly)
        self.initialize()

    def initialize(self):
        """basic ifg properties read on opening Ifg"""
        self._init_dates()
        md = self.dataset.GetMetadata()
        self.wavelength = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        self.meta_data = md
        # creating code needs to set this flag after 0 -> NaN replacement
        self.nan_converted = False

    def _init_dates(self):
        """grab master and slade date from metadata"""
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
        Converts given values in phase data to NaNs
        :param val: value to convert, default is 0
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
            msg = '{}: converted phase units ' \
                  'to millimetres'.format(self.data_path)
            log.debug(msg)
        else:  # pragma: no cover
            msg = 'Phase units are not millimetres or radians'
            raise IfgException(msg)

    @phase_data.setter
    def phase_data(self, data):
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

    def write_modified_phase(self, data=None):
        """
        Writes phase data to disk.
        For this to work, a copy of the original file
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

    def save_numpy_phase(self, numpy_file):
        """
        Parameters
        ----------
        numpy_file: str
            file path where phase data is saved
        """
        np.save(file=numpy_file, arr=self.phase_data)


class IfgPart(object):
    """
    slice of Ifg data object
    """
    # pylint: disable=missing-docstring
    # pylint: disable=too-many-instance-attributes
    def __init__(self, ifg_or_path, tile, ifg_dict=None):

        self.tile = tile
        self.r_start = self.tile.top_left_y
        self.r_end = self.tile.bottom_right_y
        self.c_start = self.tile.top_left_x
        self.c_end = self.tile.bottom_right_x
        # TODO: fix this if cond
        if ifg_dict is not None:  # should be used with MPI
            ifg = ifg_dict[ifg_or_path]
            self.nan_fraction = ifg.nan_fraction
            self.master = ifg.master
            self.slave = ifg.slave
            self.time_span = ifg.time_span
            phase_file = 'phase_data_{}_{}.npy'.format(
                basename(ifg_or_path).split('.')[0], tile.index)
            self.phase_data = np.load(join(dirname(ifg_or_path), phase_file))
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
        return self.r_end - self.r_start

    @property
    def ncols(self):
        return self.c_end - self.c_start


class Incidence(RasterBase):   # pragma: no cover
    """ Incidence class"""

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
        """DEM constructor."""
        RasterBase.__init__(self, path)
        self._band = None

    @property
    def height_band(self):
        """Returns the GDALBand for the elevation layer."""

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
    return data * ifc.MM_PER_METRE * (wavelength / (4 * math.pi))


def nanmedian(x):
    """
    :param x:
    :return:
    """
    # pylint: disable=no-member
    version = [int(i) for i in
               pkg_resources.get_distribution("numpy").version.split('.')]
    if version[0] == 1 and version[1] > 9:
        return np.nanmedian(x)
    else:   # pragma: no cover
        return np.median(x[~np.isnan(x)])


def write_geotiff(header, data_path, dest, nodata):
    # pylint: disable=too-many-statements
    """
    Writes input image data (interferograms, DEM, incidence maps etc)
    to GeoTIFF format with PyRate metadata
    """
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-locals
    is_ifg = ifc.PYRATE_WAVELENGTH_METRES in header
    is_incidence = 'FILE_TYPE' in header
    ifg_proc = header[ifc.PYRATE_INSAR_PROCESSOR]
    ncols = header[ifc.PYRATE_NCOLS]
    nrows = header[ifc.PYRATE_NROWS]

    # need to have gamma/roipac functionality here?
    if ifg_proc == ROIPAC:
        roipac.check_raw_data(is_ifg, data_path, ncols, nrows)
        roipac.check_step_mismatch(header)
    else:  # GAMMA
        gamma.check_raw_data(data_path, ncols, nrows)
        gamma.check_step_mismatch(header)

    driver = gdal.GetDriverByName("GTiff")
    dtype = gdal.GDT_Float32 if (is_ifg or is_incidence) else gdal.GDT_Int16
    ds = driver.Create(dest, ncols, nrows, 1, dtype)

    # write pyrate parameters to headers
    if is_ifg:
        for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                  ifc.PYRATE_INSAR_PROCESSOR,
                  ifc.MASTER_DATE, ifc.SLAVE_DATE,
                  ifc.DATA_UNITS, ifc.DATA_TYPE]:
            ds.SetMetadataItem(k, str(header[k]))
        if ifg_proc == GAMMA:
            for k in [ifc.MASTER_TIME, ifc.SLAVE_TIME, ifc.PYRATE_INCIDENCE_DEGREES]:
                ds.SetMetadataItem(k, str(header[k]))
    elif is_incidence:
        ds.SetMetadataItem(ifc.DATA_TYPE, ifc.INCIDENCE)
    else: # must be dem
        ds.SetMetadataItem(ifc.DATA_TYPE, ifc.DEM)

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
    else:  # pragma: no cover
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


def write_unw_from_data_or_geotiff(geotif_or_data, dest_unw, ifg_proc):
    """
    :param geotif_or_data: data or geotif to covert into unw
    :param dest_unw: destination unw file
    :param ifg_proc: processor type, GAMMA=1, ROIPAC=0
    """
    if ifg_proc != 1:
        raise NotImplementedError('only support gamma processor for now')
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


def write_output_geotiff(md, gt, wkt, data, dest, nodata):
    # pylint: disable=too-many-arguments
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
    ds.SetMetadataItem(ifc.EPOCH_DATE, str(md[ifc.EPOCH_DATE]))

    # set other metadata
    ds.SetMetadataItem('DATA_TYPE', str(md['DATA_TYPE']))
    # sequence position for time series products
    if 'SEQUENCE_POSITION' in md:
        ds.SetMetadataItem('SEQUENCE_POSITION', str(md['SEQUENCE_POSITION']))

    # write data to geotiff
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)


class GeotiffException(Exception):
    """
    Geotiff exception class
    """


def create_tiles(shape, nrows=2, ncols=2):
    """
    shape must be a 2-tuple, i.e., 2d_array.shape.
    The returned list contains nrowsXncols Tiles with each tile preserving the
    "physical" layout of original arr.

    The number of rows can be changed (increased) such that  the resulting tiles
    with float32's do not exceed 500MB in memory.

    When the array shape (rows, cols) are not divisible by (nrows, ncols) then
    some of the array dimensions can change according to numpy.array_split.

    Parameters
    ----------
    shape: tuple
        shape tuple of ifg
    nrows: int
        number of rows of tiles
    ncols: int
        number of columns of tiles
    Returns
    -------
    list of Tile class instances
    """

    if len(shape) != 2:
        raise ValueError('shape must be a length 2 tuple')

    no_y, no_x = shape

    if ncols > no_x or nrows > no_y:
        raise ValueError('nrows/cols must be greater than ifg dimensions')
    col_arr = np.array_split(range(no_x), ncols)
    row_arr = np.array_split(range(no_y), nrows)
    return [Tile(i, (r[0], c[0]), (r[-1]+1, c[-1]+1))
            for i, (r, c) in enumerate(product(row_arr, col_arr))]


class Tile:
    """
    Tile class containing part of the ifgs
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


def copytree(src, dst, symlinks=False, ignore=None):
    # pylint: disable=line-too-long
    """
    copy contents of src dir into dst dir
    stolen from: http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth?lq=1
    :param src: source dir to copy from
    :param dst: dst dir to copy to, created if does not exist
    :param symlinks: bool, whether to copy symlink or not
    :param ignore:
    """
    # pylint: disable=invalid-name
    if not os.path.exists(dst):  # pragma: no cover
        os.makedirs(dst)
    shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if symlinks and os.path.islink(s):  # pragma: no cover
            if os.path.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except AttributeError:
                pass  # lchmod not available
        elif os.path.isdir(s):  # pragma: no cover
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def pre_prepare_ifgs(ifg_paths, params):
    """
    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict
    Returns
    -------
    ifgs: list
        list of Ifg instances
    """
    ifgs = [Ifg(p) for p in ifg_paths]
    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        nan_and_mm_convert(i, params)
    log.info('Opened ifg for reading')
    return ifgs


def nan_and_mm_convert(ifg, params):
    """
    Parameters
    ----------
    ifg: Ifg class instance
    params: dict
        parameters dict
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
    Collection of geodesy/pyproj algorithms for PyRate.
    This function depends on PyProj/PROJ4 to replace llh2local.m in MATLAB Pirate.
    Converts X|Y_STEP in degrees to X & Y cell length/width in metres.
    lat: latitude in degrees
    lon: longitude in degrees
    x_step: horizontal step size in degrees
    y_step: vertical step size in degrees
    """
    if lat > 84.0 or lat < -80:
        msg = "No UTM zone for polar region: > 84 degrees N or < 80 degrees S"
        raise ValueError(msg)

    zone = utm_zone(lon)
    p0 = pyproj.Proj(proj='latlong', ellps='WGS84')
    p1 = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')
    assert p0.is_latlong()
    assert not p1.is_latlong()

    x0, y0 = pyproj.transform(p0, p1, lon, lat)
    x1, y1 = pyproj.transform(p0, p1, lon + x_step, lat + y_step)
    return tuple(abs(e) for e in (x1 - x0, y1 - y0))


def utm_zone(longitude):
    """
    Returns basic UTM zone for given longitude in degrees. Currently does NOT
    handle the sub-zoning around Scandanavian countries.
    See http://www.dmap.co.uk/utmworld.htm
    """
    if longitude == 180:
        return 60.0
    return floor((longitude + 180) / 6.0) + 1


class PrereadIfg:
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


def prepare_ifg(ifg_path, params):
    """
    Parameters
    ----------
    ifg_path: str
        ifg path
    tiles: list
        list of Tile instances
    params: dict
        config dict
    Returns
    -------
    ifg: ifg class instance
    """
    ifg = Ifg(ifg_path)
    ifg.open()
    nan_and_mm_convert(ifg, params)
    return ifg


def save_numpy_phase(ifg_paths, tiles, params):
    """
    Save ifg phase data as numpy array.

    Parameters
    ----------
    ifg_paths: list
        list of strings corresponding to ifg paths
    params: dict
        config dict
    tiles: list
        list of Shared.Tile instances
    """
    process_ifgs = mpiops.array_split(ifg_paths)
    outdir = params[cf.OUT_DIR]
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


def get_projection_info(ifg_path):
    """
    return projection information of ifg

    Parameters
    ----------
    ifg_path: str
        ifg path
    """
    ds = gdal.Open(ifg_path)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    gt = ds.GetGeoTransform()  # get geographical bounds of data
    wkt = ds.GetProjection()  # get projection of data
    ds = None  # close dataset
    return gt, md, wkt


def warp_required(xlooks, ylooks, crop):
    """
    Returns True if params show rasters need to be cropped and/or resized.

    Parameters
    ----------
    xlooks: int
        resampling/multilooking in x dir

    ylooks: int
        resampling/multilooking in y dir

    crop: int
        ifg crop option
    """

    if xlooks > 1 or ylooks > 1:
        return True

    if crop is None:
        return False

    return True


def output_tiff_filename(inpath, outpath):
    """
    return the geotiff filename for a given input filename
    """
    fname, ext = os.path.basename(inpath).split('.')
    return os.path.join(outpath, fname + '_' + ext + '.tif')
