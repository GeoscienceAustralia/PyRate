#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
from typing import List, Union, Iterable, Callable

import errno
import math
from math import floor
import os
from os.path import basename, join
from pathlib import Path
import struct
from datetime import date
from itertools import product
from enum import Enum
from joblib import Parallel, delayed
import numpy as np
from numpy import where, nan, isnan, sum as nsum, isclose
import pyproj
import pkg_resources

from osgeo import osr, gdal
from osgeo.gdalconst import GA_Update, GA_ReadOnly

import pyrate.constants as C

from pyrate.core import ifgconstants as ifc, mpiops
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


class InputTypes(Enum):
    """An enum of the types of input files & directories"""
    IFG = 'ifg'
    COH = 'coh'
    BASE = 'base'
    LT = 'lt'
    DEM = 'dem'
    HEADER = 'header'
    DIR_MAP = {
        IFG: C.INTERFEROGRAM_DIR,
        COH: C.COHERENCE_DIR,
        DEM: C.GEOMETRY_DIR,
    }


def joblib_log_level(level: str) -> int:
    """
    Convert python log level to joblib int verbosity.
    """
    if level == 'INFO':
        return 0

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


class RasterBase:
    """
    Base class for PyRate GeoTIFF based raster datasets.
    """
    # pylint: disable=missing-docstring
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path: Union[gdal.Dataset, str, Path]):
        if isinstance(path, gdal.Dataset):
            self.dataset = path  # path will be Dataset in this case
            self.data_path = self.dataset  # data_path dummy
            self.add_geographic_data()
        else:
            path = path.as_posix() if isinstance(path, Path) else path
            self.data_path = path
            self.dataset = None  # for GDAL dataset obj
            self._readonly = not os.access(path, os.R_OK | os.W_OK)

            if self._readonly is None:
                raise NotImplementedError  # os.access() has failed?

    def __str__(self):
        name = self.__class__.__name__
        return f"{name}('{self.data_path}')"

    def __repr__(self):
        name = self.__class__.__name__
        return f"{name}('{self.data_path}')"

    def open(self, readonly=None):
        """
        Opens generic raster dataset.
        """
        if self.dataset is not None:
            msg = f"open() already called for {self}"
            raise RasterException(msg)

        if not os.path.exists(self.data_path):
            raise IOError(
                f'The file {self.data_path} does not exist. Consider first running prepifg'
            )

        # unless read only, by default open files as writeable
        if readonly not in [True, False, None]:
            raise ValueError("readonly must be True, False or None")

        if readonly is False and self._readonly is True:
            raise IOError("Cannot open write protected file for writing")

        flag = GA_ReadOnly if self._readonly else GA_Update
        self.dataset = gdal.Open(self.data_path, flag)
        if self.dataset is None:
            raise RasterException(f"Error opening {self.data_path}")

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
        self.x_size, self.y_size = cell_size(
            self.lat_centre, self.long_centre, self.x_step, self.y_step
        )

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
        if not self.is_open:
            raise RasterException('Dataset not open')

        return self.dataset.RasterXSize * self.dataset.RasterYSize

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
        if self.dataset is None:
            raise RasterException(f"Raster {self.data_path} has not been opened")

        return self.dataset.GetRasterBand(band)


class Ifg(RasterBase):
    """
    Interferogram (Ifg) class objects; double as a container for
    interferometric phase raster band data and related data.
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path: Union[str, Path, gdal.Dataset]):
        """
        Interferogram constructor, for 2-band Ifg raster datasets.

        :param str path: Path to interferogram file
        """
        RasterBase.__init__(self, path)
        self._phase_band = None
        self._phase_data = None
        self.first = None
        self.second = None
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
        self.nan_converted = False  # This flag set True after NaN conversion

    def _init_dates(self):
        """
        Determine first and second image dates, and interferogram timespan
        """
        def _to_date(datestr):
            year, month, day = [int(i) for i in datestr.split('-')]
            return date(year, month, day)

        md = self.dataset.GetMetadata()
        datestrs = [md[k] for k in [ifc.FIRST_DATE, ifc.SECOND_DATE]]

        if all(datestrs):
            self.first, self.second = [_to_date(s) for s in datestrs]
            self.time_span = (self.second - self.first).days/ifc.DAYS_PER_YEAR
        else:
            msg = f'Missing first and/or second date in {self.data_path}'
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
            msg = f'{self.data_path}: ignored as previous nan conversion detected'
            log.debug(msg)
            return

        self.phase_data = where(
            isclose(self.phase_data, self._nodata_value, atol=1e-6),
            nan,
            self.phase_data)
        self.meta_data[ifc.NAN_STATUS] = ifc.NAN_CONVERTED
        self.nan_converted = True

        # self.write_modified_phase(self.phase_data)

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
        # TODO: enhance this to use x/y offset and size
        if self._phase_data is None:
            self._phase_data = self.phase_band.ReadAsArray()
        return self._phase_data

    def convert_to_mm(self):
        """
        Convert phase_data units from radians to millimetres.
        Note: converted phase_data held in memory and not written to disc
        (see shared.write_modified_phase)
        """
        if self.dataset.GetMetadataItem(ifc.DATA_UNITS) == MILLIMETRES:
            self.mm_converted = True
            msg = f'{self.data_path}: ignored as phase units are already millimetres'
            log.debug(msg)
        elif self.dataset.GetMetadataItem(ifc.DATA_UNITS) == RADIANS:
            self.phase_data = convert_radians_to_mm(self.phase_data, self.wavelength)
            self.meta_data[ifc.DATA_UNITS] = MILLIMETRES
            self.mm_converted = True
            # self.write_modified_phase()
            # otherwise NaN's don't write to bytecode properly
            # and numpy complains
            # self.dataset.FlushCache()
            msg = f'{self.data_path}: converted phase units to millimetres'
            log.debug(msg)
        else:  # pragma: no cover
            msg = 'Phase units are not millimetres or radians'
            raise IfgException(msg)

    def convert_to_radians(self):
        """
        Convert phase_data units from millimetres to radians.
        Note: converted phase_data held in memory and not written to disc
        (see shared.write_modified_phase)
        """
        if self.meta_data[ifc.DATA_UNITS] == MILLIMETRES:
            self.phase_data = convert_mm_to_radians(self.phase_data, wavelength=self.wavelength)
            self.meta_data[ifc.DATA_UNITS] = RADIANS
            self.mm_converted = False
            msg = f'{self.data_path}: converted phase units to radians'
            log.debug(msg)
        elif self.meta_data[ifc.DATA_UNITS] == RADIANS:
            self.mm_converted = False
            msg = f'{self.data_path}: ignored as phase units are already radians'
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

    def add_metadata(self, **kwargs):
        """Adds metadata to the interferogram's geotiff file"""
        if (not self.is_open) or self.is_read_only:
            raise IOError("Ifg not open or readonly. Cannot write!")

        for k, v in kwargs.items():
            self.dataset.SetMetadataItem(k, v)
        self.dataset.FlushCache()  # write to disc


class Tile:
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


class IfgPart:
    """
    Create a tile (subset) of an Ifg data object
    """
    def __init__(self, ifg_or_path, tile: Tile, ifg_dict=None, params=None):
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
            self.first = ifg.first
            self.second = ifg.second
            self.time_span = ifg.time_span
            phase_file = f'phase_data_{basename(ifg_or_path).split(".")[0]}_{tile.index}.npy'
            self.phase_data = np.load(join(params[C.TMPDIR], phase_file))
        else:
            # check if Ifg was sent.
            if isinstance(ifg_or_path, Ifg):
                ifg = ifg_or_path
            else:
                self.data_path = ifg_or_path
                ifg = Ifg(ifg_or_path)
            self.phase_data = None
            self.nan_fraction = None
            self.first = None
            self.second = None
            self.time_span = None
        if isinstance(ifg, Ifg):
            self.read_tile(ifg)

    def read_tile(self, ifg: Ifg):
        """
        Read interferogram file if not already open.
        """
        if not ifg.is_open:
            ifg.open(readonly=True)
        ifg.nodata_value = 0
        self.phase_data = ifg.phase_data[self.r_start:self.r_end,
                                         self.c_start:self.c_end]
        self.nan_fraction = ifg.nan_fraction
        self.first = ifg.first
        self.second = ifg.second
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


class TileMixin:
    """A mixin class that makes the inheriting class callable to subscript by a tile region"""
    def __call__(self, tile: Tile):
        t = tile
        return self.data[t.top_left_y:t.bottom_right_y, t.top_left_x:t.bottom_right_x]


class DEM(RasterBase, TileMixin):
    """
    Generic raster class for single band DEM/Geometry files.
    """

    def __init__(self, path):
        """
        DEM constructor.
        """
        RasterBase.__init__(self, path)
        self._band = None
        self._data = None

    @property
    def band(self):
        """
        Returns the GDALBand for the elevation layer.
        """
        if not self.is_open:
            self.open()
        if self._band is None:
            self._band = self._get_band(1)
        return self._band

    @property
    def data(self):
        """
        Returns the geometry band as an array.
        """
        if self._data is None:
            self._data = self.band.ReadAsArray()
        return self._data


class MemGeometry(TileMixin):
    """A simple class which holds geometry data in memory"""
    def __init__(self, data):
        """
        Set phase data value
        """
        self._data = data

    @property
    def data(self):
        """
        Returns the geometry band as an array.
        """
        return self._data


Geometry = DEM


class IfgException(Exception):
    """
    Generic exception class for interferogram errors.
    """


class RasterException(Exception):
    """
    Generic exception for raster errors.
    """


class EpochList:
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
        return f"EpochList: {str(self.dates)}"

    def __repr__(self):
        return f"EpochList: {repr(self.dates)}"


def convert_radians_to_mm(data, wavelength):
    """
    Function to translates phase in units of radians to units in millimetres.

    :param ndarray data: Interferogram phase data array
    :param float wavelength: Radar wavelength in metres

    :return: data: converted phase data
    :rtype: ndarray
    """
    return data * ifc.MM_PER_METRE * (wavelength / (4 * math.pi))


def convert_mm_to_radians(data, wavelength):
    """
    Function to translates phase in units of radians to units in millimetres.

    :param ndarray data: Interferogram phase data array
    :param float wavelength: Radar wavelength in metres

    :return: data: converted phase data
    :rtype: ndarray
    """
    return data / ifc.MM_PER_METRE * ((4 * math.pi) / wavelength)


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

    # pragma: no cover
    return np.median(x[~np.isnan(x)])


def _is_interferogram(hdr):
    """
    Convenience function to determine if file is interferogram
    """
    return (ifc.PYRATE_WAVELENGTH_METRES in hdr) and \
           (hdr[ifc.INPUT_TYPE] == InputTypes.IFG if ifc.INPUT_TYPE in hdr else True)


def _is_coherence(hdr):
    """
    Convenience function to determine if file is interferogram
    """
    return (ifc.PYRATE_WAVELENGTH_METRES in hdr) and \
           (hdr[ifc.INPUT_TYPE] == InputTypes.COH if ifc.INPUT_TYPE in hdr else False)


def _is_baseline(hdr):
    """
    Convenience function to determine if file is baseline file
    """
    return (ifc.PYRATE_WAVELENGTH_METRES in hdr) and \
           (hdr[ifc.INPUT_TYPE] == InputTypes.BASE if ifc.INPUT_TYPE in hdr else False)


def _is_lookuptable(hdr):
    """
    Convenience function to determine if file is lookup table file
    """
    return (ifc.PYRATE_WAVELENGTH_METRES in hdr) and \
           (hdr[ifc.INPUT_TYPE] == InputTypes.LT if ifc.INPUT_TYPE in hdr else False)


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
    bytes_per_col, fmtstr = data_format(ifg_proc, _is_interferogram(header), ncols)
    if _is_interferogram(header) and ifg_proc == ROIPAC:
        # roipac ifg has 2 bands
        _check_raw_data(bytes_per_col*2, data_path, ncols, nrows)
    else:
        _check_raw_data(bytes_per_col, data_path, ncols, nrows)

    # position and projection data
    gt = [
        header[ifc.PYRATE_LONG], header[ifc.PYRATE_X_STEP], 0,
        header[ifc.PYRATE_LAT], 0, header[ifc.PYRATE_Y_STEP]
    ]
    srs = osr.SpatialReference()
    res = srs.SetWellKnownGeogCS(header[ifc.PYRATE_DATUM])
    if res:
        msg = f'Unrecognised projection: {header[ifc.PYRATE_DATUM]}'
        raise GeotiffException(msg)

    wkt = srs.ExportToWkt()
    is_float = _is_interferogram(header)
    is_float |= _is_incidence(header)
    is_float |= _is_coherence(header)
    is_float |= _is_baseline(header)
    is_float |= _is_lookuptable(header)
    dtype = 'float32' if is_float else 'int16'

    # get subset of metadata relevant to PyRate
    md = collate_metadata(header)

    # create GDAL object
    ds = gdal_dataset(
        dest, ncols, nrows,
        driver="GTiff", bands=1, dtype=dtype,
        metadata=md, crs=wkt, geotransform=gt,
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
    md = {}

    def __common_ifg_coh_update(header, md):
        for k in [ifc.PYRATE_WAVELENGTH_METRES, ifc.PYRATE_TIME_SPAN,
                  ifc.PYRATE_INSAR_PROCESSOR,
                  ifc.FIRST_DATE, ifc.SECOND_DATE,
                  ifc.DATA_UNITS]:
            md.update({k: str(header[k])})
        if header[ifc.PYRATE_INSAR_PROCESSOR] == GAMMA:
            for k in [ifc.FIRST_TIME, ifc.SECOND_TIME,
                      ifc.PYRATE_NROWS, ifc.PYRATE_NCOLS,
                      ifc.PYRATE_INCIDENCE_DEGREES, ifc.PYRATE_HEADING_DEGREES,
                      ifc.PYRATE_AZIMUTH_DEGREES, ifc.PYRATE_RANGE_PIX_METRES,
                      ifc.PYRATE_RANGE_N, ifc.PYRATE_RANGE_LOOKS,
                      ifc.PYRATE_AZIMUTH_PIX_METRES, ifc.PYRATE_AZIMUTH_N,
                      ifc.PYRATE_AZIMUTH_LOOKS, ifc.PYRATE_PRF_HERTZ,
                      ifc.PYRATE_NEAR_RANGE_METRES, ifc.PYRATE_SAR_EARTH_METRES,
                      ifc.PYRATE_SEMI_MAJOR_AXIS_METRES, ifc.PYRATE_SEMI_MINOR_AXIS_METRES]:
                md.update({k: str(header[k])})
            if ifc.PYRATE_BASELINE_T in header:
                for k in [ifc.PYRATE_BASELINE_T, ifc.PYRATE_BASELINE_C,
                          ifc.PYRATE_BASELINE_N, ifc.PYRATE_BASELINE_RATE_T,
                          ifc.PYRATE_BASELINE_RATE_C, ifc.PYRATE_BASELINE_RATE_N]:
                    md.update({k: str(header[k])})

    if _is_coherence(header):
        __common_ifg_coh_update(header, md)
        md.update({ifc.DATA_TYPE: ifc.COH})
    elif _is_baseline(header):
        __common_ifg_coh_update(header, md)
        md.update({ifc.DATA_TYPE: ifc.BASE})
    elif _is_interferogram(header):
        __common_ifg_coh_update(header, md)
        md.update({ifc.DATA_TYPE: ifc.ORIG})
    elif _is_lookuptable(header):
        md.update({ifc.DATA_TYPE: ifc.LT})
    elif _is_incidence(header):
        md.update({ifc.DATA_TYPE: ifc.INCIDENCE})
    else:  # must be dem
        md.update({ifc.DATA_TYPE: ifc.DEM})

    return md


def data_format(ifg_proc, is_ifg, ncols):
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
        msg = f'Unrecognised InSAR Processor: {ifg_proc}'
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

    # set data type metadata
    ds.SetMetadataItem(ifc.DATA_TYPE, str(md[ifc.DATA_TYPE]))

    # set other metadata
    metadata = [
        ifc.SEQUENCE_POSITION, ifc.PYRATE_REFPIX_X, ifc.PYRATE_REFPIX_Y, ifc.PYRATE_REFPIX_LAT,
        ifc.PYRATE_REFPIX_LON, ifc.PYRATE_MEAN_REF_AREA, ifc.PYRATE_STDDEV_REF_AREA,
        ifc.EPOCH_DATE, C.LOS_PROJECTION.upper(), C.SIGNAL_POLARITY.upper(),
        C.VELERROR_NSIG.upper()
    ]

    for k in metadata:
        if k in md:
            ds.SetMetadataItem(k, str(md[k]))

    # write data to geotiff
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(data, 0, 0)

    del ds


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
        _, _, _ = data.shape
    elif data.ndim == 2:
        _, _ = data.shape
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
    cols = np.array_split(range(no_x), ncols)
    rows = np.array_split(range(no_y), nrows)
    return [
        Tile(i, (r[0], c[0]), (r[-1]+1, c[-1]+1)) for i, (r, c) in enumerate(product(rows, cols))
    ]


def get_tiles(ifg_path, rows, cols) -> List[Tile]:
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


def nan_and_mm_convert(ifg, params):
    """
    Perform millimetre and nan conversion on interferogram data

    :param Ifg instance ifg: Interferogram class instance
    :param dict params: Dictionary of parameters

    :return: None, data modified internally
    """
    nan_conversion = params[C.NAN_CONVERSION]
    if nan_conversion:  # nan conversion happens here in networkx mst
        # if not ifg.nan_converted:
        ifg.nodata_value = params[C.NO_DATA_VALUE]
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


class PrereadIfg:
    """
    Convenience class for handling pre-calculated ifg params
    """
    # pylint: disable=too-many-arguments
    # pylint: disable=too-many-instance-attributes
    def __init__(self, path, tmp_path, nan_fraction, first, second, time_span,
                 nrows, ncols, metadata):
        self.path = path
        self.tmp_path = tmp_path
        self.nan_fraction = nan_fraction
        self.first = first
        self.second = second
        self.time_span = time_span
        self.nrows = nrows
        self.ncols = ncols
        self.shape = (nrows, ncols)
        self.metadata = metadata


def save_numpy_phase(ifg_paths, params):
    """
    Split interferogram phase data in to tiles (if they exist in the params
    dict) and save as numpy array files on disk.

    :param list ifg_paths: List of strings for interferogram paths
    :param dict params: Dictionary of configuration parameters

    :return: None, numpy file saved to disk
    """
    tiles = params['tiles']
    outdir = params[C.TMPDIR]
    if not os.path.exists(outdir):
        mkdir_p(outdir)
    for ifg_path in mpiops.array_split(ifg_paths):
        ifg = Ifg(ifg_path)
        ifg.open()
        phase_data = ifg.phase_data
        bname = basename(ifg_path).split('.')[0]
        for t in tiles:
            p_data = phase_data[t.top_left_y:t.bottom_right_y,
                                t.top_left_x:t.bottom_right_x]
            phase_file = f'phase_data_{bname}_{t.index}.npy'
            np.save(file=join(outdir, phase_file),
                    arr=p_data)
        ifg.close()
    mpiops.comm.barrier()
    log.debug(f'Finished writing phase_data to numpy files in {outdir}')


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

    if not all(flags) and any(flags):
        log.debug('Detected mix of corrected and uncorrected interferograms')
        for flag in flags:
            if flag:
                msg = '{i.data_path}: correction detected'
            else:
                msg = '{i.data_path}: correction NOT detected'
            log.debug(msg)
            close_all(ifgs)
            raise CorrectionStatusError(msg)

    log.debug('Calculating corrections')
    close_all(ifgs)
    return False


class CorrectionStatusError(Exception):
    """
    Generic class for correction status errors.
    """


def extract_epochs_from_filename(filename_with_epochs: str) -> List[str]:
    """Extract the 8 or 6 digit epochs from IFG filenames"""
    src_epochs = re.findall(r"(\d{8})", str(filename_with_epochs))
    if not len(src_epochs) > 0:
        src_epochs = re.findall(r"(\d{6})", str(filename_with_epochs))
    return src_epochs


def mpi_vs_multiprocess_logging(step, params):
    """Logs the state of job processing (MPI vs. parallel)"""
    if mpiops.size > 1:  # Over-ride input options if this is an MPI job
        log.info(f"Running '{step}' step with MPI using {mpiops.size} processes")
        log.warning("Disabling joblib parallel processing (setting parallel = 0)")
        params[C.PARALLEL] = 0
    else:
        if params[C.PARALLEL] == 1:
            log.info(f"Running '{step}' step in parallel using {params[C.PROCESSES]} processes")
        else:
            log.info(f"Running '{step}' step in serial")


def dem_or_ifg(data_path: str) -> Union[Ifg, DEM]:
    """
    Returns an Ifg or DEM class object from input geotiff file.

    :param data_path: file path name

    :return: Interferogram or DEM object from input file
    :rtype: Ifg or DEM class object
    """
    ds = gdal.Open(data_path)
    md = ds.GetMetadata()
    is_ifg = ifc.FIRST_DATE in md
    ds = None  # close the dataset

    return Ifg(data_path) if is_ifg else DEM(data_path)


def join_dicts(dicts: List[dict]) -> dict:
    """
    Function to concatenate a list of dictionaries of distinct keys.
    """
    if dicts is None:  # pragma: no cover
        return {}
    assembled_dict = {k: v for D in dicts for k, v in D.items()}
    return assembled_dict


def iterable_split(func: Callable, iterable: Iterable, params: dict, *args, **kwargs) -> np.ndarray:
    """
    # TODO: a faster version using buffer-provider objects via the uppercase communication method
    A faster version of iterable/tiles_split is possible when the return values from each process
    is of the same size and will be addressed in future. In this case a buffer-provider object can
    be sent between processes using the uppercase communication (like Gather instead of gather)
    methods which can be significantly faster.
    """
    if params[C.PARALLEL]:
        ret_combined = {}
        rets = Parallel(
            n_jobs=params[C.PROCESSES],
            verbose=joblib_log_level(C.LOG_LEVEL)
        )(delayed(func)(t, params, *args, **kwargs) for t in iterable)
        for i, r in enumerate(rets):
            ret_combined[i] = r
    else:
        iterable_with_index = list(enumerate(iterable))
        process_iterables = mpiops.array_split(iterable_with_index)
        ret_combined = {}
        for i, t in process_iterables:
            ret_combined[i] = func(t, params, *args, **kwargs)
        ret_combined = join_dicts(mpiops.comm.allgather(ret_combined))
    ret = np.array([v[1] for v in ret_combined.items()], dtype=object)
    mpiops.comm.barrier()
    return ret


def tiles_split(func: Callable, params: dict, *args, **kwargs) -> np.ndarray:
    """
    Function to pass tiles of a full array to an array processing function call.
    :param func: Name of function to pass tiles to.
    :param params: Dictionary of PyRate configuration parameters; must contain a 'tiles' list
    """
    tiles = params[C.TILES]
    return iterable_split(func, tiles, params, *args, **kwargs)


def output_tiff_filename(inpath: str, outpath: str) -> str:
    """
    Output geotiff filename for a given input filename.

    :param inpath: Path of input file location.
    :param outpath: Path of output file location.
    :return: name: Geotiff filename for the given file.
    """
    fname, ext = os.path.basename(inpath).split('.')
    outpath = os.path.dirname(inpath) if outpath is None else outpath
    if ext == 'tif':
        name = os.path.join(outpath, fname + '.tif')
    else:
        name = os.path.join(outpath, fname + '_' + ext + '.tif')
    return name


def remove_file_if_exists(filename: str) -> None:
    """
    Function to remove a file if it already exists.
    :param filename: Name of file to be removed.
    """
    try:
        os.remove(filename)
    except OSError:
        pass
