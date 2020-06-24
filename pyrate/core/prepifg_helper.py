#   This Python module is part of the PyRate software package
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
"""
This Python module converts interferogram input files to a common geotiff
format with PyRate specific metadata headers. The module also implements
multilooking/downsampling and cropping operations to reduce the size of
the computational problem.
"""
# pylint: disable=too-many-arguments,invalid-name
from collections import namedtuple
from math import modf
from numbers import Number
from decimal import Decimal
from typing import List, Tuple, Union
from numpy import array, nan, isnan, nanmean, float32, zeros, sum as nsum
from osgeo import gdal

from pyrate.core.gdal_python import crop_resample_average
from pyrate.core import config as cf
from pyrate.core.shared import output_tiff_filename, dem_or_ifg, Ifg, DEM

CustomExts = namedtuple('CustExtents', ['xfirst', 'yfirst', 'xlast', 'ylast'])


# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
ALREADY_SAME_SIZE = 4
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE]

GRID_TOL = 1e-6


def get_analysis_extent(crop_opt: int, rasters: List[Union[Ifg, DEM]], xlooks: int, ylooks: int,
                        user_exts: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    """
    Function checks prepifg parameters and returns extents/bounding box.

    :param int crop_opt: Cropping option
    :param list rasters: List of either Ifg or DEM class objects
    :param int xlooks: Number of multi-looks in x
    :param int ylooks: Number of multi-looks in y
    :param tuple user_exts: Tuple of user defined cropping coordinates

    :return: extents: tuple of four bounding coordinates
    :rtype: tuple
    """

    if crop_opt not in CROP_OPTIONS:
        raise PreprocessError("Unrecognised crop option: %s" % crop_opt)

    if crop_opt == CUSTOM_CROP:
        if not user_exts:
            raise PreprocessError('No custom cropping extents specified')
        elif len(user_exts) != 4:  # check for required numbers
            raise PreprocessError('Custom extents must have all 4 values')
        elif len(user_exts) == 4:  # check for non floats
            if not all([_is_number(z) for z in user_exts]):
                raise PreprocessError('Custom extents must be 4 numbers')

    for raster in rasters:
        if not raster.is_open:
            raster.open()

    _check_looks(xlooks, ylooks)
    _check_resolution(rasters)

    return _get_extents(rasters, crop_opt, user_exts)


def _is_number(s):
    """
    Check whether string can be converted to float
    """
    try:
        float(s)
        return True
    except TypeError:  # for example if 'None' is sent
        return False
    except ValueError:  # for example if a 'string' is sent
        return False


def _check_looks(xlooks, ylooks):
    """
    Convenience function to verify that looks parameters are valid.
    """

    if not (isinstance(xlooks, Number) and
            isinstance(ylooks, Number)):  # pragma: no cover
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise PreprocessError(msg)

    if not (xlooks > 0 and ylooks > 0):  # pragma: no cover
        msg = "Invalid looks parameter(s), x: %s, y: %s. " \
              "Looks must be an integer greater than zero" % (xlooks, ylooks)
        raise PreprocessError(msg)


def _check_resolution(ifgs):
    """
    Convenience function to verify Ifg resolutions are equal.
    """

    for var in ['x_step', 'y_step']:
        values = array([getattr(i, var) for i in ifgs])
        if not (values == values[0]).all():  # pragma: no cover
            msg = "Grid resolution does not match for %s" % var
            raise PreprocessError(msg)


def _get_extents(ifgs, crop_opt, user_exts=None):
    """
    Convenience function that returns extents/bounding box.
    MINIMUM_CROP = 1
    MAXIMUM_CROP = 2
    CUSTOM_CROP = 3
    ALREADY_SAME_SIZE = 4
    """
    if crop_opt == MINIMUM_CROP:
        extents = _min_bounds(ifgs)
    elif crop_opt == MAXIMUM_CROP:
        extents = _max_bounds(ifgs)
    elif crop_opt == CUSTOM_CROP:
        extents = _custom_bounds(ifgs, *user_exts)
        # only need to check crop coords when custom bounds are supplied
        _check_crop_coords(ifgs, *extents)
    else:
        extents = _get_same_bounds(ifgs)

    return extents


def prepare_ifg(raster_path, xlooks, ylooks, exts, thresh, crop_opt, header, write_to_disk=True, out_path=None,
                coherence_path=None, coherence_thresh=None):
    """
    Open, resample, crop and optionally save to disk an interferogram or DEM.
    Returns are only given if write_to_disk=False

    :param str raster_path: Input raster file path name
    :param int xlooks: Number of multi-looks in x; 5 is 5 times smaller,
        1 is no change
    :param int ylooks: Number of multi-looks in y
    :param tuple exts: Tuple of user defined georeferenced extents for
        new file: (xfirst, yfirst, xlast, ylast)cropping coordinates
    :param float thresh: see thresh in prepare_ifgs()
    :param int crop_opt: Crop option
    :param bool write_to_disk: Write new data to disk
    :param str out_path: Path for output file
    :param dict header: dictionary of metadata from header file

    :return: resampled_data: output cropped and resampled image
    :rtype: ndarray
    :return: out_ds: destination gdal dataset object
    :rtype: gdal.Dataset
    """
    do_multilook = xlooks > 1 or ylooks > 1
    # resolution=None completes faster for non-multilooked layers in gdalwarp
    resolution = [None, None]
    raster = dem_or_ifg(raster_path)
    if not raster.is_open:
        raster.open()
    if do_multilook:
        resolution = [xlooks * raster.x_step, ylooks * raster.y_step]

    # cut, average, resample the final output layers
    op = output_tiff_filename(raster.data_path, out_path)
    looks_path = cf.mlooked_path(op, ylooks, crop_opt)

    if xlooks != ylooks:
        raise ValueError('X and Y looks mismatch')

    #     # Add missing/updated metadata to resampled ifg/DEM
    #     new_lyr = type(ifg)(looks_path)
    #     new_lyr.open(readonly=True)
    #     # for non-DEMs, phase bands need extra metadata & conversions
    #     if hasattr(new_lyr, "phase_band"):
    #         # TODO: LOS conversion to vertical/horizontal (projection)
    #         # TODO: push out to workflow
    #         #if params.has_key(REPROJECTION_FLAG):
    #         #    reproject()
    driver_type = 'GTiff' if write_to_disk else 'MEM'
    resampled_data, out_ds = crop_resample_average(
        input_tif=raster.data_path, extents=exts, new_res=resolution, output_file=looks_path, thresh=thresh,
        out_driver_type=driver_type, hdr=header, coherence_path=coherence_path, coherence_thresh=coherence_thresh
    )

    return resampled_data, out_ds


# TODO: crop options 0 = no cropping? get rid of same size
def prepare_ifgs(raster_data_paths, crop_opt, xlooks, ylooks, headers, thresh=0.5, user_exts=None, write_to_disc=True,
                 out_path=None):
    """
    Wrapper function to prepare a sequence of interferogram files for
    PyRate analysis. See prepifg.prepare_ifg() for full description of
    inputs and returns.
    
    Note: function need refining for crop options

    :param list raster_data_paths: List of interferogram file paths
    :param int crop_opt: Crop option
    :param int xlooks: Number of multi-looks in x; 5 is 5 times smaller, 1 is no change
    :param int ylooks: Number of multi-looks in y
    :param float thresh: see thresh in prepare_ifgs()
    :param tuple user_exts: Tuple of user defined georeferenced extents for
        new file: (xfirst, yfirst, xlast, ylast)cropping coordinates
    :param bool write_to_disk: Write new data to disk

    :return: resampled_data: output cropped and resampled image
    :rtype: ndarray
    :return: out_ds: destination gdal dataset object
    :rtype: List[gdal.Dataset]
    """
    # use metadata check to check whether it's a dem or ifg
    rasters = [dem_or_ifg(r) for r in raster_data_paths]
    exts = get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts)
    return [prepare_ifg(d, xlooks, ylooks, exts, thresh, crop_opt, h, write_to_disc, out_path) for d, h
            in zip(raster_data_paths, headers)]


# TODO: Not being used. Remove in future?
def _resample(data, xscale, yscale, thresh):
    """
    Resamples/averages 'data' to return an array from the averaging of blocks
    of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

    :param data: source array to resample to different size
    :param xscale: number of cells to average along X axis
    :param yscale: number of Y axis cells to average
    :param thresh: minimum allowable
        proportion of NaN cells (range from 0.0-1.0), eg. 0.25 = 1/4 or
        more as NaNs results in a NaN value for the output cell.
    """
    # TODO: make more efficient
    if thresh < 0 or thresh > 1:
        raise ValueError("threshold must be >= 0 and <= 1")

    xscale = int(xscale)
    yscale = int(yscale)
    ysize, xsize = data.shape
    xres, yres = int(xsize / xscale), int(ysize / yscale)
    dest = zeros((yres, xres), dtype=float32) * nan
    tile_cell_count = xscale * yscale

    # calc mean without nans (fractional threshold ignores tiles
    # with excess NaNs)
    for x in range(xres):
        for y in range(yres):
            tile = data[y * yscale: (y+1) * yscale, x * xscale: (x+1) * xscale]
            nan_fraction = nsum(isnan(tile)) / float(tile_cell_count)
            if nan_fraction < thresh or (nan_fraction == 0 and thresh == 0):
                dest[y, x] = nanmean(tile)
    return dest


def _min_bounds(ifgs: List[Ifg]) -> Tuple[float, float, float, float]:
    """
    Returns bounds for overlapping area of the given interferograms.
    """

    xmin = max([i.x_first for i in ifgs])
    ymax = min([i.y_first for i in ifgs])
    xmax = min([i.x_last for i in ifgs])
    ymin = max([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def _max_bounds(ifgs: List[Ifg]) -> Tuple[float, float, float, float]:
    """
    Returns bounds for the total area covered by the given interferograms.
    """

    xmin = min([i.x_first for i in ifgs])
    ymax = max([i.y_first for i in ifgs])
    xmax = max([i.x_last for i in ifgs])
    ymin = min([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def _get_same_bounds(ifgs: List[Ifg]) -> Tuple[float, float, float, float]:
    """
    Check and return bounding box for ALREADY_SAME_SIZE option.
    """

    tfs = [i.dataset.GetGeoTransform() for i in ifgs]

    equal = []

    for t in tfs[1:]:
        for i,tf in enumerate(tfs[0]):

            if round(Decimal (tf),4) == round(Decimal (t[i]),4):
                equal.append(True)
            else:
                equal.append(False)

    if not all(equal):
        msg = 'Ifgs do not have the same bounding box for crop option: %s'
        raise PreprocessError(msg % ALREADY_SAME_SIZE)
    ifg = ifgs[0]
    xmin, xmax = ifg.x_first, ifg.x_last
    ymin, ymax = ifg.y_first, ifg.y_last

    # swap y_first & y_last when using southern hemisphere -ve coords
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    return xmin, ymin, xmax, ymax


def _custom_bounds(ifgs, xw, ytop, xe, ybot):
    """
    Check and modify input custom crop bounds to line up with grid interval
    """
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-branches
    msg = 'Cropped image bounds are outside the original image bounds'
    i = ifgs[0]

    if ytop < ybot:
        raise PreprocessError('ERROR Custom crop bounds: '
                              'ifgyfirst must be greater than ifgylast')

    if xe < xw:
        raise PreprocessError('ERROR Custom crop bounds: '
                              'ifgxfirst must be greater than ifgxlast')

    for par, crop, orig, step in zip(['x_first', 'x_last', 'y_first', 'y_last'],
                                     [xw, xe, ytop, ybot],
                                     [i.x_first, i.x_last, i.y_first, i.y_last],
                                     [i.x_step, i.x_step, i.y_step, i.y_step]):
        diff = crop - orig
        nint = round(diff / step)

        if par == 'x_first':
            if diff < 0:
                raise PreprocessError(msg)
            xmin = orig + (nint * step)

        elif par == 'x_last':
            if diff > 0:
                raise PreprocessError(msg)
            xmax = orig + (nint * step)

        elif par == 'y_first':
            if diff > 0:
                raise PreprocessError(msg)
            y1 = orig + (nint * step)

        elif par == 'y_last':
            if diff < 0:
                raise PreprocessError(msg)
            y2 = orig + (nint * step)
        else:
            raise ValueError('Value error in supplied custom bounds')

    if y2 > y1:
        ymin = y1
        ymax = y2
    else:
        ymin = y2
        ymax = y1

    return xmin, ymin, xmax, ymax


def _check_crop_coords(ifgs, xmin, ymin, xmax, ymax):
    """
    Ensures cropping coords line up with grid system within tolerance.
    """

    # NB: assumption is the first Ifg is correct, so only test against it
    i = ifgs[0]

    for par, crop, step in zip(['x_first', 'x_last', 'y_first', 'y_last'],
                               [xmin, xmax, ymax, ymin],
                               [i.x_step, i.x_step, i.y_step, i.y_step]):

        # is diff of the given extent from grid a multiple of X|Y_STEP ?
        param = getattr(i, par)
        diff = abs(crop - param)
        remainder = abs(modf(diff / step)[0])

        # handle cases where division gives remainder near zero, or just < 1
        if (remainder > GRID_TOL) and (remainder < (1 - GRID_TOL)):  # pragma: no cover
            msg = "%s crop extent not within %s of grid coordinate"
            raise PreprocessError(msg % (par, GRID_TOL))


class PreprocessError(Exception):
    """
    Preprocess exception
    """
