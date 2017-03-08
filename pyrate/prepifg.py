#   This Python module is part of the PyRate software package
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
This Python module converts interferogram input files to a common geotiff
format with PyRate specific metadata headers. The module also implements
multilooking/downsampling and cropping operations to reduce the size of
the computational problem.
"""
# pylint: disable=too-many-arguments,invalid-name
import os
import shutil
from collections import namedtuple
from math import modf
from numbers import Number
from subprocess import check_call
from tempfile import mkstemp

from numpy import array, where, nan, isnan, nanmean, float32, zeros, \
    sum as nsum
from osgeo import gdal

from pyrate import config as cf
from pyrate import gdal_python as gdalwarp
from pyrate import ifgconstants as ifc
from pyrate.shared import Ifg, DEM

CustomExts = namedtuple('CustExtents', ['xfirst', 'yfirst', 'xlast', 'ylast'])


# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
ALREADY_SAME_SIZE = 4
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE]

GRID_TOL = 1e-6


def is_number(s):
    """ check whether string can be converted to float"""
    try:
        float(s)
        return True
    except TypeError:  # for example if 'None' is sent
        return False
    except ValueError:  # for example if a 'string' is sent
        return False


def get_analysis_extent(
        crop_opt,
        rasters,
        xlooks,
        ylooks,
        user_exts):
    """
    Returns extents/bounding box args for gdalwarp as strings.
    """

    if crop_opt not in CROP_OPTIONS:
        raise PreprocessError("Unrecognised crop option: %s" % crop_opt)

    if crop_opt == CUSTOM_CROP:
        if not user_exts:
            raise PreprocessError('No custom cropping extents specified')
        elif len(user_exts) != 4:  # check for required numbers
            raise PreprocessError('Custom extents must have all 4 values')
        elif len(user_exts) == 4:  # check for non floats
            if not all([is_number(z) for z in user_exts]):
                raise PreprocessError('Custom extents must be 4 numbers')

    for raster in rasters:
        if not raster.is_open:
            raster.open()

    check_looks(xlooks, ylooks)
    check_resolution(rasters)

    return get_extents(rasters, crop_opt, user_exts)


def prepare_ifg(
        raster_path,
        xlooks,
        ylooks,
        exts,
        thresh,
        crop_opt,
        write_to_disc=True):
    """
    resample and crop ifgs/dems
    """

    do_multilook = xlooks > 1 or ylooks > 1
    # resolution=None completes faster for non-multilooked layers in gdalwarp
    resolution = [None, None]
    raster = dem_or_ifg(raster_path)
    if not raster.is_open:
        raster.open()
    if do_multilook:
        resolution = [xlooks * raster.x_step, ylooks * raster.y_step]

    if not do_multilook and crop_opt == ALREADY_SAME_SIZE:
        renamed_path = \
            cf.mlooked_path(raster.data_path, looks=xlooks, crop_out=crop_opt)
        shutil.copy(raster.data_path, renamed_path)
        # set metadata to indicated has been cropped and multilooked
        # copy file with mlooked path
        return dummy_warp(renamed_path)

    return warp(raster, xlooks, ylooks, exts, resolution, thresh,
                crop_opt, write_to_disc)


def dummy_warp(renamed_path):
    """ convenience dummy operation """
    ifg = dem_or_ifg(renamed_path)
    ifg.open()
    ifg.dataset.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)
    data = ifg.dataset.ReadAsArray()
    return data, ifg.dataset


# TODO: crop options 0 = no cropping? get rid of same size
# (but it is in explained file)
def prepare_ifgs(
        raster_data_paths,
        crop_opt,
        xlooks,
        ylooks,
        thresh=0.5,
        user_exts=None,
        write_to_disc=True):
    """
    Produces multilooked/resampled data files for PyRate analysis.

    :param raster_data_paths: sequence of Ifg data paths
            (Currently DEMs are not supported)
    :param crop_opt: integer cropping type option (see config)
    :param xlooks: multilooking factor for the X axis
    :param ylooks: Y axis multilooking factor
    :param float thresh: (0.0, 1.0). Controls NaN handling when resampling to
     coarser grids. Value is the proportion above which the number of NaNs in
     an area is considered invalid. thresh=0 resamples to NaN if 1 or more
     contributing cells are NaNs. At 0.25, it resamples to NaN if 1/4 or
     more contributing cells are NaNs. At 1.0, areas are resampled to NaN
     only if all contributing cells are NaNs.
    :param user_exts: CustomExts tuple with user sepcified lat long corners
    :param verbose: Controls level of gdalwarp output
    :param write_to_disc: bool, whether to write to disc during warp
    """
    # use metadata check to check whether it's a dem or ifg
    rasters = [dem_or_ifg(r) for r in raster_data_paths]
    exts = get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts)

    return [prepare_ifg(d, xlooks, ylooks, exts, thresh, crop_opt,
                        write_to_disc)
            for d in raster_data_paths]


def dem_or_ifg(data_path):
    """ whether tif is a dem or an ifg """
    ds = gdal.Open(data_path)
    md = ds.GetMetadata()
    if 'DATE' in md:  # ifg
        return Ifg(data_path)
    else:
        return DEM(data_path)


def get_extents(ifgs, crop_opt, user_exts=None):
    """
    Returns extents/bounding box args for gdalwarp as strings.
    """

    if crop_opt == MINIMUM_CROP:
        extents = min_bounds(ifgs)
    elif crop_opt == MAXIMUM_CROP:
        extents = max_bounds(ifgs)
    elif crop_opt == CUSTOM_CROP:
        extents = custom_bounds(ifgs, *user_exts)
    else:
        extents = get_same_bounds(ifgs)

    check_crop_coords(ifgs, *extents)
    return extents


# TODO: push out to a shared module
def _file_ext(raster):
    """
    Returns file ext string based on type of raster.
    """
    if isinstance(raster, Ifg):
        return "tif"
    elif isinstance(raster, DEM):
        return "dem"
    else:
        # TODO: several possible file types to implement:
        # Coherence file: single band
        # LOS file:  has 2 bands: beam incidence angle & ground azimuth)
        # Baseline file: perpendicular baselines (single band?)
        raise NotImplementedError("Missing raster types for LOS, Coherence and baseline")


def _resample_ifg(ifg, cmd, x_looks, y_looks, thresh, md=None):
    """
    Convenience function to resample data from a given Ifg (more coarse).
    """

    # HACK: create tmp ifg, extract data array for manual resampling as gdalwarp
    # lacks Pirate's averaging method
    fp, tmp_path = mkstemp(suffix='.tif')
    check_call(cmd + [ifg.data_path, tmp_path])

    # now write the metadata from the input to the output
    if md is not None:
        new_lyr = gdal.Open(tmp_path)
        for k, v in md.iteritems():
            new_lyr.SetMetadataItem(k, v)
        new_lyr = None  # manually close

    tmp = type(ifg)(tmp_path)  # dynamically handle Ifgs & Rasters
    tmp.open()

    if isinstance(ifg, Ifg):
        # TODO: add an option to retain amplitude band (resample this if reqd)
        data = tmp.phase_band.ReadAsArray()
        data = where(data == 0, nan, data)  # flag incoherent cells as NaNs
    elif isinstance(ifg, DEM):
        data = tmp.height_band.ReadAsArray()
    else:
        # TODO: need to handle resampling of LOS and baseline files
        raise NotImplementedError("Resampling LOS & baseline not implemented")

    tmp.close()  # manual close
    os.close(fp)
    os.remove(tmp_path)
    return resample(data, x_looks, y_looks, thresh)


def warp(ifg, x_looks, y_looks, extents, resolution, thresh, crop_out,
         write_to_disc=True):
    """
    Resamples 'ifg' and returns a new Ifg obj.

    :param xlooks: int
        factor to scale X axis by, 5 is 5x smaller, 1 is no change.
    :param ylooks: int
        as xlooks, but for Y axis
    :param extents: tuple
        georeferenced extents for new file: (xfirst, yfirst, xlast, ylast)
    :param resolution: [xres, yres] or None.
        Sets resolution output Ifg metadata.
        Use *None* if raster size is not being changed.
    :param thresh: see thresh in prepare_ifgs().
    :param verbose: True to print gdalwarp output to stdout
    :param write_to_disc: bool, whether to write to disc during warp
    """
    if x_looks != y_looks:
        raise ValueError('X and Y looks mismatch')

    # cut, average, resample the final output layers
    looks_path = cf.mlooked_path(ifg.data_path, y_looks, crop_out)

    #     # Add missing/updated metadata to resampled ifg/DEM
    #     new_lyr = type(ifg)(looks_path)
    #     new_lyr.open(readonly=True)
    #     # for non-DEMs, phase bands need extra metadata & conversions
    #     if hasattr(new_lyr, "phase_band"):
    #         # TODO: LOS conversion to vertical/horizontal (projection)
    #         # TODO: push out to workflow
    #         #if params.has_key(REPROJECTION_FLAG):
    #         #    reproject()
    driver_type = 'GTiff' if write_to_disc else 'MEM'
    resampled_data, out_ds = gdalwarp.crop_resample_average(
        input_tif=ifg.data_path,
        extents=extents,
        new_res=resolution,
        output_file=looks_path,
        thresh=thresh,
        out_driver_type=driver_type)

    if not write_to_disc:
        return resampled_data, out_ds


def resample(data, xscale, yscale, thresh):
    """
    # TODO: make more efficient
    This is the slowest step in prepifg
    Resamples/averages 'data' to return an array from the averaging of blocks
    of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

    :param data: source array to resample to different size
    :param xscale: number of cells to average along X axis
    :param yscale: number of Y axis cells to average
    :param thresh: minimum allowable
        proportion of NaN cells (range from 0.0-1.0), eg. 0.25 = 1/4 or
        more as NaNs results in a NaN value for the output cell.
    """
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


def check_resolution(ifgs):
    """
    Verifies Ifg resolutions are equal for the given grids.
    """

    for var in ['x_step', 'y_step']:
        values = array([getattr(i, var) for i in ifgs])
        if not (values == values[0]).all():  # pragma: no cover
            msg = "Grid resolution does not match for %s" % var
            raise PreprocessError(msg)


def check_looks(xlooks, ylooks):
    """
    Verifies looks parameters are valid.
    """

    if not (isinstance(xlooks, Number) and
            isinstance(ylooks, Number)):  # pragma: no cover
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise PreprocessError(msg)

    if not (xlooks > 0 and ylooks > 0):  # pragma: no cover
        msg = "Invalid looks parameter(s), x: %s, y: %s. " \
              "Looks must be an integer greater than zero" % (xlooks, ylooks)
        raise PreprocessError(msg)


def min_bounds(ifgs):
    """
    Returns bounds for overlapping area of the given interferograms.
    """

    xmin = max([i.x_first for i in ifgs])
    ymax = min([i.y_first for i in ifgs])
    xmax = min([i.x_last for i in ifgs])
    ymin = max([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def max_bounds(ifgs):
    """
    Returns bounds for the total area covered by the given interferograms.
    """

    xmin = min([i.x_first for i in ifgs])
    ymax = max([i.y_first for i in ifgs])
    xmax = max([i.x_last for i in ifgs])
    ymin = min([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def get_same_bounds(ifgs):
    """
    Check and return bounding box for ALREADY_SAME_SIZE option.
    """

    tfs = [i.dataset.GetGeoTransform() for i in ifgs]
    equal = [t == tfs[0] for t in tfs[1:]]
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


def custom_bounds(ifgs, xw, ytop, xe, ybot):
    """
    Check and modify input custom crop bounds to line up with grid interval
    """
    # pylint: disable=too-many-locals
    # pylint: disable=too-many-branches
    msg = 'Cropped image bounds exceed original image bounds'
    i = ifgs[0]

    if ytop < ybot:
        raise PreprocessError('ERROR Custom crop bounds: '
                              'ifgyfirst must be greater than ifgylast')

    if xe < xw:
        raise PreprocessError('ERROR Custom crop bounds: '
                              'ifgxfirst must be greater than ifgxlast')

    for par, crop, orig, step in zip(['x_first', 'x_last',
                                      'y_first', 'y_last'],
                                     [xw, xe, ytop, ybot],
                                     [i.x_first, i.x_last,
                                      i.y_first, i.y_last],
                                     [i.x_step, i.x_step,
                                      i.y_step, i.y_step]):
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


def check_crop_coords(ifgs, xmin, ymin, xmax, ymax):
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
        if (remainder > GRID_TOL) and \
                (remainder < (1 - GRID_TOL)):  # pragma: no cover
            msg = "%s crop extent not within %s of grid coordinate"
            raise PreprocessError(msg % (par, GRID_TOL))


class PreprocessError(Exception):
    """
    Preprocess exception
    """


def extents_from_params(params):
    """ Custom extents from supplied parameters """
    keys = (cf.IFG_XFIRST, cf.IFG_YFIRST, cf.IFG_XLAST, cf.IFG_YLAST)
    return CustomExts(*[params[k] for k in keys])
