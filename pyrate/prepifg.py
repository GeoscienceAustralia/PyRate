"""
Prepares input files and associated data for the PyRate work flow.

Input rasters often may cropping, scaling, and multilooking/downsampling to
coarser grids before being processed. This module uses gdalwarp to handle these
operations.

The rasters need to be in GeoTIFF format with PyRate specific metadata headers.

Created on 23/10/2012

.. codeauthor:: Ben Davies
"""

# TODO: check new average option for gdalwarp (GDAL 1.10.x +)

import os
from math import modf
from numbers import Number
from tempfile import mkstemp
from itertools import product
from subprocess import check_call
from collections import namedtuple
from os.path import splitext
import shutil

import numpy as np
from numpy import array, where, nan, isnan, nanmean, float32, zeros, sum as nsum

from osgeo import gdal, osr, gdalconst
import pyrate.config as cfg
from pyrate.shared import Ifg, DEM, PHASE_BAND
from pyrate import gdal_python as gdalwarp

CustomExts = namedtuple('CustExtents', ['xfirst', 'yfirst', 'xlast', 'ylast'])


# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
ALREADY_SAME_SIZE = 4
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE]

GRID_TOL = 1e-6

def getAnalysisExtent(
        cropOpt,
        rasters,
        xlooks,
        ylooks,
        userExts):

    if cropOpt not in CROP_OPTIONS:
        raise PreprocessError("Unrecognised crop option: %s" % cropOpt)

    if cropOpt == CUSTOM_CROP and not userExts:
        raise PreprocessError('No custom cropping extents specified')

    for raster in rasters:
        if not raster.is_open:
            raster.open()

    check_looks(xlooks, ylooks)
    check_resolution(rasters)

    return get_extents(rasters, cropOpt, userExts)


def prepare_ifg(
        raster_path,
        xlooks,
        ylooks,
        exts,
        thresh,
        crop_opt,
        write_to_disc=True):
    # Determine cmd line args for gdalwarp calls for each ifg (gdalwarp has no
    # API. For resampling, gdalwarp is called 2x. 1st to subset the source data
    # for Pirate style averaging/resampling, 2nd to generate the final dataset
    # with correct extents/shape/cell count. Without resampling, gdalwarp is
    # only needed to cut out the required segment.

    do_multilook = xlooks > 1 or ylooks > 1
    # resolution=None completes faster for non-multilooked layers in gdalwarp
    resolution = [None, None]
    raster = Ifg(raster_path)
    if do_multilook:
        if not raster.is_open:
            raster.open()
        resolution = [xlooks * raster.x_step, ylooks * raster.y_step]

    if not do_multilook and crop_opt == ALREADY_SAME_SIZE:
        renamed_path = \
            mlooked_path(raster.data_path, looks=xlooks, crop_out=crop_opt)
        # copy file with mlooked path
        shutil.copy(raster.data_path, renamed_path)
        return None

    if not raster.is_open:
        raster.open()

    return warp(raster, xlooks, ylooks, exts, resolution, thresh,
                crop_opt, write_to_disc)


# TODO: crop options 0 = no cropping? get rid of same size (but it is in explained file)
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
    :param float thresh: (0.0, 1.0). Controls NaN handling when resampling to coarser grids.
         Value is the proportion above which the number of NaNs in an area is
         considered invalid. thresh=0 resamples to NaN if 1 or more contributing
         cells are NaNs. At 0.25, it resamples to NaN if 1/4 or more contributing
         cells are NaNs. At 1.0, areas are resampled to NaN only if all
         contributing cells are NaNs.
    :param user_exts: CustomExts tuple with user sepcified lat long corners
    :param verbose: Controls level of gdalwarp output
    """
    # TODO: make dems work in prep_ifgs again
    rasters = [Ifg(r) for r in raster_data_paths]
    exts = getAnalysisExtent(crop_opt, rasters, xlooks, ylooks, user_exts)

    return [prepare_ifg(d, xlooks, ylooks, exts, thresh, crop_opt, write_to_disc)
            for d in raster_data_paths]


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


def mlooked_path(path, looks, crop_out):
    """
    Adds suffix to path, for creating a new path for mlooked files.
    """
    base, ext = splitext(path)
    return "{base}_{looks}rlks_{crop_out}cr{ext}".format(
        base=base, looks=looks, crop_out=crop_out, ext=ext)


def warp_old(ifg, x_looks, y_looks, extents, resolution, thresh, crop_out, verbose,
         ret_ifg=True):
    """
    Resamples 'ifg' and returns a new Ifg obj.

    :param xlooks: integer factor to scale X axis by, 5 is 5x smaller, 1 is no change.
    :param ylooks: as xlooks, but for Y axis
    :param extents: georeferenced extents for new file: (xfirst, yfirst, xlast, ylast)
    :param resolution: [xres, yres] or None. Sets resolution output Ifg metadata.
         Use *None* if raster size is not being changed.
    :param thresh: see thresh in prepare_ifgs().
    :param verbose: True to print gdalwarp output to stdout
    """
    if x_looks != y_looks:
        raise ValueError('X and Y looks mismatch')

    # dynamically build command for call to gdalwarp
    cmd = ["gdalwarp", "-overwrite", "-srcnodata", "None", "-te"] + extents
    if not verbose:
        cmd.append("-q")

    # It appears that before vrions 1.10 gdal-warp did not copy meta-data
    # ... so we need to. Get the keys from the input here
    if sum((int(v)*i for v, i in zip(gdal.__version__.split('.')[:2], [10, 1]))) < 20:
        fl = ifg.data_path
        dat = gdal.Open(fl)
        md = {k:v for k, v in dat.GetMetadata().iteritems()}
    else:
        md = None

    # HACK: if resampling, cut segment with gdalwarp & manually average tiles
    data = None
    if resolution[0]:
        data = _resample_ifg(ifg, cmd, x_looks, y_looks, thresh, md)
        cmd += ["-tr"] + [str(r) for r in resolution] # change res of final output

    # use GDAL to cut (and resample) the final output layers
    looks_path = mlooked_path(ifg.data_path, y_looks, crop_out)
    cmd += [ifg.data_path, looks_path]

    check_call(cmd)
    # now write the metadata from the input to the output
    if md is not None:
        new_lyr = gdal.Open(looks_path)
        for k, v in md.iteritems():
            new_lyr.SetMetadataItem(k, v)

    # Add missing/updated metadata to resampled ifg/DEM
    new_lyr = type(ifg)(looks_path)
    new_lyr.open(readonly=False)
    # for non-DEMs, phase bands need extra metadata & conversions
    if hasattr(new_lyr, "phase_band"):
        if data is None:  # data wasn't resampled, so flag incoherent cells
            data = new_lyr.phase_band.ReadAsArray()
            data = np.where(np.isclose(data, 0.0, atol=1e-6), np.nan, data)

        # TODO: LOS conversion to vertical/horizontal (projection)
        # TODO: push out to workflow
        #if params.has_key(REPROJECTION_FLAG):
        #    reproject()

        # tricky: write either resampled or the basic cropped data to new layer
        new_lyr.phase_band.SetNoDataValue(nan)
        new_lyr.phase_band.WriteArray(data)
        new_lyr.nan_converted = True

    if ret_ifg:
        return data, looks_path
    else:
        return


def warp(ifg, x_looks, y_looks, extents, resolution, thresh, crop_out,
         write_to_disc=True):
    """
    Resamples 'ifg' and returns a new Ifg obj.

    :param xlooks: integer factor to scale X axis by, 5 is 5x smaller, 1 is no change.
    :param ylooks: as xlooks, but for Y axis
    :param extents: georeferenced extents for new file: (xfirst, yfirst, xlast, ylast)
    :param resolution: [xres, yres] or None. Sets resolution output Ifg metadata.
         Use *None* if raster size is not being changed.
    :param thresh: see thresh in prepare_ifgs().
    :param verbose: True to print gdalwarp output to stdout
    """
    if x_looks != y_looks:
        raise ValueError('X and Y looks mismatch')

    # cut, average, resample the final output layers
    looks_path = mlooked_path(ifg.data_path, y_looks, crop_out)

    #     # Add missing/updated metadata to resampled ifg/DEM
    #     new_lyr = type(ifg)(looks_path)
    #     new_lyr.open(readonly=True)
    #     # for non-DEMs, phase bands need extra metadata & conversions
    #     if hasattr(new_lyr, "phase_band"):
    #         # TODO: LOS conversion to vertical/horizontal (projection)
    #         # TODO: push out to workflow
    #         #if params.has_key(REPROJECTION_FLAG):
    #         #    reproject()
    if not write_to_disc:
        driver_type = 'MEM'
    else:
        driver_type = 'GTiff'

    out_ds = gdalwarp.crop_resample_average(input_tif=ifg.data_path,
                                            extents=extents,
                                            new_res=resolution,
                                            output_file=looks_path,
                                            thresh=thresh,
                                            out_driver_type=driver_type)[1]
    if not write_to_disc:
        return out_ds  # this is used in orbfit correction when looks!=1


def resample(data, xscale, yscale, thresh):
    """
    # TODO: make more efficient
    This is the slowest step in prepifg
    Resamples/averages 'data' to return an array from the averaging of blocks
    of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

    :param data: source array to resample to different size
    :param xscale: number of cells to average along X axis
    :param yscale: number of Y axis cells to average
    :param thresh: minimum allowable proportion of NaN cells (range from 0.0-1.0),
        eg. 0.25 = 1/4 or more as NaNs results in a NaN value for the output cell.
    """
    if thresh < 0 or thresh > 1:
        raise ValueError("threshold must be >= 0 and <= 1")

    xscale = int(xscale)
    yscale = int(yscale)
    ysize, xsize = data.shape
    xres, yres = (xsize / xscale), (ysize / yscale)
    dest = zeros((yres, xres), dtype=float32) * nan
    tile_cell_count = xscale * yscale

    # calc mean without nans (fractional threshold ignores tiles with excess NaNs)
    for x in xrange(xres):
        for y in xrange(yres):
            tile = data[y * yscale: (y+1) * yscale, x * xscale: (x+1) * xscale]
            nan_fraction = nsum(isnan(tile)) / float(tile_cell_count)
            if nan_fraction < thresh or (nan_fraction == 0 and thresh == 0):
                dest[y, x] = nanmean(tile)
    return dest


def reproject():
    raise NotImplementedError("TODO: Reprojection LOS/Horiz/Vert")


def check_resolution(ifgs):
    """
    Verifies Ifg resolutions are equal for the given grids.
    """

    for var in ['x_step', 'y_step']:
        values = array([getattr(i, var) for i in ifgs])
        if not (values == values[0]).all():
            msg = "Grid resolution does not match for %s" % var
            raise PreprocessError(msg)


def check_looks(xlooks, ylooks):
    """
    Verifies looks parameters are valid.
    """

    if not (isinstance(xlooks, Number) and isinstance(ylooks, Number)):
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise PreprocessError(msg)

    if not (xlooks > 0 and ylooks > 0):
        msg = "Invalid looks parameter(s), x: %s, y: %s. Looks must be an integer greater than zero" % (xlooks, ylooks)
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

    xmin, xmax = i.x_first, i.x_last
    ymin, ymax = i.y_first, i.y_last

    # swap y_first & y_last when using southern hemisphere -ve coords
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    return xmin, ymin, xmax, ymax


def custom_bounds(ifgs, xw, ytop, xe, ybot):
    """
    Check and modify input custom crop bounds to line up with grid interval
    """
    msg = 'Cropped image bounds exceed original image bounds'
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
        if (remainder > GRID_TOL) and (remainder < (1 - GRID_TOL)):
            msg = "%s crop extent not within %s of grid coordinate"
            raise PreprocessError(msg % (par, GRID_TOL))


class PreprocessError(Exception):
    pass


def extents_from_params(params):
    keys = (cfg.IFG_XFIRST, cfg.IFG_YFIRST, cfg.IFG_XLAST, cfg.IFG_YLAST)
    return CustomExts(*[params[k] for k in keys])
