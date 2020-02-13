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
"""
This Python module converts interferogram input files to a common geotiff
format with PyRate specific metadata headers. The module also implements
multilooking/downsampling and cropping operations to reduce the size of
the computational problem.
"""
import os
import shutil
from collections import namedtuple
from math import modf
from numbers import Number
from subprocess import check_call
from tempfile import mkstemp

from numpy import array, where, nan, isnan, nanmean, float32, zeros, sum as nsum
from osgeo import gdal

from core import ifgconstants as ifc
from core.gdal_python import crop_resample_average
from core.logger import pyratelogger as log
from core.shared import Ifg, DEM

CustomExts = namedtuple("CustExtents", ["xfirst", "yfirst", "xlast", "ylast"])

# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
ALREADY_SAME_SIZE = 4
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE]

GRID_TOL = 1e-6


def get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts):
    """Function checks prepifg parameters and returns extents/bounding box.

    Args:
      crop_opt(int): Cropping option
      rasters(list): List of either Ifg or DEM class objects
      xlooks(int): Number of multi-looks in x
      ylooks(int): Number of multi-looks in y
      user_exts(tuple): Tuple of user defined cropping coordinates

    Returns:
      tuple: extents: tuple of four bounding coordinates

    """

    if crop_opt not in CROP_OPTIONS:
        raise PreprocessError("Unrecognised crop option: %s" % crop_opt)

    if crop_opt == CUSTOM_CROP:
        if not user_exts:
            raise PreprocessError("No custom cropping extents specified")
        elif len(user_exts) != 4:  # check for required numbers
            raise PreprocessError("Custom extents must have all 4 values")
        elif len(user_exts) == 4:  # check for non floats
            if not all([_is_number(z) for z in user_exts]):
                raise PreprocessError("Custom extents must be 4 numbers")

    for raster in rasters:
        if not raster.is_open:
            raster.open()

    _check_looks(xlooks, ylooks)
    _check_resolution(rasters)

    return _get_extents(rasters, crop_opt, user_exts)


def _is_number(s):
    """Check whether string can be converted to float

    Args:
      s: 

    Returns:

    """
    try:
        float(s)
        return True
    except TypeError:  # for example if 'None' is sent
        return False
    except ValueError:  # for example if a 'string' is sent
        return False


def _check_looks(xlooks, ylooks):
    """Convenience function to verify that looks parameters are valid.

    Args:
      xlooks: param ylooks:
      ylooks: 

    Returns:

    """

    if not (isinstance(xlooks, Number) and isinstance(ylooks, Number)):  # pragma: no cover
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise PreprocessError(msg)

    if not (xlooks > 0 and ylooks > 0):  # pragma: no cover
        msg = "Invalid looks parameter(s), x: %s, y: %s. " "Looks must be an integer greater than zero" % (
            xlooks, ylooks)
        raise PreprocessError(msg)


def _check_resolution(ifgs):
    """Convenience function to verify Ifg resolutions are equal.

    Args:
      ifgs: 

    Returns:

    """

    for var in ["x_step", "y_step"]:
        values = array([getattr(i, var) for i in ifgs])
        if not (values == values[0]).all():  # pragma: no cover
            msg = "Grid resolution does not match for %s" % var
            raise PreprocessError(msg)


def _get_extents(ifgs, crop_opt, user_exts=None):
    """Convenience function that returns extents/bounding box. MINIMUM_CROP = 1
    MAXIMUM_CROP = 2 CUSTOM_CROP = 3 ALREADY_SAME_SIZE = 4

    Args:
      ifgs: param crop_opt:
      user_exts: Default value = None)
      crop_opt: 

    Returns:

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


def prepare_ifg(input_path, output_path, xlooks, ylooks, extents, thresh, crop_out, header=None, coherence_path=None,
                coherence_thresh=None):
    """Open, resample, crop and optionally save to disk an interferogram or DEM.
    Returns are only given if write_to_disk=False

    Args:
      input_path(str): Input raster file path name
      output_path(str): Path for output file
      xlooks(int): Number of multi-looks in x; 5 is 5 times smaller, 1 is no
    change
      ylooks(int): Number of multi-looks in y
      extents(tuple): Tuple of user defined georeferenced extents for new
    file: (xfirst, yfirst, xlast, ylast)cropping coordinates
      thresh(float): see thresh in prepare_ifgs()
      crop_out(int): Crop option
      header(dict, optional): dictionary of metadata from header file (Default value = None)
      coherence_path: Default value = None)
      coherence_thresh: Default value = None)

    Returns:
      ndarray: resampled_data: output cropped and resampled image

    """

    do_multilook = xlooks > 1 or ylooks > 1
    # resolution=None completes faster for non-multilooked layers in gdalwarp
    resolution = [None, None]
    input_raster = dem_or_ifg(input_path)

    if not input_raster.is_open:
        input_raster.open()

    if do_multilook:
        log.debug("xlooks: " + str(xlooks))
        log.debug("input_raster.x_step: " + str(input_raster.x_step))
        log.debug("ylooks: " + str(ylooks))
        log.debug("input_raster.y_step" + str(input_raster.y_step))

        resolution = [xlooks * input_raster.x_step, ylooks * input_raster.y_step]

    if not do_multilook and crop_out == ALREADY_SAME_SIZE:
        shutil.copy(input_raster.data_path, output_path)
        # set metadata to indicated has been cropped and multilooked
        # copy file with mlooked path
        _dummy_warp(output_path)

    if xlooks != ylooks:
        raise ValueError("X and Y looks mismatch")

    #     # Add missing/updated metadata to resampled ifg/DEM
    #     new_lyr = type(ifg)(looks_path)
    #     new_lyr.open(readonly=True)
    #     # for non-DEMs, phase bands need extra metadata & conversions
    #     if hasattr(new_lyr, "phase_band"):
    #         # TODO: LOS conversion to vertical/horizontal (projection)
    #         # TODO: push out to workflow
    #         #if params.has_key(REPROJECTION_FLAG):
    #         #    reproject()

    crop_resample_average(
        input_raster=input_raster.data_path,
        extents=extents,
        resolution=resolution,
        output_file=output_path,
        thresh=thresh,
        header=header,
        coherence_path=coherence_path,
        coherence_thresh=coherence_thresh,
    )


# TODO: crop options 0 = no cropping? get rid of same size
def prepare_ifgs(raster_data_paths, crop_opt, xlooks, ylooks, thresh=0.5, user_exts=None, write_to_disc=True,
                 out_path=None):
    """Wrapper function to prepare a sequence of interferogram files for PyRate
    analysis. See prepifg.prepare_ifg() for full description of inputs and
    returns.
    
    Note: function need refining for crop options

    Args:
      raster_data_paths(list): List of interferogram file paths
      crop_opt(int): Crop option
      xlooks(int): Number of multi-looks in x; 5 is 5 times smaller, 1 is no
    change
      ylooks(int): Number of multi-looks in y
      thresh(float, optional): see thresh in prepare_ifgs() (Default value = 0.5)
      user_exts(tuple, optional): Tuple of user defined georeferenced extents for new
    file: (xfirst, yfirst, xlast, ylast)cropping coordinates (Default value = None)
      write_to_disc: Default value = True)
      out_path: Default value = None)

    Returns:
      ndarray: resampled_data: output cropped and resampled image

    """
    # use metadata check to check whether it's a dem or ifg
    rasters = [dem_or_ifg(r) for r in raster_data_paths]
    exts = get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts)

    return [prepare_ifg(d, xlooks, ylooks, exts, thresh, crop_opt, write_to_disc, out_path) for d in raster_data_paths]


def dem_or_ifg(data_path):
    """Returns an Ifg or DEM class object from input geotiff file.

    Args:
      data_path(str): file path name

    Returns:
      Ifg or DEM class object: Interferogram or DEM object from input file

    """
    ds = gdal.Open(data_path)
    md = ds.GetMetadata()
    if ifc.MASTER_DATE in md:  # ifg
        return Ifg(data_path)
    else:
        return DEM(data_path)


# TODO: Not currently used; remove in future?
def _file_ext(raster):
    """Returns file ext string based on type of raster.

    Args:
      raster: 

    Returns:

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
        raise NotImplementedError("Missing raster types for LOS, " "Coherence and baseline")


def _dummy_warp(renamed_path):
    """Convenience dummy operation for when no multi-looking or cropping
    required

    Args:
      renamed_path: 

    Returns:

    """
    ifg = dem_or_ifg(renamed_path)
    ifg.open()
    ifg.dataset.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)
    data = ifg.dataset.ReadAsArray()
    return data, ifg.dataset


# TODO: Not being used. Remove in future?
def _resample(data, xscale, yscale, thresh):
    """Resamples/averages 'data' to return an array from the averaging of blocks
    of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

    Args:
      data: source array to resample to different size
      xscale: number of cells to average along X axis
      yscale: number of Y axis cells to average
      thresh: minimum allowable proportion of NaN cells (range from 0.0-1.0),
    eg. 0.25 = 1/4 or more as NaNs results in a NaN value for the output
    cell.

    Returns:

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
            tile = data[y * yscale: (y + 1) * yscale, x * xscale: (x + 1) * xscale]
            nan_fraction = nsum(isnan(tile)) / float(tile_cell_count)
            if nan_fraction < thresh or (nan_fraction == 0 and thresh == 0):
                dest[y, x] = nanmean(tile)
    return dest


# TODO: Not being used. Remove in future?
def _resample_ifg(ifg, cmd, x_looks, y_looks, thresh, md=None):
    """Convenience function to resample data from a given Ifg (more coarse).

    Args:
      ifg: param cmd:
      x_looks: param y_looks:
      thresh: param md:  (Default value = None)
      cmd: 
      y_looks: 
      md: (Default value = None)

    Returns:

    """

    fp, tmp_path = mkstemp(suffix=".tif")
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
    return _resample(data, x_looks, y_looks, thresh)


def _min_bounds(ifgs):
    """Returns bounds for overlapping area of the given interferograms.

    Args:
      ifgs: 

    Returns:

    """

    xmin = max([i.x_first for i in ifgs])
    ymax = min([i.y_first for i in ifgs])
    xmax = min([i.x_last for i in ifgs])
    ymin = max([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def _max_bounds(ifgs):
    """Returns bounds for the total area covered by the given interferograms.

    Args:
      ifgs: 

    Returns:

    """

    xmin = min([i.x_first for i in ifgs])
    ymax = max([i.y_first for i in ifgs])
    xmax = max([i.x_last for i in ifgs])
    ymin = min([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


from decimal import Decimal


def _get_same_bounds(ifgs):
    """Check and return bounding box for ALREADY_SAME_SIZE option.

    Args:
      ifgs: 

    Returns:

    """

    tfs = [i.dataset.GetGeoTransform() for i in ifgs]

    equal = []

    for t in tfs[1:]:
        for i, tf in enumerate(tfs[0]):

            if round(Decimal(tf), 4) == round(Decimal(t[i]), 4):
                equal.append(True)
            else:
                equal.append(False)

    if not all(equal):
        msg = "Ifgs do not have the same bounding box for crop option: %s"
        raise PreprocessError(msg % ALREADY_SAME_SIZE)
    ifg = ifgs[0]
    xmin, xmax = ifg.x_first, ifg.x_last
    ymin, ymax = ifg.y_first, ifg.y_last

    # swap y_first & y_last when using southern hemisphere -ve coords
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    return xmin, ymin, xmax, ymax


def _custom_bounds(ifgs, xw, ytop, xe, ybot):
    """Check and modify input custom crop bounds to line up with grid interval

    Args:
      ifgs: param xw:
      ytop: param xe:
      ybot: 
      xw: 
      xe: 

    Returns:

    """
    msg = "Cropped image bounds are outside the original image bounds"
    i = ifgs[0]

    if ytop < ybot:
        raise PreprocessError("ERROR Custom crop bounds: " "ifgyfirst must be greater than ifgylast")

    if xe < xw:
        raise PreprocessError("ERROR Custom crop bounds: " "ifgxfirst must be greater than ifgxlast")

    for par, crop, orig, step in zip(
            ["x_first", "x_last", "y_first", "y_last"],
            [xw, xe, ytop, ybot],
            [i.x_first, i.x_last, i.y_first, i.y_last],
            [i.x_step, i.x_step, i.y_step, i.y_step],
    ):
        diff = crop - orig
        nint = round(diff / step)

        if par == "x_first":
            if diff < 0:
                raise PreprocessError(msg)
            xmin = orig + (nint * step)

        elif par == "x_last":
            if diff > 0:
                raise PreprocessError(msg)
            xmax = orig + (nint * step)

        elif par == "y_first":
            if diff > 0:
                raise PreprocessError(msg)
            y1 = orig + (nint * step)

        elif par == "y_last":
            if diff < 0:
                raise PreprocessError(msg)
            y2 = orig + (nint * step)
        else:
            raise ValueError("Value error in supplied custom bounds")

    if y2 > y1:
        ymin = y1
        ymax = y2
    else:
        ymin = y2
        ymax = y1

    return xmin, ymin, xmax, ymax


def _check_crop_coords(ifgs, xmin, ymin, xmax, ymax):
    """Ensures cropping coords line up with grid system within tolerance.

    Args:
      ifgs: param xmin:
      ymin: param xmax:
      ymax: 
      xmin: 
      xmax: 

    Returns:

    """

    # NB: assumption is the first Ifg is correct, so only test against it
    i = ifgs[0]

    for par, crop, step in zip(["x_first", "x_last", "y_first", "y_last"], [xmin, xmax, ymax, ymin],
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
    """Preprocess exception"""
