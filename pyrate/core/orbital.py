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
This Python module implements residual orbital corrections for interferograms.
"""
from collections import OrderedDict

import numpy as np
from numpy import dot, vstack, zeros, meshgrid
from numpy import empty, isnan, reshape, float32, squeeze
from numpy.linalg import pinv
# from joblib import Parallel, delayed
from scipy.linalg import lstsq


from core.gdal_python import _crop_resample_setup, _setup_source, gdal_average, _alignment, coherence_masking


from core import shared, ifgconstants as ifc, config as cf, prepifg_helper, mst
from core.algorithm import master_slave_ids, get_all_epochs
from core.logger import pyratelogger as log
from core.shared import nanmedian, Ifg
from numpy import array, where, nan, isnan, nanmean, float32, zeros, \
    sum as nsum
from osgeo import gdal
import shutil
from core import ifgconstants as ifc, config as cf
from core.shared import Ifg, DEM, output_tiff_filename
import os
import shutil
from collections import namedtuple
from math import modf
from numbers import Number
from subprocess import check_call
from tempfile import mkstemp
import logging
# Constants
MINIMUM_CROP = 1
MAXIMUM_CROP = 2
CUSTOM_CROP = 3
ALREADY_SAME_SIZE = 4
CROP_OPTIONS = [MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE]

GRID_TOL = 1e-6


def crop_resample_average(
        input_tif, extents, new_res, output_file, thresh,
        out_driver_type='GTiff',
        match_pyrate=False, hdr=None, coherence_path=None, coherence_thresh=None):
    """
    Crop, resample, and average a geotiff image.
    :param str input_tif: Path to input geotiff to resample/crop
    :param tuple extents: Cropping extents (xfirst, yfirst, xlast, ylast)
    :param list new_res: [xres, yres] resolution of output image
    :param str output_file: Path to output resampled/cropped geotiff
    :param float thresh: NaN fraction threshold
    :param str out_driver_type: The output driver; `MEM` or `GTiff` (optional)
    :param bool match_pyrate: Match Legacy output (optional)
    :param dict hdr: dictionary of metadata
    :return: resampled_average: output cropped and resampled image
    :rtype: ndarray
    :return: out_ds: destination gdal dataset object
    :rtype: gdal.Dataset
    """
    dst_ds, _, _, _ = _crop_resample_setup(
        extents, input_tif, new_res, output_file,
        out_bands=2, dst_driver_type='MEM')

    # make a temporary copy of the dst_ds for PyRate style prepifg
    tmp_ds = gdal.GetDriverByName('MEM').CreateCopy('', dst_ds) \
        if (match_pyrate and new_res[0]) else None

    src_ds, src_ds_mem = _setup_source(input_tif)

    if coherence_path and coherence_thresh:
        coherence_raster = prepifg_helper.dem_or_ifg(coherence_path)
        coherence_raster.open()
        coherence_ds = coherence_raster.dataset
        coherence_masking(src_ds_mem, coherence_ds, coherence_thresh)
    elif coherence_path and not coherence_thresh:
        raise ValueError(f"Coherence file provided without a coherence "
                         f"threshold. Please ensure you provide 'cohthresh' "
                         f"in your config if coherence masking is enabled.")

    resampled_average, src_ds_mem = \
        gdal_average(dst_ds, src_ds, src_ds_mem, thresh)
    src_dtype = src_ds_mem.GetRasterBand(1).DataType
    src_gt = src_ds_mem.GetGeoTransform()

    # required to match Legacy output
    if tmp_ds:
        _alignment(input_tif, new_res, resampled_average, src_ds_mem,
                          src_gt, tmp_ds)

    # grab metadata from existing geotiff
    gt = dst_ds.GetGeoTransform()
    wkt = dst_ds.GetProjection()

    # TEST HERE IF EXISTING FILE HAS PYRATE METADATA. IF NOT ADD HERE
    if not ifc.DATA_TYPE in dst_ds.GetMetadata() and hdr is not None:
        md = shared.collate_metadata(hdr)
    else:
        md = dst_ds.GetMetadata()

    # update metadata for output
    for k, v in md.items():
        if k == ifc.DATA_TYPE:
            # update data type metadata
            if v == ifc.ORIG and coherence_path:
                md.update({ifc.DATA_TYPE:ifc.COHERENCE})
            elif v == ifc.ORIG and not coherence_path:
                md.update({ifc.DATA_TYPE:ifc.MULTILOOKED})
            elif v == ifc.DEM:
                md.update({ifc.DATA_TYPE:ifc.MLOOKED_DEM})
            elif v == ifc.INCIDENCE:
                md.update({ifc.DATA_TYPE:ifc.MLOOKED_INC})
            elif v == ifc.COHERENCE and coherence_path:
                pass
            elif v == ifc.MULTILOOKED and coherence_path:
                md.update({ifc.DATA_TYPE:ifc.COHERENCE})
            elif v == ifc.MULTILOOKED and not coherence_path:
                pass
            else:
                raise TypeError('Data Type metadata not recognised')

    # In-memory GDAL driver doesn't support compression so turn it off.
    creation_opts = ['compress=packbits'] if out_driver_type != 'MEM' else []
    out_ds = shared.gdal_dataset(output_file, dst_ds.RasterXSize, dst_ds.RasterYSize,
                                 driver=out_driver_type, bands=1, dtype=src_dtype, metadata=md, crs=wkt,
                                 geotransform=gt, creation_opts=creation_opts)

    shared.write_geotiff(resampled_average, out_ds, np.nan)

    return resampled_average, out_ds

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

def _min_bounds(ifgs):
    """
    Returns bounds for overlapping area of the given interferograms.
    """

    xmin = max([i.x_first for i in ifgs])
    ymax = min([i.y_first for i in ifgs])
    xmax = min([i.x_last for i in ifgs])
    ymin = max([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax


def _max_bounds(ifgs):
    """
    Returns bounds for the total area covered by the given interferograms.
    """

    xmin = min([i.x_first for i in ifgs])
    ymax = max([i.y_first for i in ifgs])
    xmax = max([i.x_last for i in ifgs])
    ymin = min([i.y_last for i in ifgs])
    return xmin, ymin, xmax, ymax

from decimal import Decimal
def _get_same_bounds(ifgs):
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

def get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts):
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


def _warp(ifg, x_looks, y_looks, extents, resolution, thresh, crop_out,
          write_to_disk=True, out_path=None, header=None,
          coherence_path=None, coherence_thresh=None):
    """
    Convenience function for calling GDAL functionality
    """
    if x_looks != y_looks:
        raise ValueError('X and Y looks mismatch')

    # cut, average, resample the final output layers
    op = output_tiff_filename(ifg.data_path, out_path)
    looks_path = cf.mlooked_path(op, y_looks, crop_out)

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
    resampled_data, out_ds = crop_resample_average( input_tif=ifg.data_path, extents=extents, new_res=resolution, output_file=looks_path, thresh=thresh, out_driver_type=driver_type, hdr=header, coherence_path=coherence_path, coherence_thresh=coherence_thresh)
    if not write_to_disk:
        return resampled_data, out_ds

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

def _dummy_warp(renamed_path):
    """
    Convenience dummy operation for when no multi-looking or cropping
    required
    """
    ifg = dem_or_ifg(renamed_path)
    ifg.open()
    ifg.dataset.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)
    data = ifg.dataset.ReadAsArray()
    return data, ifg.dataset

def prepare_ifg(raster_path, xlooks, ylooks, exts, thresh, crop_opt, write_to_disk=True, out_path=None, header=None, coherence_path=None, coherence_thresh=None):
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
    if not do_multilook and crop_opt == ALREADY_SAME_SIZE:
        renamed_path = cf.mlooked_path(raster.data_path, looks=xlooks, crop_out=crop_opt)
        shutil.copy(raster.data_path, renamed_path)
        # set metadata to indicated has been cropped and multilooked
        # copy file with mlooked path
        return _dummy_warp(renamed_path)

    return _warp(raster, xlooks, ylooks, exts, resolution, thresh, crop_opt, write_to_disk, out_path, header, coherence_path, coherence_thresh)

# TODO: crop options 0 = no cropping? get rid of same size
def prepare_ifgs(raster_data_paths, crop_opt, xlooks, ylooks, thresh=0.5, user_exts=None, write_to_disc=True,
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
    :param tuple user_exts: Tuple of user defined georeferenced extents for new file: (xfirst, yfirst, xlast, ylast)cropping coordinates
    :param bool write_to_disk: Write new data to disk
    :return: resampled_data: output cropped and resampled image
    :rtype: ndarray
    :return: out_ds: destination gdal dataset object
    :rtype: gdal.Dataset
    """
    # use metadata check to check whether it's a dem or ifg
    rasters = [dem_or_ifg(r) for r in raster_data_paths]
    exts = get_analysis_extent(crop_opt, rasters, xlooks, ylooks, user_exts)

    return [prepare_ifg(d, xlooks, ylooks, exts, thresh, crop_opt, write_to_disc, out_path) for d in raster_data_paths]


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

# Orbital correction tasks
#
# TODO: options for multilooking
# 1) do the 2nd stage mlook at prepifg.py/generate up front, then delete in
#    workflow afterward
# 2) refactor prep_ifgs() call to take input filenames and params & generate
#    mlooked versions from that
#    this needs to be more generic to call at any point in the runtime.


# Design notes:
# The orbital correction code includes several enhancements. PyRate creates
# sparse arrays for the linear inversion, which contain many empty cells.
# This is unnecessary for the independent method, and temporarily wastes
# potentially a lot of memory.
#
# For the independent method, PyRate makes individual small design matrices and
# corrects the Ifgs one by one. If required in the correction, the offsets
# option adds an extra column of ones to include in the inversion.
#
# Network method design matrices are mostly empty, and offsets are handled
# differently. Individual design matrices (== independent method DMs) are
# placed in the sparse network design matrix. Offsets are not included in the
# smaller DMs to prevent unwanted cols being inserted. This is why some funcs
# appear to ignore the offset parameter in the networked method. Network DM
# offsets are cols of 1s in a diagonal line on the LHS of the sparse array.

# ORBITAL ERROR correction constants
INDEPENDENT_METHOD = cf.INDEPENDENT_METHOD
NETWORK_METHOD = cf.NETWORK_METHOD

PLANAR = cf.PLANAR
QUADRATIC = cf.QUADRATIC
PART_CUBIC = cf.PART_CUBIC


def remove_orbital_error(ifgs, params, preread_ifgs=None):
    """Wrapper function for PyRate orbital error removal functionality.
    
    NB: the ifg data is modified in situ, rather than create intermediate
    files. The network method assumes the given ifgs have already been reduced
    to a minimum spanning tree network.
    
    Args:
        ifgs (list): List of interferograms class objects
    
    Args:
      preread_ifgs: dict (Default value = None)
      MPI: jobs
      ifgs: param params:

    Args:
      ifgs: 
      params: 
      preread_ifgs:  (Default value = None)

    Returns:
      None - interferogram phase data is updated and saved to disk

    """

    ifg_paths = [i.data_path for i in ifgs] if isinstance(ifgs[0], Ifg) else ifgs

    mlooked = None
    # mlooking is not necessary for independent correction
    # can use multiple procesing if write_to_disc=True
    if params[cf.ORBITAL_FIT_METHOD] == NETWORK_METHOD:
        mlooked_dataset = prepare_ifgs(
            ifg_paths,
            crop_opt=prepifg_helper.ALREADY_SAME_SIZE,
            xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
            ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
            thresh=params[cf.NO_DATA_AVERAGING_THRESHOLD],
            write_to_disc=False,
        )
        mlooked = [Ifg(m[1]) for m in mlooked_dataset]

        for m in mlooked:
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]
            m.convert_to_nans()
            m.convert_to_mm()

    _orbital_correction(ifgs, params, mlooked=mlooked, preread_ifgs=preread_ifgs)


def _orbital_correction(ifgs_or_ifg_paths, params, mlooked=None, offset=True, preread_ifgs=None):
    """Convenience function to perform orbital correction.
    
    Args:
        ifgs_or_ifg_paths:
    
    Args:
      mlooked: Default value = None)
      offset: Default value = True)
      preread_ifgs: Default value = None)
      ifgs_or_ifg_paths: param params:

    Args:
      ifgs_or_ifg_paths: 
      params: 
      mlooked:  (Default value = None)
      offset:  (Default value = True)
      preread_ifgs:  (Default value = None)

    Returns:
      

    """
    degree = params[cf.ORBITAL_FIT_DEGREE]
    method = params[cf.ORBITAL_FIT_METHOD]
    # parallel = params[cf.PARALLEL]  # not implemented

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        msg = "Invalid degree of %s for orbital correction" % cf.ORB_DEGREE_NAMES.get(degree)
        raise OrbitalError(msg)

    log.info(
        "Removing orbital error using {} correction method" " and degree={}".format(cf.ORB_METHOD_NAMES.get(method),
                                                                                    cf.ORB_DEGREE_NAMES.get(degree))
    )

    if method == NETWORK_METHOD:
        if mlooked is None:
            network_orbital_correction(ifgs_or_ifg_paths, degree, offset, params, m_ifgs=mlooked,
                                       preread_ifgs=preread_ifgs)
        else:
            _validate_mlooked(mlooked, ifgs_or_ifg_paths)
            network_orbital_correction(ifgs_or_ifg_paths, degree, offset, params, mlooked, preread_ifgs)

    elif method == INDEPENDENT_METHOD:
        # not running in parallel
        # raises swig object pickle error
        # Parallel(n_jobs=params[cf.PROCESSES],
        #          verbose=joblib_log_level(cf.LOG_LEVEL))(
        #     delayed(_independent_correction)(ifg, degree, offset, params)
        #     for ifg in ifgs)
        for ifg in ifgs_or_ifg_paths:
            independent_orbital_correction(ifg, degree, offset, params)
    else:
        msg = "Unknown method: '%s', need INDEPENDENT or NETWORK method"
        raise OrbitalError(msg % method)


def _validate_mlooked(mlooked, ifgs):
    """Basic sanity checking of the multilooked ifgs.

    Args:
      mlooked: param ifgs:
      ifgs: 

    Returns:

    """

    if len(mlooked) != len(ifgs):
        msg = "Mismatching # ifgs and # multilooked ifgs"
        raise OrbitalError(msg)

    if not all([hasattr(i, "phase_data") for i in mlooked]):
        msg = "Mismatching types in multilooked ifgs arg:\n%s" % mlooked
        raise OrbitalError(msg)


def _get_num_params(degree, offset=None):
    """Returns number of model parameters from string parameter

    Args:
      degree: param offset:  (Default value = None)
      offset: (Default value = None)

    Returns:

    """

    if degree == PLANAR:
        nparams = 2
    elif degree == QUADRATIC:
        nparams = 5
    elif degree == PART_CUBIC:
        nparams = 6
    else:
        msg = "Invalid orbital model degree: %s" % cf.ORB_DEGREE_NAMES.get(degree)
        raise OrbitalError(msg)

    # NB: independent method only, network method handles offsets separately
    if offset is True:
        nparams += 1  # eg. y = mx + offset
    return nparams


def independent_orbital_correction(ifg, degree, offset, params):
    """Calculates and removes an orbital error surface from a single independent
    interferogram.
    
    Warning: This will write orbital error corrected phase_data to the ifg.
    
    Args:
        ifg (Ifg class instance): the interferogram to be corrected
        degree (str): model to fit (PLANAR / QUADRATIC / PART_CUBIC)
        offset (bool): True to calculate the model using an offset
    
    Args:
      ifg: param degree:
      offset: param params:
      degree:

    Args:
      ifg: 
      degree: 
      offset: 
      params: 

    Returns:
      None - interferogram phase data is updated and saved to disk

    """
    ifg = shared.Ifg(ifg) if isinstance(ifg, str) else ifg
    if not ifg.is_open:
        ifg.open()
    shared.nan_and_mm_convert(ifg, params)
    # vectorise, keeping NODATA
    vphase = reshape(ifg.phase_data, ifg.num_cells)
    dm = get_design_matrix(ifg, degree, offset)

    # filter NaNs out before getting model
    clean_dm = dm[~isnan(vphase)]
    data = vphase[~isnan(vphase)]
    model = lstsq(clean_dm, data)[0]  # first arg is the model params

    # calculate forward model & morph back to 2D
    if offset:
        fullorb = np.reshape(np.dot(dm[:, :-1], model[:-1]), ifg.phase_data.shape)
    else:
        fullorb = np.reshape(np.dot(dm, model), ifg.phase_data.shape)
    offset_removal = nanmedian(np.ravel(ifg.phase_data - fullorb))
    # subtract orbital error from the ifg
    ifg.phase_data -= fullorb - offset_removal
    # set orbfit meta tag and save phase to file
    _save_orbital_error_corrected_phase(ifg)
    if ifg.open():
        ifg.close()


def network_orbital_correction(ifgs, degree, offset, params, m_ifgs=None, preread_ifgs=None):
    """This algorithm implements a network inversion to determine orbital
    corrections for a set of interferograms forming a connected network.
    
    Warning: This will write orbital error corrected phase_data to the ifgs.
    
    Args:
        ifgs (list): List of Ifg class objects reduced to a minimum spanning
            tree network
        degree (str): model to fit (PLANAR / QUADRATIC / PART_CUBIC)
        offset (bool): True to calculate the model using offsets
    
    Args:
      m_ifgs: list (Default value = None)
      multilooked: versions of
      preread_ifgs: dict (Default value = None)
      MPI: jobs
      ifgs: param degree:
      offset: param params:
      degree:

    Args:
      ifgs: 
      degree: 
      offset: 
      params: 
      m_ifgs:  (Default value = None)
      preread_ifgs:  (Default value = None)

    Returns:
      None - interferogram phase data is updated and saved to disk

    """
    src_ifgs = ifgs if m_ifgs is None else m_ifgs
    src_ifgs = mst.mst_from_ifgs(src_ifgs)[3]  # use networkx mst

    vphase = vstack([i.phase_data.reshape((i.num_cells, 1)) for i in src_ifgs])
    vphase = squeeze(vphase)

    B = get_network_design_matrix(src_ifgs, degree, offset)

    # filter NaNs out before getting model
    B = B[~isnan(vphase)]
    orbparams = dot(pinv(B, 1e-6), vphase[~isnan(vphase)])

    ncoef = _get_num_params(degree)
    if preread_ifgs:
        temp_ifgs = OrderedDict(sorted(preread_ifgs.items())).values()
        ids = master_slave_ids(get_all_epochs(temp_ifgs))
    else:
        ids = master_slave_ids(get_all_epochs(ifgs))
    coefs = [orbparams[i: i + ncoef] for i in range(0, len(set(ids)) * ncoef, ncoef)]

    # create full res DM to expand determined coefficients into full res
    # orbital correction (eg. expand coarser model to full size)

    if preread_ifgs:
        temp_ifg = Ifg(ifgs[0])  # ifgs here are paths
        temp_ifg.open()
        dm = get_design_matrix(temp_ifg, degree, offset=False)
        temp_ifg.close()
    else:
        dm = get_design_matrix(ifgs[0], degree, offset=False)

    for i in ifgs:
        # open if not Ifg instance
        if isinstance(i, str):  # pragma: no cover
            # are paths
            i = Ifg(i)
            i.open(readonly=False)
            shared.nan_and_mm_convert(i, params)
        _remove_network_orb_error(coefs, dm, i, ids, offset)


def _remove_network_orb_error(coefs, dm, ifg, ids, offset):
    """Convenience function to remove network orbital error from input
    interferograms

    Args:
      coefs: param dm:
      ifg: param ids:
      offset: 
      dm: 
      ids: 

    Returns:

    """
    orb = dm.dot(coefs[ids[ifg.slave]] - coefs[ids[ifg.master]])
    orb = orb.reshape(ifg.shape)
    # offset estimation
    if offset:
        # bring all ifgs to same base level
        orb -= nanmedian(np.ravel(ifg.phase_data - orb))
    # subtract orbital error from the ifg
    ifg.phase_data -= orb
    # set orbfit meta tag and save phase to file
    _save_orbital_error_corrected_phase(ifg)


def _save_orbital_error_corrected_phase(ifg):
    """Convenience function to update metadata and save latest phase after
    orbital fit correction

    Args:
      ifg: 

    Returns:

    """
    # set orbfit tags after orbital error correction
    ifg.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR, ifc.ORB_REMOVED)
    ifg.write_modified_phase()
    ifg.close()


# TODO: subtract reference pixel coordinate from x and y
def get_design_matrix(ifg, degree, offset, scale=100.0):
    """Returns orbital error design matrix with columns for model parameters.

    Args:
      ifg(Ifg class instance): interferogram to get design matrix for
      degree(str): model to fit (PLANAR / QUADRATIC / PART_CUBIC)
      offset(bool): True to include offset column, otherwise False.
      scale(float, optional): Scale factor to divide cell size by in order to improve
    inversion robustness (Default value = 100.0)

    Returns:
      ndarray: dm: design matrix

    """

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        raise OrbitalError("Invalid degree argument")

    # scaling required with higher degree models to help with param estimation
    xsize = ifg.x_size / scale if scale else ifg.x_size
    ysize = ifg.y_size / scale if scale else ifg.y_size

    # mesh needs to start at 1, otherwise first cell resolves to 0 and ignored
    xg, yg = [g + 1 for g in meshgrid(range(ifg.ncols), range(ifg.nrows))]
    x = xg.reshape(ifg.num_cells) * xsize
    y = yg.reshape(ifg.num_cells) * ysize

    # TODO: performance test this vs np.concatenate (n by 1 cols)??
    dm = empty((ifg.num_cells, _get_num_params(degree, offset)), dtype=float32)

    # apply positional parameter values, multiply pixel coordinate by cell size
    # to get distance (a coord by itself doesn't tell us distance from origin)
    if degree == PLANAR:
        dm[:, 0] = x
        dm[:, 1] = y
    elif degree == QUADRATIC:
        dm[:, 0] = x ** 2
        dm[:, 1] = y ** 2
        dm[:, 2] = x * y
        dm[:, 3] = x
        dm[:, 4] = y
    elif degree == PART_CUBIC:
        dm[:, 0] = x * (y ** 2)
        dm[:, 1] = x ** 2
        dm[:, 2] = y ** 2
        dm[:, 3] = x * y
        dm[:, 4] = x
        dm[:, 5] = y
    if offset is True:
        dm[:, -1] = np.ones(ifg.num_cells)

    return dm


def get_network_design_matrix(ifgs, degree, offset):
    """Returns larger-format design matrix for network error correction. The
    network design matrix includes rows which relate to those of NaN cells.

    Args:
      ifgs(list): List of Ifg class objects
      degree(str): model to fit (PLANAR / QUADRATIC / PART_CUBIC)
      offset(bool): True to include offset cols, otherwise False.

    Returns:
      ndarray: netdm: network design matrix

    """

    if degree not in [PLANAR, QUADRATIC, PART_CUBIC]:
        raise OrbitalError("Invalid degree argument")

    nifgs = len(ifgs)
    if nifgs < 1:
        # can feasibly do correction on a single Ifg/2 epochs
        raise OrbitalError("Invalid number of Ifgs: %s" % nifgs)

    # init sparse network design matrix
    nepochs = len(set(get_all_epochs(ifgs)))

    # no offsets: they are made separately below
    ncoef = _get_num_params(degree)
    shape = [ifgs[0].num_cells * nifgs, ncoef * nepochs]

    if offset:
        shape[1] += nifgs  # add extra block for offset cols

    netdm = zeros(shape, dtype=float32)

    # calc location for individual design matrices
    dates = [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]
    ids = master_slave_ids(dates)
    offset_col = nepochs * ncoef  # base offset for the offset cols
    tmpdm = get_design_matrix(ifgs[0], degree, offset=False)

    # iteratively build up sparse matrix
    for i, ifg in enumerate(ifgs):
        rs = i * ifg.num_cells  # starting row
        m = ids[ifg.master] * ncoef  # start col for master
        s = ids[ifg.slave] * ncoef  # start col for slave
        netdm[rs: rs + ifg.num_cells, m: m + ncoef] = -tmpdm
        netdm[rs: rs + ifg.num_cells, s: s + ncoef] = tmpdm

        # offsets are diagonal cols across the extra array block created above
        if offset:
            netdm[rs: rs + ifg.num_cells, offset_col + i] = 1  # init offset cols

    return netdm


class OrbitalError(Exception):
    """Generic class for errors in orbital correction."""
