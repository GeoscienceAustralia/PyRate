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
This Python script applies optional multilooking and cropping to input
interferogram geotiff files.
"""
# -*- coding: utf-8 -*-
import os
import shutil
from decimal import Decimal
from math import modf

import numpy as np
from osgeo import gdal
from osgeo import gdalconst


from constants import (
    CROP_OPTIONS,
    GRID_TOL,
    GDAL_X_CELLSIZE,
    GDAL_Y_CELLSIZE,
    GDAL_X_FIRST,
    GDAL_Y_FIRST,
    MINIMUM_CROP,
    MAXIMUM_CROP,
    CUSTOM_CROP,
    ALREADY_SAME_SIZE,
)
from core import shared, config as cf, gamma, roipac, ifgconstants as ifc
from core.gdal_python import _crop_resample_setup, _setup_source, coherence_masking, gdal_average, _alignment
from core.logger import pyratelogger as log
from core.mpiops import rank, comm, size, chunks
from core.shared import Ifg, DEM

GAMMA = 1
ROIPAC = 0


def main(params):
    """Main workflow function for preparing interferograms for PyRate.
    Args:
      params: 

    Returns:
    """
    if rank == 0:
        log.info("Collecting jobs: prepifg")
        # a job is list of parameters passed to multiprocessing function
        jobs = []

        user_extents = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST], params[cf.IFG_XLAST], params[cf.IFG_YLAST])

        datasets_to_calcualte_extents = []
        for interferogram_file in params["interferogram_files"]:
            if not os.path.isfile(interferogram_file.converted_path):
                raise Exception("Can not find geotiff: " + str(interferogram_file.converted_path) +
                                ". Ensure you have converted your interferograms to geotiffs.")
            datasets_to_calcualte_extents.append(interferogram_file.converted_path)

        # optional DEM conversion
        if params["dem_file"] is not None:
            if not os.path.isfile(params["dem_file"].converted_path):
                raise Exception("Can not find geotiff: " + str(params["dem_file"].converted_path) +
                                ". Ensure you have converted your interferograms to geotiffs.")
            datasets_to_calcualte_extents.append(params["dem_file"].converted_path)

        extents = get_analysis_extent(datasets_to_calcualte_extents, params["ifgcropopt"], params["ifglksx"], params["ifglksy"], user_exts=user_extents)

        log.info("Preparing interferograms by cropping/multilooking")
        for interferogram_file in params["interferogram_files"]:
            jobs.append((interferogram_file.converted_path, interferogram_file.sampled_path, extents, params, "interferogram"))

        # optional DEM conversion
        if params["dem_file"] is not None:
            jobs.append((params["dem_file"].converted_path, params["dem_file"].sampled_path, extents, params, "dem"))
        jobs = chunks(jobs, size)
    else:
        jobs = None

    jobs = comm.scatter(jobs, root=0)

    for job in jobs:
        _prepifg_multiprocessing(*job)


def _prepifg_multiprocessing(input_path, output_path, extents, params, tag):

    xlooks = params["ifglksx"]
    ylooks = params["ifglksy"]
    crop_out = params["ifgcropopt"]
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    processor = params[cf.PROCESSOR]  # roipac or gamma

    if processor == GAMMA:
        header = gamma.gamma_header(input_path, params)
    elif processor == ROIPAC:
        header = roipac.roipac_header(input_path, params)
    else:
        raise Exception("Processor must be ROI_PAC (0) or GAMMA (1)")

    do_multilook = xlooks > 1 or ylooks > 1
    # resolution=None completes faster for non-multilooked layers in gdalwarp
    resolution = [None, None]

    if "dem" in tag:
        input_raster = DEM(input_path)
    else:
        input_raster = Ifg(input_path)

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
        input_raster.dataset.SetMetadataItem(ifc.DATA_TYPE, ifc.MULTILOOKED)


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
    input_raster = input_raster.data_path
    output_file = output_path

    out_driver_type = "GTiff"
    match_pyrate = False
    dst_ds, _, _, _ = _crop_resample_setup(extents, input_raster, resolution, output_file, out_bands=2, dst_driver_type="MEM")

    # make a temporary copy of the dst_ds for PyRate style prepifg
    tmp_ds = gdal.GetDriverByName("MEM").CreateCopy("", dst_ds) if (match_pyrate and resolution[0]) else None

    src_ds, src_ds_mem = _setup_source(input_raster)

    if params["cohmask"] and "interferogram" in tag:
        if params["coherence_file_paths"] and params["cohthresh"]:
            coherence_masking(src_ds_mem, input_path, params["coherence_file_paths"], params["cohthresh"])

    resampled_average, src_ds_mem = gdal_average(dst_ds, src_ds, src_ds_mem, thresh)
    src_dtype = src_ds_mem.GetRasterBand(1).DataType
    src_gt = src_ds_mem.GetGeoTransform()

    # required to match Legacy output
    if tmp_ds:
        _alignment(input_raster, resolution, resampled_average, src_ds_mem, src_gt, tmp_ds)

    # grab metadata from existing geotiff
    gt = dst_ds.GetGeoTransform()
    wkt = dst_ds.GetProjection()

    # TEST HERE IF EXISTING FILE HAS PYRATE METADATA. IF NOT ADD HERE
    if not ifc.DATA_TYPE in dst_ds.GetMetadata() and header is not None:
        md = shared.collate_metadata(header)
    else:
        md = dst_ds.GetMetadata()

    # update metadata for output

    for k, v in md.items():
        if k == ifc.DATA_TYPE:
            # update data type metadata
            if v == ifc.ORIG and params["coherence_file_paths"]:
                md.update({ifc.DATA_TYPE: ifc.COHERENCE})
            elif v == ifc.ORIG and not params["coherence_file_paths"]:
                md.update({ifc.DATA_TYPE: ifc.MULTILOOKED})
            elif v == ifc.DEM:
                md.update({ifc.DATA_TYPE: ifc.MLOOKED_DEM})
            elif v == ifc.INCIDENCE:
                md.update({ifc.DATA_TYPE: ifc.MLOOKED_INC})
            elif v == ifc.COHERENCE and params["coherence_file_paths"]:
                pass
            elif v == ifc.MULTILOOKED and params["coherence_file_paths"]:
                md.update({ifc.DATA_TYPE: ifc.COHERENCE})
            elif v == ifc.MULTILOOKED and not params["coherence_file_paths"]:
                pass
            else:
                raise TypeError("Data Type metadata not recognised")

    # In-memory GDAL driver doesn't support compression so turn it off.
    creation_opts = ["compress=packbits"] if out_driver_type != "MEM" else []
    out_ds = shared.gdal_dataset(
        output_file,
        dst_ds.RasterXSize,
        dst_ds.RasterYSize,
        driver=out_driver_type,
        bands=1,
        dtype=src_dtype,
        metadata=md,
        crs=wkt,
        geotransform=gt,
        creation_opts=creation_opts,
    )

    shared.write_geotiff(resampled_average, out_ds, np.nan)

def get_analysis_extent(ifgs, crop_opt, xlooks, ylooks, user_exts):
    """
    Args:
      ifgs: param crop_opt:
      xlooks: param ylooks:
      user_exts: 
      crop_opt: 
      ylooks: 

    Returns:

    """
    if crop_opt not in CROP_OPTIONS:
        raise ValueError("Unrecognised crop option: %s" % crop_opt)

    if crop_opt == CUSTOM_CROP:
        if not user_exts:
            raise ValueError("No custom cropping extents specified")
        elif len(user_exts) != 4:  # check for required numbers
            raise ValueError("Custom extents must have all 4 values")
        elif len(user_exts) == 4:  # check for non floats
            if not all([isinstance(z, (int, float)) for z in user_exts]):
                raise ValueError("Custom extents must be 4 numbers")

    if not isinstance(xlooks, (int, float)) and isinstance(ylooks, (int, float)):
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise ValueError(msg)

    if not (xlooks > 0 and ylooks > 0):
        msg = "Invalid looks parameter(s), x: %s, y: %s. " "Looks must be an integer greater than zero" % (
        xlooks, ylooks)
        raise ValueError(msg)

    x_step_values = []
    y_step_values = []

    xmin = []
    ymax = []
    xmax = []
    ymin = []

    tfs = []

    x_first = None
    x_last = None
    y_first = None
    y_last = None

    x_step = None
    y_step = None

    for data_path in ifgs:
        dataset = gdal.Open(data_path, gdalconst.GA_ReadOnly)

        x_step = float(dataset.GetGeoTransform()[GDAL_X_CELLSIZE])
        x_step_values.append(x_step)
        y_step = float(dataset.GetGeoTransform()[GDAL_Y_CELLSIZE])
        y_step_values.append(y_step)

        tfs.append(dataset.GetGeoTransform())

        ncols = dataset.RasterXSize
        nrows = dataset.RasterYSize

        x_first = float(dataset.GetGeoTransform()[GDAL_X_FIRST])
        y_first = float(dataset.GetGeoTransform()[GDAL_Y_FIRST])
        x_last = x_first + (x_step * ncols)
        y_last = y_first + (y_step * nrows)

        xmin.append(x_first)
        ymax.append(y_first)
        xmax.append(x_last)
        ymin.append(y_last)

        dataset = None
        del dataset

    if not len(set(x_step_values)) < 2:
        raise ValueError("Grid resolution does not match for x_step: " + str(set(x_step_values)))

    if not len(set(y_step_values)) < 2:
        raise ValueError("Grid resolution does not match for y_step: " + str(set(y_step_values)))

    if crop_opt == MINIMUM_CROP:
        extents = max(xmin), max(ymin), min(xmax), min(ymax)
    elif crop_opt == MAXIMUM_CROP:
        extents = min(xmin), min(ymin), max(xmax), max(ymax)
    elif crop_opt == CUSTOM_CROP:
        xw, ytop, xe, ybot = user_exts
        if ytop < ybot:
            raise ValueError("ERROR Custom crop bounds: ifgyfirst must be greater than ifgylast")
        if xe < xw:
            raise ValueError("ERROR Custom crop bounds: ifgxfirst must be greater than ifgxlast")

        y2 = None
        y1 = None

        for par, crop, orig, step in zip(
                ["x_first", "x_last", "y_first", "y_last"], [xw, xe, ytop, ybot], [x_first, x_last, y_first, y_last],
                [x_step, x_step, y_step, y_step]
        ):
            diff = crop - orig
            nint = round(diff / step)

            if par == "x_first":
                if diff < 0:
                    raise ValueError("Cropped image bounds are outside the original image bounds")
                xmin = orig + (nint * step)

            elif par == "x_last":
                if diff > 0:
                    raise ValueError("Cropped image bounds are outside the original image bounds")
                xmax = orig + (nint * step)

            elif par == "y_first":
                if diff > 0:
                    raise ValueError("Cropped image bounds are outside the original image bounds")
                y1 = orig + (nint * step)

            elif par == "y_last":
                if diff < 0:
                    raise ValueError("Cropped image bounds are outside the original image bounds")
                y2 = orig + (nint * step)
            else:
                raise ValueError("Value error in supplied custom bounds")

        if y2 > y1:
            ymin = y1
            ymax = y2
        else:
            ymin = y2
            ymax = y1

        extents = xmin, ymin, xmax, ymax

        # only need to check crop coords when custom bounds are supplied
        for par, crop, step, param in zip(
                ["x_first", "x_last", "y_first", "y_last"], [xmin, xmax, ymax, ymin], [x_step, x_step, y_step, y_step],
                [x_first, x_last, y_first, y_last]
        ):

            # is diff of the given extent from grid a multiple of X|Y_STEP ?
            diff = abs(crop - param)
            remainder = abs(modf(diff / step)[0])

            # handle cases where division gives remainder near zero, or just < 1
            if (remainder > GRID_TOL) and (remainder < (1 - GRID_TOL)):  # pragma: no cover
                msg = "%s crop extent not within %s of grid coordinate"
                raise ValueError(msg % (par, GRID_TOL))
    else:

        equal = []

        for t in tfs[1:]:
            for i, tf in enumerate(tfs[0]):

                if round(Decimal(tf), 4) == round(Decimal(t[i]), 4):
                    equal.append(True)
                else:
                    equal.append(False)

        if not all(equal):
            msg = "Ifgs do not have the same bounding box for crop option: %s"
            raise Exception(msg % ALREADY_SAME_SIZE)

        xmin, xmax = x_first, x_last
        ymin, ymax = y_first, y_last

        # swap y_first & y_last when using southern hemisphere -ve coords
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        extents = xmin, ymin, xmax, ymax

    return extents
