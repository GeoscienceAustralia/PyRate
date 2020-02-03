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
import logging
import os

from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from osgeo import gdalconst
from osgeo import gdal_array
from gdalconst import GA_Update, GA_ReadOnly
from decimal import Decimal
from math import modf
import numpy as np
from core import shared, mpiops, config as cf, prepifg_helper, gamma, roipac
from constants import CROP_OPTIONS, CUSTOM_CROP, GRID_TOL, GDAL_X_CELLSIZE, GDAL_Y_CELLSIZE, GDAL_X_FIRST, GDAL_Y_FIRST, MINIMUM_CROP, MAXIMUM_CROP, CUSTOM_CROP, ALREADY_SAME_SIZE
from core.prepifg_helper import PreprocessError

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

GAMMA = 1
ROIPAC = 0


def main(params):
    """
    Main workflow function for preparing interferograms for PyRate.

    :param dict params: Parameters dictionary read in from the config file
    """
    # TODO: looks like base_ifg_paths are ordered according to ifg list
    # This probably won't be a problem because input list won't be reordered
    # and the original gamma generated list is ordered) this may not affect
    # the important pyrate stuff anyway, but might affect gen_thumbs.py.
    # Going to assume base_ifg_paths is ordered correctly

    xlooks, ylooks, crop = cf.transform_params(params)
    user_extents = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST], params[cf.IFG_XLAST], params[cf.IFG_YLAST])
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]

    datasets_to_calcualte_extents = []
    for interferogram_file in params["interferogram_files"]:
        if not os.path.isfile(interferogram_file.converted_path):
            raise Exception("Can not find geotiff: " + str(interferogram_file.converted_path) + ". Ensure you have converted your interferograms to geotiffs.")
        datasets_to_calcualte_extents.append(interferogram_file.converted_path)

    # optional DEM conversion
    if params["dem_file"] is not None:
        if not os.path.isfile(params["dem_file"].converted_path):
            raise Exception("Can not find geotiff: " + str(params["dem_file"].converted_path) + ". Ensure you have converted your interferograms to geotiffs.")
        datasets_to_calcualte_extents.append(params["dem_file"].converted_path)

    extents = get_analysis_extent(datasets_to_calcualte_extents, crop, xlooks, ylooks, user_exts=user_extents)


    # processor = params[cf.PROCESSOR]  # roipac or gamma
    # if processor == GAMMA: # Incidence/elevation only supported for GAMMA
    #     if params[cf.APS_INCIDENCE_MAP]:
    #         base_ifg_paths.append(params[cf.APS_INCIDENCE_MAP])
    #     if params[cf.APS_ELEVATION_MAP]:
    #         base_ifg_paths.append(params[cf.APS_ELEVATION_MAP])

    log.info("Preparing interferograms by cropping/multilooking")
    for interferogram_file in params["interferogram_files"]:
        _prepifg_multiprocessing(interferogram_file.converted_path, interferogram_file.sampled_path, xlooks, ylooks, extents, thresh, crop, params)

    # optional DEM conversion
    if params["dem_file"] is not None:
        _prepifg_multiprocessing(params["dem_file"].converted_path, params["dem_file"].sampled_path, xlooks, ylooks, extents, thresh, crop, params)

def _prepifg_multiprocessing(input_path, output_path, xlooks, ylooks, extents, thresh, crop_opt, params):
    """
    Multiprocessing wrapper for prepifg
    """
    processor = params[cf.PROCESSOR]  # roipac or gamma
    if processor == GAMMA:
        header = gamma.gamma_header(input_path, params)
    elif processor == ROIPAC:
        header = roipac.roipac_header(input_path, params)
    else:
        raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
    # If we're performing coherence masking, find the coherence file for this IFG.
    # TODO: Refactor _is_interferogram to be unprotected (remove '_')
    if params[cf.COH_MASK] and shared._is_interferogram(header):
        coherence_path = cf.coherence_paths_for(input_path, params, tif=True)[0]
        coherence_thresh = params[cf.COH_THRESH]
    else:
        coherence_path = None
        coherence_thresh = None

    prepifg_helper.prepare_ifg(input_path, output_path, xlooks, ylooks, extents, thresh, crop_opt, header, coherence_path, coherence_thresh)


def get_analysis_extent(ifgs, crop_opt, xlooks, ylooks, user_exts):

    if crop_opt not in CROP_OPTIONS:
        raise ValueError("Unrecognised crop option: %s" % crop_opt)

    if crop_opt == CUSTOM_CROP:
        if not user_exts:
            raise ValueError('No custom cropping extents specified')
        elif len(user_exts) != 4:  # check for required numbers
            raise ValueError('Custom extents must have all 4 values')
        elif len(user_exts) == 4:  # check for non floats
            if not all([isinstance(z, (int, float)) for z in user_exts]):
                raise ValueError('Custom extents must be 4 numbers')

    if not isinstance(xlooks, (int, float)) and isinstance(ylooks, (int, float)):
        msg = "Non-numeric looks parameter(s), x: %s, y: %s" % (xlooks, ylooks)
        raise ValueError(msg)

    if not (xlooks > 0 and ylooks > 0):
        msg = "Invalid looks parameter(s), x: %s, y: %s. " \
              "Looks must be an integer greater than zero" % (xlooks, ylooks)
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
        dataset = gdal.Open(data_path, GA_ReadOnly)

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
        raise ValueError("Grid resolution does not match for x_step: "+str(set(x_step_values)))

    if not len(set(y_step_values)) < 2:
        raise ValueError("Grid resolution does not match for y_step: "+str(set(y_step_values)))

    if crop_opt == MINIMUM_CROP:
        extents = max(xmin), max(ymin),  min(xmax), min(ymax)
    elif crop_opt == MAXIMUM_CROP:
        extents = min(xmin), min(ymin), max(xmax), max(ymax)
    elif crop_opt == CUSTOM_CROP:
        xw, ytop, xe, ybot = user_exts
        if ytop < ybot:
            raise ValueError('ERROR Custom crop bounds: ifgyfirst must be greater than ifgylast')
        if xe < xw:
            raise ValueError('ERROR Custom crop bounds: ifgxfirst must be greater than ifgxlast')

        y2 = None
        y1 = None

        for par, crop, orig, step in zip(['x_first', 'x_last', 'y_first', 'y_last'],
                                         [xw, xe, ytop, ybot],
                                         [x_first, x_last, y_first, y_last],
                                         [x_step, x_step, y_step, y_step]):
            diff = crop - orig
            nint = round(diff / step)

            if par == 'x_first':
                if diff < 0:
                    raise ValueError('Cropped image bounds are outside the original image bounds')
                xmin = orig + (nint * step)

            elif par == 'x_last':
                if diff > 0:
                    raise ValueError('Cropped image bounds are outside the original image bounds')
                xmax = orig + (nint * step)

            elif par == 'y_first':
                if diff > 0:
                    raise ValueError('Cropped image bounds are outside the original image bounds')
                y1 = orig + (nint * step)

            elif par == 'y_last':
                if diff < 0:
                    raise ValueError('Cropped image bounds are outside the original image bounds')
                y2 = orig + (nint * step)
            else:
                raise ValueError('Value error in supplied custom bounds')

        if y2 > y1:
            ymin = y1
            ymax = y2
        else:
            ymin = y2
            ymax = y1

        extents = xmin, ymin, xmax, ymax

        # only need to check crop coords when custom bounds are supplied
        for par, crop, step in zip(['x_first', 'x_last', 'y_first', 'y_last'],
                                   [xmin, xmax, ymax, ymin],
                                   [x_step, x_step, y_step, y_step]):

            # is diff of the given extent from grid a multiple of X|Y_STEP ?
            param = x_first, x_last, y_first, y_last
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
            msg = 'Ifgs do not have the same bounding box for crop option: %s'
            raise PreprocessError(msg % ALREADY_SAME_SIZE)

        xmin, xmax = x_first, x_last
        ymin, ymax = y_first, y_last

        # swap y_first & y_last when using southern hemisphere -ve coords
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        extents = xmin, ymin, xmax, ymax

    return extents