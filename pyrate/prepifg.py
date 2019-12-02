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
This Python script applies optional multilooking and cropping to input
interferogram geotiff files.
"""
# -*- coding: utf-8 -*-
import sys
import logging
import os
from joblib import Parallel, delayed
import numpy as np
from core import shared, mpiops, config as cf, prepifg_helper, gamma, roipac
from core.prepifg_helper import PreprocessError

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
    # Going to assume base_ifg_paths is ordered correcly
    # pylint: disable=too-many-branches
    usage = 'Usage: pyrate prepifg <config_file>'
    if mpiops.size > 1:  # Over-ride input options if this is an MPI job
        params[cf.PARALLEL] = False

    base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST], params[cf.OBS_DIR])

    if params[cf.DEM_FILE] is not None: # optional DEM conversion
        base_ifg_paths.append(params[cf.DEM_FILE])

    processor = params[cf.PROCESSOR]  # roipac or gamma
    if processor == GAMMA: # Incidence/elevation only supported for GAMMA
        if params[cf.APS_INCIDENCE_MAP]:
            base_ifg_paths.append(params[cf.APS_INCIDENCE_MAP])
        if params[cf.APS_ELEVATION_MAP]:
            base_ifg_paths.append(params[cf.APS_ELEVATION_MAP])

    shared.mkdir_p(params[cf.OUT_DIR]) # create output dir

    process_base_ifgs_paths = np.array_split(base_ifg_paths, mpiops.size)[mpiops.rank]
    gtiff_paths = [shared.output_tiff_filename(f, params[cf.OBS_DIR]) for f in process_base_ifgs_paths]
    do_prepifg(gtiff_paths, params)
    log.info("Finished prepifg")


def do_prepifg(gtiff_paths, params):
    """
    Prepare interferograms by applying multilooking/cropping operations.

    :param list gtiff_paths: List of full-res geotiffs
    :param dict params: Parameters dictionary corresponding to config file
    """
    # pylint: disable=expression-not-assigned
    log.info("Preparing interferograms by cropping/multilooking")
    parallel = params[cf.PARALLEL]

    for f in gtiff_paths:
        if not os.path.isfile(f):
            raise Exception("Can not find geotiff: " + str(f) + ". Ensure you have converted your interferograms to geotiffs.")

    ifgs = [prepifg_helper.dem_or_ifg(p) for p in gtiff_paths]
    xlooks, ylooks, crop = cf.transform_params(params)
    user_exts = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST], params[cf.IFG_XLAST], params[cf.IFG_YLAST])
    exts = prepifg_helper.get_analysis_extent(crop, ifgs, xlooks, ylooks, user_exts=user_exts)
    log.debug("Extents (xmin, ymin, xmax, ymax): "+str(exts))
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    if parallel:
        Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(_prepifg_multiprocessing)(p, xlooks, ylooks, exts, thresh, crop, params) for p in gtiff_paths
        )
    else:
        [_prepifg_multiprocessing(p, xlooks, ylooks, exts, thresh, crop, params) for p in gtiff_paths]

def _prepifg_multiprocessing(path, xlooks, ylooks, exts, thresh, crop, params):
    """
    Multiprocessing wrapper for prepifg
    """
    processor = params[cf.PROCESSOR]  # roipac or gamma
    if processor == GAMMA:
        header = gamma.gamma_header(path, params)
    elif processor == ROIPAC:
        header = roipac.roipac_header(path, params)
    else:
        raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
    # If we're performing coherence masking, find the coherence file for this IFG.
    # TODO: Refactor _is_interferogram to be unprotected (remove '_')
    if params[cf.COH_MASK] and shared._is_interferogram(header):
        coherence_path = cf.coherence_paths_for(path, params, tif=True)[0]
        coherence_thresh = params[cf.COH_THRESH]
    else:
        coherence_path = None
        coherence_thresh = None

    prepifg_helper.prepare_ifg(path, xlooks, ylooks, exts, thresh, crop, out_path=params[cf.OUT_DIR], header=header, coherence_path=coherence_path, coherence_thresh=coherence_thresh)

