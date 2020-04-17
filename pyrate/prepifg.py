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
import os
from typing import List
from joblib import Parallel, delayed
import numpy as np
from pyrate.core import shared, mpiops, config as cf, prepifg_helper, gamma, roipac
from pyrate.core.prepifg_helper import PreprocessError
from pyrate.core.logger import pyratelogger as log


GAMMA = 1
ROIPAC = 0
GEOTIF = 2


def main(params):
    """
    Main workflow function for preparing interferograms for PyRate.

    :param dict params: Parameters dictionary read in from the config file
    """
    # TODO: looks like ifg_paths are ordered according to ifg list
    # This probably won't be a problem because input list won't be reordered
    # and the original gamma generated list is ordered) this may not affect
    # the important pyrate stuff anyway, but might affect gen_thumbs.py.
    # Going to assume ifg_paths is ordered correcly
    # pylint: disable=too-many-branches

    if mpiops.size > 1:  # Over-ride input options if this is an MPI job
        log.debug("Running using mpi processing. Disabling parallel processing.")
        params[cf.PARALLEL] = False

    ifg_paths = params[cf.INTERFEROGRAM_FILES]
    if params[cf.DEM_FILE] is not None:  # optional DEM conversion
        ifg_paths.append(params[cf.DEM_FILE_PATH])

    shared.mkdir_p(params[cf.OUT_DIR])  # create output dir

    process_base_ifgs_paths = np.array_split(ifg_paths, mpiops.size)[mpiops.rank]
    gtiff_paths = [p.converted_path for p in process_base_ifgs_paths]
    do_prepifg(gtiff_paths, params)
    mpiops.comm.barrier()
    log.info("Finished prepifg")


def do_prepifg(gtiff_paths: List[str], params: dict) -> None:
    """
    Prepare interferograms by applying multilooking/cropping operations.

    :param list gtiff_paths: List of full-res geotiffs
    :param dict params: Parameters dictionary corresponding to config file
    """
    # pylint: disable=expression-not-assigned
    parallel = params[cf.PARALLEL]

    for f in gtiff_paths:
        if not os.path.isfile(f):
            raise Exception("Can not find geotiff: " + str(f) +
                            ". Ensure you have converted your interferograms to geotiffs.")

    ifgs = [prepifg_helper.dem_or_ifg(p) for p in gtiff_paths]
    xlooks, ylooks, crop = cf.transform_params(params)
    user_exts = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST], params[cf.IFG_XLAST], params[cf.IFG_YLAST])
    exts = prepifg_helper.get_analysis_extent(crop, ifgs, xlooks, ylooks, user_exts=user_exts)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    if parallel:
        Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(_prepifg_multiprocessing)(p, xlooks, ylooks, exts, thresh, crop, params) for p in gtiff_paths
        )
    else:
        for gtiff_path in gtiff_paths:
            _prepifg_multiprocessing(gtiff_path, xlooks, ylooks, exts, thresh, crop, params)


def _prepifg_multiprocessing(path, xlooks, ylooks, exts, thresh, crop, params):
    """
    Multiprocessing wrapper for prepifg
    """
    processor = params[cf.PROCESSOR]  # roipac, gamma or geotif
    if (processor == GAMMA) or (processor == GEOTIF):
        header = gamma.gamma_header(path, params)
    elif processor == ROIPAC:
        header = roipac.roipac_header(path, params)
    else:
        raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')

    # If we're performing coherence masking, find the coherence file for this IFG.
    if params[cf.COH_MASK] and shared._is_interferogram(header):
        coherence_path = cf.coherence_paths_for(path, params, tif=True)[0]
        coherence_thresh = params[cf.COH_THRESH]
    else:
        coherence_path = None
        coherence_thresh = None

    prepifg_helper.prepare_ifg(path, xlooks, ylooks, exts, thresh, crop, out_path=params[cf.OUT_DIR], header=header,
                               coherence_path=coherence_path, coherence_thresh=coherence_thresh)

