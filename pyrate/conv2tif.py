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
This Python script converts ROI_PAC or GAMMA format input interferograms 
into geotiff format files
"""
# -*- coding: utf-8 -*-
import sys
import os
import logging
from joblib import Parallel, delayed
import numpy as np

from core.prepifg_helper import PreprocessError
from core import shared, mpiops, config as cf, gamma, roipac

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

GAMMA = 1
ROIPAC = 0


def main(params=None):
    """
    Parse parameters and prepare files for conversion.

    :param dict params: Parameters dictionary read in from the config file
    """
    # TODO: looks like base_ifg_paths are ordered according to ifg list
    # This probably won't be a problem because input list won't be reordered
    # and the original gamma generated list is ordered) this may not affect
    # the important pyrate stuff anyway, but might affect gen_thumbs.py.
    # Going to assume base_ifg_paths is ordered correcly
    # pylint: disable=too-many-branches

    base_ifg_paths = []
    for interferogram_file in params["interferogram_files"]:
        base_ifg_paths.append(interferogram_file.unwrapped_path)

    if params[cf.COH_MASK]:
        base_ifg_paths.extend(cf.coherence_paths(params))

    # optional DEM conversion
    if params[cf.DEM_FILE] is not None:
        base_ifg_paths.append(params[cf.DEM_FILE])

    processor = params[cf.PROCESSOR]  # roipac or gamma
    if processor == GAMMA:  # Incidence/elevation only supported for GAMMA
        if params[cf.APS_INCIDENCE_MAP]:
            base_ifg_paths.append(params[cf.APS_INCIDENCE_MAP])
        if params[cf.APS_ELEVATION_MAP]:
            base_ifg_paths.append(params[cf.APS_ELEVATION_MAP])

    process_ifgs_paths = np.array_split(base_ifg_paths, mpiops.size)[mpiops.rank]

    log.info("Converting input interferograms to geotiff")

    log.info("Running geotiff conversion in serial")
    for b in process_ifgs_paths:
        _geotiff_multiprocessing(b, params)

    log.info("Finished conv2tif")


def _geotiff_multiprocessing(unw_path, params):
    """
    Multiprocessing wrapper for full-res geotiff conversion
    """
    # TODO: Need a more robust method for identifying coherence files.
    if params[cf.COH_FILE_DIR] and unw_path.endswith('.cc'):
        # If the user has provided a dir for coherence files, place 
        #  converted coherence files in that directory.
        dest = shared.output_tiff_filename(unw_path, params[cf.COH_FILE_DIR])
    else:
        dest = shared.output_tiff_filename(unw_path, params[cf.OBS_DIR])
    processor = params[cf.PROCESSOR]  # roipac or gamma

    # Create full-res geotiff if not already on disk
    if not os.path.exists(dest):
        if processor == GAMMA:
            header = gamma.gamma_header(unw_path, params)
        elif processor == ROIPAC:
            header = roipac.roipac_header(unw_path, params)
        else:
            raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
        shared.write_fullres_geotiff(header, unw_path, dest, nodata=params[cf.NO_DATA_VALUE])
        return dest
    else:
        log.info("Full-res geotiff already exists")
        return None

