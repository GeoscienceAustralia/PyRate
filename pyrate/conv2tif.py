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
import os
from typing import Tuple, List
from joblib import Parallel, delayed
import numpy as np

from pyrate.core.prepifg_helper import PreprocessError
from pyrate.core import shared, mpiops, config as cf, gamma, roipac
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import MultiplePaths
from pyrate.core.shared import mpi_vs_multiprocess_logging

GAMMA = 1
ROIPAC = 0
GEOTIF = 2


def main(params):
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

    if params[cf.PROCESSOR] == 2:  # if geotif
        log.warning("conv2tif not required for geotiff!")
        return

    mpi_vs_multiprocess_logging("conv2tif", params)

    base_ifg_paths = params[cf.INTERFEROGRAM_FILES]

    if params[cf.COH_MASK]:
        base_ifg_paths.extend(params[cf.COHERENCE_FILE_PATHS])

    if params[cf.DEM_FILE] is not None:  # optional DEM conversion
        base_ifg_paths.append(params[cf.DEM_FILE_PATH])

    process_base_ifgs_paths = np.array_split(base_ifg_paths, mpiops.size)[mpiops.rank]
    gtiff_paths = do_geotiff(process_base_ifgs_paths, params)
    mpiops.comm.barrier()
    log.info("Finished conv2tif")
    return gtiff_paths


def do_geotiff(unw_paths: List[MultiplePaths], params: dict) -> List[str]:
    """
    Convert input interferograms to geotiff format.
    """
    # pylint: disable=expression-not-assigned
    log.info("Converting input interferograms to geotiff")
    parallel = params[cf.PARALLEL]

    if parallel:
        log.info("Running geotiff conversion in parallel with {} processes".format(params[cf.PROCESSES]))
        dest_base_ifgs = Parallel(n_jobs=params[cf.PROCESSES], verbose=shared.joblib_log_level(cf.LOG_LEVEL))(
            delayed(_geotiff_multiprocessing)(p, params) for p in unw_paths)
    else:
        log.info("Running geotiff conversion in serial")
        dest_base_ifgs = [_geotiff_multiprocessing(b, params) for b in unw_paths]
    return dest_base_ifgs


def _geotiff_multiprocessing(unw_path: MultiplePaths, params: dict) -> Tuple[str, bool]:
    """
    Multiprocessing wrapper for full-res geotiff conversion
    """
    # TODO: Need a more robust method for identifying coherence files.
    dest = unw_path.converted_path
    processor = params[cf.PROCESSOR]  # roipac or gamma

    # Create full-res geotiff if not already on disk
    if not os.path.exists(dest):
        if processor == GAMMA:
            header = gamma.gamma_header(unw_path.unwrapped_path, params)
        elif processor == ROIPAC:
            log.info("Warning: ROI_PAC support will be deprecated in a future PyRate release")
            header = roipac.roipac_header(unw_path.unwrapped_path, params)
        else:
            raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
        shared.write_fullres_geotiff(header, unw_path.unwrapped_path, dest, nodata=params[cf.NO_DATA_VALUE])
        return dest, True
    else:
        log.warning(f"Full-res geotiff already exists in {dest}! Returning existing geotiff!")
        return dest, False
