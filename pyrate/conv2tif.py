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

from core import shared, config as cf, gamma, roipac
from core.logger import pyratelogger as log
from core.mpiops import rank, comm, size, chunks
from core.prepifg_helper import PreprocessError

GAMMA = 1
ROIPAC = 0


def main(params):
    """
    Args:
        params:
    """
    # TODO: looks like base_ifg_paths are ordered according to ifg list
    # This probably won't be a problem because input list won't be reordered
    # and the original gamma generated list is ordered) this may not affect
    # the important pyrate stuff anyway, but might affect gen_thumbs.py.
    # Going to assume base_ifg_paths is ordered correcly
    if rank == 0:
        log.info("Collecting jobs: conv2tif")
        # a job is list of parameters passed to multiprocessing function
        jobs = []
        for interferogram_file in params["interferogram_files"]:
            jobs.append((interferogram_file.unwrapped_path, interferogram_file.converted_path, params))

        if params[cf.COH_MASK]:
            for coherence_file in params["coherence_file_paths"]:
                jobs.append((coherence_file.unwrapped_path, coherence_file.converted_path, params))

        # optional DEM conversion
        if params[cf.DEM_FILE] is not None:
            jobs.append((params["dem_file"].unwrapped_path, params["dem_file"].converted_path, params))
        jobs = chunks(jobs, size)
    else:
        jobs = None

    jobs = comm.scatter(jobs, root=0)

    for job in jobs:
        _geotiff_multiprocessing(*job)


def _geotiff_multiprocessing(input_file_name, output_file_name, params):
    """Multiprocessing wrapper for full-res geotiff conversion :param
    input_file_name: :type input_file_name: :param output_file_name: :type
    output_file_name: :param params: :type params:

    Args:
        input_file_name:
        output_file_name:
        params:
    """
    processor = params[cf.PROCESSOR]  # roipac or gamma
    if not os.path.exists(output_file_name):
        log.info("Started processing: " + str(input_file_name))
        if processor == GAMMA:
            header = gamma.gamma_header(input_file_name, params)
        elif processor == ROIPAC:
            header = roipac.roipac_header(input_file_name, params)
        else:
            raise PreprocessError("Processor must be ROI_PAC (0) or GAMMA (1)")
        shared.write_fullres_geotiff(header, input_file_name, output_file_name, nodata=params[cf.NO_DATA_VALUE])
        log.info("Finished creating: " + str(output_file_name))
    else:
        log.info("Full resolution GeoTIFF already exists: " + str(output_file_name))
