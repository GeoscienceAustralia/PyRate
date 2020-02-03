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
import logging

from core.prepifg_helper import PreprocessError
from core import shared, config as cf, gamma, roipac

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

    for interferogram_file in params["interferogram_files"]:
        _geotiff_multiprocessing(interferogram_file.unwrapped_path, interferogram_file.converted_path, params)

    if params[cf.COH_MASK]:
        for coherence_file in params["coherence_file_paths"]:
            _geotiff_multiprocessing(coherence_file.unwrapped_path, coherence_file.converted_path, params)

    # optional DEM conversion
    if params[cf.DEM_FILE] is not None:
        _geotiff_multiprocessing(params["dem_file"].unwrapped_path, params["dem_file"].converted_path, params)

    log.info("Finished conv2tif")


def _geotiff_multiprocessing(input_file_name, output_file_name, params):
    """
    Multiprocessing wrapper for full-res geotiff conversion
    """
    # TODO: Need a more robust method for identifying coherence files.
    if params[cf.COH_FILE_DIR] and input_file_name.endswith('.cc'):
        # If the user has provided a dir for coherence files, place
        # converted coherence files in that directory.
        output_file_name = shared.output_tiff_filename(input_file_name, params[cf.COH_FILE_DIR])
    else:
        output_file_name = shared.output_tiff_filename(input_file_name, params[cf.OBS_DIR])
    processor = params[cf.PROCESSOR]  # roipac or gamma

    # Create full-res geotiff if not already on disk
    if not os.path.exists(output_file_name):
        if processor == GAMMA:
            header = gamma.gamma_header(input_file_name, params)
        elif processor == ROIPAC:
            header = roipac.roipac_header(input_file_name, params)
        else:
            raise PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
        shared.write_fullres_geotiff(header, input_file_name, output_file_name, nodata=params[cf.NO_DATA_VALUE])
        return output_file_name
    else:
        log.info("Full-res geotiff already exists")
        return None
