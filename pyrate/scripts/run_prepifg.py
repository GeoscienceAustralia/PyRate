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
This Python script converts ROI_PAC or GAMMA format unwrapped
interferograms into geotiffs and applies optional multilooking and cropping.
"""
# -*- coding: utf-8 -*-
import sys
import logging
import os
from os.path import join
import re
import glob2
from joblib import Parallel, delayed
import numpy as np
from pyrate import prepifg
from pyrate import config as cf
from pyrate import roipac
from pyrate import gamma
from pyrate.shared import write_geotiff, mkdir_p, output_tiff_filename
import pyrate.ifgconstants as ifc
from pyrate import mpiops

log = logging.getLogger(__name__)

ROI_PAC_HEADER_FILE_EXT = 'rsc'
GAMMA = 1
ROIPAC = 0


def main(params=None):
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

    if params:
        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])
    else:
        # if params not provided read from config file
        if (not params) and (len(sys.argv) < 3):
            print(usage)
            return
        base_ifg_paths, _, params = cf.get_ifg_paths(sys.argv[2])

    base_ifg_paths.append(params[cf.DEM_FILE])
    processor = params[cf.PROCESSOR]  # roipac or gamma
    if processor == GAMMA: # Incidence/elevation only supported for GAMMA
        if params[cf.APS_INCIDENCE_MAP]:
            base_ifg_paths.append(params[cf.APS_INCIDENCE_MAP])
        if params[cf.APS_ELEVATION_MAP]:
            base_ifg_paths.append(params[cf.APS_ELEVATION_MAP])

    mkdir_p(params[cf.OUT_DIR]) # create output dir

    process_base_ifgs_paths = np.array_split(base_ifg_paths, mpiops.size)[mpiops.rank]
    if processor == ROIPAC:
        roipac_prepifg(process_base_ifgs_paths, params)
    elif processor == GAMMA:
        gamma_prepifg(process_base_ifgs_paths, params)
    else:
        raise prepifg.PreprocessError('Processor must be ROI_PAC (0) or GAMMA (1)')
    log.info("Finished prepifg")


PTN = re.compile(r'\d{8}')  # match 8 digits for the dates

def get_header_paths(input_file, slc_dir=None):
    """
    Function that matches input GAMMA file names with GAMMA header file names

    :param str input_file: input GAMMA .unw file.
    :param str slc_dir: GAMMA SLC header file directory

    :return: list of matching header files
    :rtype: list
    """
    if slc_dir:
        dir_name = slc_dir
        _, file_name = os.path.split(input_file)
    else:  # header file must exist in the same dir as that of .unw
        dir_name, file_name = os.path.split(input_file)
    matches = PTN.findall(file_name)
    return [glob2.glob(join(dir_name, '**/*%s*slc.par' % m))[0] for m in matches]

def roipac_prepifg(base_ifg_paths, params):
    """
    Prepare ROI_PAC interferograms which combines both conversion to geotiff
    and multilooking/cropping operations.

    :param list base_ifg_paths: List of unwrapped interferograms
    :param dict params: Parameters dictionary corresponding to config file
    """
    log.info("Preparing ROI_PAC format interferograms")
    parallel = params[cf.PARALLEL]

    if parallel:
        log.info("Parallel prepifg is not implemented for ROI_PAC")

    log.info("Running prepifg in serial")
    xlooks, ylooks, crop = cf.transform_params(params)
    rsc_file = os.path.join(params[cf.DEM_HEADER_FILE])
    if rsc_file is not None:
        projection = roipac.parse_header(rsc_file)[ifc.PYRATE_DATUM]
    else:
        raise roipac.RoipacException('No DEM resource/header file is '
                                     'provided')
    dest_base_ifgs = [os.path.join(params[cf.OUT_DIR],
                                   os.path.basename(q).split('.')[0] + '_' +
                                   os.path.basename(q).split('.')[1] + '.tif')
                      for q in base_ifg_paths]

    for b, d in zip(base_ifg_paths, dest_base_ifgs):
        header_file = "%s.%s" % (b, ROI_PAC_HEADER_FILE_EXT)
        header = roipac.manage_header(header_file, projection)
        write_geotiff(header, b, d, nodata=params[cf.NO_DATA_VALUE])
    prepifg.prepare_ifgs(
        dest_base_ifgs, crop_opt=crop, xlooks=xlooks, ylooks=ylooks)


def gamma_prepifg(base_unw_paths, params):
    """
    Prepare GAMMA interferograms which combines both conversion to geotiff
    and multilooking/cropping operations.

    :param list base_unw_paths: List of unwrapped interferograms
    :param dict params: Parameters dictionary corresponding to config file
    """
    # pylint: disable=expression-not-assigned
    log.info("Preparing GAMMA format interferograms")
    parallel = params[cf.PARALLEL]

    # dest_base_ifgs: location of geo_tif's
    if parallel:
        log.info("Running prepifg in parallel with {} "
                 "processes".format(params[cf.PROCESSES]))
        dest_base_ifgs = Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(_gamma_multiprocessing)(p, params)
            for p in base_unw_paths)
    else:
        log.info("Running prepifg in serial")
        dest_base_ifgs = [_gamma_multiprocessing(b, params)
                          for b in base_unw_paths]

    ifgs = [prepifg.dem_or_ifg(p) for p in dest_base_ifgs]
    xlooks, ylooks, crop = cf.transform_params(params)
    user_exts = (params[cf.IFG_XFIRST], params[cf.IFG_YFIRST],
                 params[cf.IFG_XLAST], params[cf.IFG_YLAST])
    exts = prepifg.get_analysis_extent(crop, ifgs, xlooks, ylooks,
                                       user_exts=user_exts)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    if parallel:
        Parallel(n_jobs=params[cf.PROCESSES], verbose=50)(
            delayed(prepifg.prepare_ifg)(p, xlooks, ylooks, exts, thresh, crop)
            for p in dest_base_ifgs)
    else:
        [prepifg.prepare_ifg(i, xlooks, ylooks, exts, thresh, crop) for i in dest_base_ifgs]


def _gamma_multiprocessing(unw_path, params):
    """
    Multiprocessing wrapper for GAMMA full-res geotiff conversion
    """
    dem_hdr_path = params[cf.DEM_HEADER_FILE]
    slc_dir = params[cf.SLC_DIR]
    header_paths = get_header_paths(unw_path, slc_dir=slc_dir)
    combined_headers = gamma.manage_headers(dem_hdr_path, header_paths)
    dest = output_tiff_filename(unw_path, params[cf.OUT_DIR])

    if os.path.basename(unw_path).split('.')[1] == (params[cf.APS_INCIDENCE_EXT] or params[cf.APS_ELEVATION_EXT]):
        # TODO: implement incidence class here
        combined_headers['FILE_TYPE'] = 'Incidence'

    # Create full-res geotiff if not already on disk
    if not os.path.exists(dest):
        write_geotiff(combined_headers, unw_path, dest,
                      nodata=params[cf.NO_DATA_VALUE])
    else:
        log.info("Full-res geotiff already exists")

    return dest
