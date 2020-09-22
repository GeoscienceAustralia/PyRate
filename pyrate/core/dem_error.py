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
This Python module implements the calculation of correction for residual topographic effects (a.k.a. as DEM errors).
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import numpy as np
from typing import Optional, List, Dict, Iterable
from os.path import join

from pyrate.core import geometry, shared, ifgconstants as ifc, mpiops, config as cf
from pyrate.core.logger import pyratelogger as log
from pyrate.core.shared import Ifg, Geometry


def remove_dem_error(ifgs: List, params: dict) -> None:
    """
    Wrapper function for PyRate DEM error removal functionality.
    """
    # read first IFG
    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    ifg0_path = ifg_paths[0]
    ifg0 = Ifg(ifg0_path)
    ifg0.open(readonly=True)

    # not currently implemented for ROIPAC data which breaks some tests
    # if statement can be deleted once ROIPAC is deprecated from PyRate
    if not ifg0.meta_data[ifc.PYRATE_INSAR_PROCESSOR] == 'ROIPAC':

        # read radar azimuth and range tif files
        rdc_az_file = join(params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        az = geom_az.geometry_data
        rdc_rg_file = join(params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        rg = geom_rg.geometry_data
        # read lon and lat values of multi-looked ifg (first ifg only)
        lon, lat = geometry.get_lonlat_coords(ifg0)

        # loop over all Ifgs in stack
        ifgs = [shared.Ifg(p) for p in ifg_paths] if isinstance(ifgs[0], str) else ifgs
        process_ifgs = mpiops.array_split(ifgs)
        for ifg in process_ifgs:
            ifg.open(readonly=True)
            # calculate look angle for interferograms (using the Near Range of the primary SLC)
            look_angle = geometry.calc_local_geometry(ifg, None, rg, lon, lat, params)
            bperp = geometry.calc_local_baseline(ifg, az, look_angle, params)
    # todo: add code to invert for DEM error using the per-pixel baseline values saved in bperp and the unwrapped phase


def dem_error_calc_wrapper(params: dict) -> None:
    """
    MPI wrapper for DEM error correction
    """
    if params[cf.BASE_FILE_LIST] is None:
        log.info('No baseline files supplied: DEM error correction not computed')
        return

    if params[cf.LT_FILE] is None:
        log.info('No lookup table file supplied: DEM error correction not computed')
        return

    log.info('Calculating DEM error correction')
    multi_paths = params[cf.INTERFEROGRAM_FILES]
    ifg_paths = [p.tmp_sampled_path for p in multi_paths]
    remove_dem_error(ifg_paths, params)
    mpiops.comm.barrier()
    shared.save_numpy_phase(ifg_paths, params)
    log.debug('Finished DEM error correction')


