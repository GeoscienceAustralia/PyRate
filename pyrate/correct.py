#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
# coding: utf-8
"""
This Python module runs the main PyRate correction workflow
"""
import shutil
import os
from pathlib import Path
import pickle as cp
from typing import List
import sys

import pyrate.constants as C
from pyrate.core import (shared, algorithm, mpiops)
from pyrate.core.aps import spatio_temporal_filter
from pyrate.core.covariance import maxvar_vcm_calc_wrapper
from pyrate.core.mst import mst_calc_wrapper
from pyrate.core.orbital import orb_fit_calc_wrapper
from pyrate.core.dem_error import dem_error_calc_wrapper
from pyrate.core.phase_closure.closure_check import iterative_closure_check, mask_pixels_with_unwrapping_errors, \
    update_ifg_list
from pyrate.core.ref_phs_est import ref_phase_est_wrapper
from pyrate.core.refpixel import ref_pixel_calc_wrapper
from pyrate.core.shared import PrereadIfg, Ifg, get_tiles, mpi_vs_multiprocess_logging, join_dicts, \
        nan_and_mm_convert, save_numpy_phase
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import Configuration, MultiplePaths, ConfigException

MAIN_PROCESS = 0


def _create_ifg_dict(params):
    """
    Save the preread_ifgs dict with information about the ifgs that are
    later used for fast loading of Ifg files in IfgPart class

    :param list dest_tifs: List of destination tifs
    :param dict params: Config dictionary
    :param list tiles: List of all Tile instances

    :return: preread_ifgs: Dictionary containing information regarding
                interferograms that are used later in workflow
    :rtype: dict
    """
    dest_tifs = [ifg_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    ifgs_dict = {}
    process_tifs = mpiops.array_split(dest_tifs)
    for d in process_tifs:
        ifg = Ifg(d.tmp_sampled_path) # get the writable copy
        ifg.open()
        nan_and_mm_convert(ifg, params)
        ifgs_dict[d.tmp_sampled_path] = PrereadIfg(
            path=d.sampled_path,
            tmp_path=d.tmp_sampled_path,
            nan_fraction=ifg.nan_fraction,
            first=ifg.first,
            second=ifg.second,
            time_span=ifg.time_span,
            nrows=ifg.nrows,
            ncols=ifg.ncols,
            metadata=ifg.meta_data
        )
        ifg.write_modified_phase() # update phase converted to mm
        ifg.close()
    ifgs_dict = join_dicts(mpiops.comm.allgather(ifgs_dict))

    ifgs_dict = mpiops.run_once(__save_ifgs_dict_with_headers_and_epochs, dest_tifs, ifgs_dict, params, process_tifs)

    params[C.PREREAD_IFGS] = ifgs_dict
    return ifgs_dict


def __save_ifgs_dict_with_headers_and_epochs(dest_tifs, ifgs_dict, params, process_tifs):
    tmpdir = params[C.TMPDIR]
    if not os.path.exists(tmpdir):
        shared.mkdir_p(tmpdir)

    preread_ifgs_file = Configuration.preread_ifgs(params)
    nifgs = len(dest_tifs)
    # add some extra information that's also useful later
    gt, md, wkt = shared.get_geotiff_header_info(process_tifs[0].tmp_sampled_path)
    epochlist = algorithm.get_epochs(ifgs_dict)[0]
    log.info('Found {} unique epochs in the {} interferogram network'.format(len(epochlist.dates), nifgs))
    ifgs_dict['epochlist'] = epochlist
    ifgs_dict['gt'] = gt
    ifgs_dict['md'] = md
    ifgs_dict['wkt'] = wkt
    # dump ifgs_dict file for later use
    cp.dump(ifgs_dict, open(preread_ifgs_file, 'wb'))

    for k in ['gt', 'epochlist', 'md', 'wkt']:
        ifgs_dict.pop(k)

    return ifgs_dict


def _copy_mlooked(params):
    """
    Make a copy of the multi-looked files in the 'tmp_sampled_path'
    for manipulation during correct steps
    """
    log.info("Copying input files into tempdir for manipulation during 'correct' steps")
    mpaths = params[C.INTERFEROGRAM_FILES]
    process_mpaths = mpiops.array_split(mpaths)
    for p in process_mpaths:
        shutil.copy(p.sampled_path, p.tmp_sampled_path)
        Path(p.tmp_sampled_path).chmod(0o664)  # assign write permission as prepifg output is readonly


def main(config):
    """
    Top level function to perform PyRate workflow on given interferograms

    :param dict params: Dictionary of configuration parameters

    :return: refpt: tuple of reference pixel x and y position
    :rtype: tuple
    :return: maxvar: array of maximum variance values of interferograms
    :rtype: ndarray
    :return: vcmt: Variance-covariance matrix array
    :rtype: ndarray
    """
    params = config.__dict__
    mpi_vs_multiprocess_logging("correct", params)

    _copy_mlooked(params)

    return correct_ifgs(config)


def _update_params_with_tiles(params: dict) -> None:
    ifg_path = params[C.INTERFEROGRAM_FILES][0].tmp_sampled_path
    rows, cols = params["rows"], params["cols"]
    tiles = mpiops.run_once(get_tiles, ifg_path, rows, cols)
    # add tiles to params
    params[C.TILES] = tiles


def phase_closure_wrapper(params: dict, config: Configuration) -> dict:
    """
    This wrapper will run the iterative phase closure check to return a stable
    list of checked interferograms, and then mask pixels in interferograms that
    exceed the unwrapping error threshold.
    :param params: Dictionary of PyRate configuration parameters. 
    :param config: Configuration class instance.
    :return: params: Updated dictionary of PyRate configuration parameters.
    """

    if not params[C.PHASE_CLOSURE]:
        log.info("Phase closure correction is not required!")
        return

    rets = iterative_closure_check(config)
    if rets is None:
        log.info("Zero loops are returned from the iterative closure check.")
        log.warning("Abandoning phase closure correction without modifying the interferograms.")
        return
        
    ifg_files, ifgs_breach_count, num_occurences_each_ifg = rets

    # update params with closure checked ifg list
    params[C.INTERFEROGRAM_FILES] = \
        mpiops.run_once(update_ifg_list, ifg_files, params[C.INTERFEROGRAM_FILES])

    if mpiops.rank == 0:
        with open(config.phase_closure_filtered_ifgs_list(params), 'w') as f:
            lines = [p.converted_path + '\n' for p in params[C.INTERFEROGRAM_FILES]]
            f.writelines(lines)

    # mask ifgs with nans where phase unwrap threshold is breached
    if mpiops.rank == 0:
        mask_pixels_with_unwrapping_errors(ifgs_breach_count, num_occurences_each_ifg, params)

    _create_ifg_dict(params) # update the preread_ifgs dict

    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    # update/save the phase_data in the tiled numpy files
    save_numpy_phase(ifg_paths, params)

    return params


correct_steps = {
    'orbfit': orb_fit_calc_wrapper,
    'refphase': ref_phase_est_wrapper,
    'phase_closure': phase_closure_wrapper,
    'demerror': dem_error_calc_wrapper,
    'mst': mst_calc_wrapper,
    'apscorrect': spatio_temporal_filter,
    'maxvar': maxvar_vcm_calc_wrapper,
}


def correct_ifgs(config: Configuration) -> None:
    """
    Top level function to perform PyRate workflow on given interferograms
    """
    params = config.__dict__
    __validate_correct_steps(params)

    # work out the tiling and add to params dict
    _update_params_with_tiles(params)

    # create the preread_ifgs dict for use with tiled data
    _create_ifg_dict(params)

    ifg_paths = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]

    # create initial tiled phase_data numpy files on disc
    save_numpy_phase(ifg_paths, params)

    params[C.REFX_FOUND], params[C.REFY_FOUND] = ref_pixel_calc_wrapper(params)

    # run through the correct steps in user specified sequence
    for step in params['correct']:
        if step == 'phase_closure':
            correct_steps[step](params, config)
        else:
            correct_steps[step](params)
    log.info("Finished 'correct' step")


def __validate_correct_steps(params):
    for step in params['correct']:
        if step not in correct_steps.keys():
            raise ConfigException(f"{step} is not a supported 'correct' step. \n"
                                  f"Supported steps are {correct_steps.keys()}")
