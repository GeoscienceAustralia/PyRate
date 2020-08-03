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
# coding: utf-8
"""
This Python module runs the main PyRate processing workflow
"""
import os
import shutil
from os.path import join
from pathlib import Path
import pickle as cp
from pyrate.core import (shared, algorithm, mpiops, config as cf)
from pyrate.core.config import ConfigException
from pyrate.core.aps import wrap_spatio_temporal_filter
from pyrate.core.covariance import maxvar_vcm_calc_wrapper
from pyrate.core.mst import mst_calc_wrapper
from pyrate.core.orbital import orb_fit_calc_wrapper
from pyrate.core.ref_phs_est import ref_phase_est_wrapper
from pyrate.core.refpixel import ref_pixel_calc_wrapper
from pyrate.core.shared import PrereadIfg, get_tiles, mpi_vs_multiprocess_logging
from pyrate.core.logger import pyratelogger as log
from pyrate.core.stack import stack_calc_wrapper
from pyrate.core.timeseries import timeseries_calc_wrapper


def _join_dicts(dicts):
    """
    Function to concatenate dictionaries
    """
    if dicts is None:  # pragma: no cover
        return
    assembled_dict = {k: v for D in dicts for k, v in D.items()}
    return assembled_dict


def _create_ifg_dict(params):
    """
    1. Convert ifg phase data into numpy binary files.
    2. Save the preread_ifgs dict with information about the ifgs that are
    later used for fast loading of Ifg files in IfgPart class

    :param list dest_tifs: List of destination tifs
    :param dict params: Config dictionary
    :param list tiles: List of all Tile instances

    :return: preread_ifgs: Dictionary containing information regarding
                interferograms that are used later in workflow
    :rtype: dict
    """
    dest_tifs = [ifg_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    ifgs_dict = {}
    process_tifs = mpiops.array_split(dest_tifs)
    for d in process_tifs:
        ifg = shared._prep_ifg(d.sampled_path, params)
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
        ifg.close()
    ifgs_dict = _join_dicts(mpiops.comm.allgather(ifgs_dict))

    ifgs_dict = mpiops.run_once(__save_ifgs_dict_with_headers_and_epochs, dest_tifs, ifgs_dict, params, process_tifs)

    params[cf.PREREAD_IFGS] = ifgs_dict
    log.debug('Finished converting phase_data to numpy in process {}'.format(mpiops.rank))
    return ifgs_dict


def __save_ifgs_dict_with_headers_and_epochs(dest_tifs, ifgs_dict, params, process_tifs):
    tmpdir = params[cf.TMPDIR]
    if not os.path.exists(tmpdir):
        shared.mkdir_p(tmpdir)

    preread_ifgs_file = join(params[cf.TMPDIR], 'preread_ifgs.pk')
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
    log.info("Copying input files into tempdir for manipulation during process steps")
    mpaths = params[cf.INTERFEROGRAM_FILES]
    process_mpaths = mpiops.array_split(mpaths)
    for p in process_mpaths:
        shutil.copy(p.sampled_path, p.tmp_sampled_path)
        Path(p.tmp_sampled_path).chmod(0o664)  # assign write permission as prepifg output is readonly


def main(params):
    """
    Top level function to perform PyRate workflow on given interferograms

    :return: refpt: tuple of reference pixel x and y position
    :rtype: tuple
    :return: maxvar: array of maximum variance values of interferograms
    :rtype: ndarray
    :return: vcmt: Variance-covariance matrix array
    :rtype: ndarray
    """
    mpi_vs_multiprocess_logging("process", params)

    if mpiops.size > 1:  # turn of multiprocessing during mpi jobs
        params[cf.PARALLEL] = False

    # Make a copy of the multi-looked files for manipulation during process steps
    _copy_mlooked(params)

    return process_ifgs(params)


def _update_params_with_tiles(params: dict) -> None:
    ifg_path = params[cf.INTERFEROGRAM_FILES][0].sampled_path
    rows, cols = params["rows"], params["cols"]
    tiles = mpiops.run_once(get_tiles, ifg_path, rows, cols)
    # add tiles to params
    params[cf.TILES] = tiles


process_steps = {
    'orbfit': orb_fit_calc_wrapper,
    'refphase': ref_phase_est_wrapper,
    'mst': mst_calc_wrapper,
    'apscorrect': wrap_spatio_temporal_filter,
    'maxvar': maxvar_vcm_calc_wrapper,
    'timeseries': timeseries_calc_wrapper,
    'stack': stack_calc_wrapper
}


def process_ifgs(params: dict) -> None:
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

    __validate_process_steps(params)

    # house keeping
    _update_params_with_tiles(params)
    _create_ifg_dict(params)
    ref_pixel_calc_wrapper(params)

    # run through the process steps in user specified sequence
    for step in params['process']:
        process_steps[step](params)

    log.info('Finished process workflow steps')


def __validate_process_steps(params):
    for step in params['process']:
        if step not in process_steps.keys():
            raise ConfigException(f"{step} is not a supported process step. \n"
                                  f"Supported steps are {process_steps.keys()}")
