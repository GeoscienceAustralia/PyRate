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
from os.path import join
import pickle as cp
from typing import List, Tuple
import numpy as np

from pyrate.core import (shared, algorithm, orbital, ref_phs_est as rpe, 
                         ifgconstants as ifc, mpiops, config as cf, 
                         timeseries, mst, covariance as vcm_module, 
                         stack, refpixel)
from pyrate.core.aps import wrap_spatio_temporal_filter
from pyrate.core.shared import Ifg, PrereadIfg, get_tiles, mpi_vs_multiprocess_logging
from pyrate.core.logger import pyratelogger as log
from pyrate.prepifg import find_header

MASTER_PROCESS = 0


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
    dest_tifs = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    ifgs_dict = {}
    process_tifs = mpiops.array_split(dest_tifs)
    for d in process_tifs:
        ifg = shared._prep_ifg(d, params)
        ifgs_dict[d] = PrereadIfg(path=d,
                                  nan_fraction=ifg.nan_fraction,
                                  master=ifg.master,
                                  slave=ifg.slave,
                                  time_span=ifg.time_span,
                                  nrows=ifg.nrows,
                                  ncols=ifg.ncols,
                                  metadata=ifg.meta_data)
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
    gt, md, wkt = shared.get_geotiff_header_info(process_tifs[0])
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


def _mst_calc(params):
    """
    MPI wrapper function for MST calculation
    """
    tiles = params['tiles']
    preread_ifgs = params['preread_ifgs']
    dest_tifs = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating minimum spanning tree matrix')

    def _save_mst_tile(tile, i, preread_ifgs):
        """
        Convenient inner loop for mst tile saving
        """
        mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs, params)
        # locally save the mst_mat
        mst_file_process_n = join(params[cf.TMPDIR], 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t in process_tiles:
        _save_mst_tile(t, t.index, preread_ifgs)
    log.debug('Finished mst calculation for process {}'.format(mpiops.rank))
    mpiops.comm.barrier()


def _ref_pixel_calc(params: dict) -> Tuple[int, int]:
    """
    Wrapper for reference pixel calculation
    """
    ifg_paths = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    lon = params[cf.REFX]
    lat = params[cf.REFY]

    ifg = Ifg(ifg_paths[0])
    ifg.open(readonly=True)
    # assume all interferograms have same projection and will share the same transform
    transform = ifg.dataset.GetGeoTransform()

    if lon == -1 or lat == -1:

        log.info('Searching for best reference pixel location')

        half_patch_size, thresh, grid = refpixel.ref_pixel_setup(ifg_paths, params)
        process_grid = mpiops.array_split(grid)
        refpixel.save_ref_pixel_blocks(process_grid, half_patch_size, ifg_paths, params)
        mean_sds = refpixel._ref_pixel_mpi(process_grid, half_patch_size, ifg_paths, thresh, params)
        mean_sds = mpiops.comm.gather(mean_sds, root=0)
        if mpiops.rank == MASTER_PROCESS:
            mean_sds = np.hstack(mean_sds)

        refpixel_returned = mpiops.run_once(refpixel.find_min_mean, mean_sds, grid)

        if isinstance(refpixel_returned, ValueError):
            from pyrate.core.refpixel import RefPixelError
            raise RefPixelError(
                "Reference pixel calculation returned an all nan slice!\n"
                "Cannot continue downstream computation. Please change reference pixel algorithm used before "
                "continuing.")
        refy, refx = refpixel_returned   # row first means first value is latitude
        log.info('Selected reference pixel coordinate (x, y): ({}, {})'.format(refx, refy))
        lon, lat = refpixel.convert_pixel_value_to_geographic_coordinate(refx, refy, transform)
        log.info('Selected reference pixel coordinate (lon, lat): ({}, {})'.format(lon, lat))

    else:
        log.info('Using reference pixel from config file (lon, lat): ({}, {})'.format(lon, lat))
        log.warning("Ensure user supplied reference pixel values are in lon/lat")
        refx, refy = refpixel.convert_geographic_coordinate_to_pixel_value(lon, lat, transform)
        log.info('Converted reference pixel coordinate (x, y): ({}, {})'.format(refx, refy))

    refpixel.update_refpix_metadata(ifg_paths, refx, refy, transform, params)

    log.debug("refpx, refpy: "+str(refx) + " " + str(refy))
    ifg.close()
    params[cf.REFX], params[cf.REFY] = int(refx), int(refy)
    return int(refx), int(refy)


def _orb_fit_calc(params: dict) -> None:
    """
    MPI wrapper for orbital fit correction
    """
    preread_ifgs = params['preread_ifgs']
    multi_paths = params[cf.INTERFEROGRAM_FILES]
    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required!')
        print('Orbital correction not required!')
        return
    log.info('Calculating orbital correction')

    ifg_paths = [p.sampled_path for p in multi_paths]
    if preread_ifgs:  # don't check except for mpi tests
        # perform some general error/sanity checks
        log.debug('Checking Orbital error correction status')
        if mpiops.run_once(shared.check_correction_status, ifg_paths, ifc.PYRATE_ORBITAL_ERROR):
            log.debug('Orbital error correction not required as all ifgs are already corrected!')
            return  # return if True condition returned

    if params[cf.ORBITAL_FIT_METHOD] == 1:
        prcs_ifgs = mpiops.array_split(ifg_paths)
        orbital.remove_orbital_error(prcs_ifgs, params, preread_ifgs)
    else:
        # Here we do all the multilooking in one process, but in memory
        # can use multiple processes if we write data to disc during
        # remove_orbital_error step
        # A performance comparison should be made for saving multilooked
        # files on disc vs in memory single process multilooking
        if mpiops.rank == MASTER_PROCESS:
            headers = [find_header(p, params) for p in multi_paths]
            orbital.remove_orbital_error(ifg_paths, params, headers, preread_ifgs=preread_ifgs)
    mpiops.comm.barrier()
    shared.save_numpy_phase(ifg_paths, params)
    log.debug('Finished Orbital error correction')


def _ref_phase_estimation(params):
    """
    Wrapper for reference phase estimation.
    """
    refpx, refpy = params[cf.REFX], params[cf.REFY]
    return _ref_phase_estimation_wrapper(params, refpx, refpy)


def _ref_phase_estimation_wrapper(params, refpx, refpy):
    """
    Wrapper for reference phase estimation.
    """
    ifg_paths = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    log.info("Calculating reference phase")
    if len(ifg_paths) < 2:
        raise rpe.ReferencePhaseError(
            "At least two interferograms required for reference phase correction ({len_ifg_paths} "
            "provided).".format(len_ifg_paths=len(ifg_paths))
        )

    if mpiops.run_once(shared.check_correction_status, ifg_paths, ifc.PYRATE_REF_PHASE):
        log.debug('Finished reference phase correction')
        return

    if params[cf.REF_EST_METHOD] == 1:
        ref_phs = rpe.est_ref_phase_method1(ifg_paths, params)
    elif params[cf.REF_EST_METHOD] == 2:
        ref_phs = rpe.est_ref_phase_method2(ifg_paths, params, refpx, refpy)
    else:
        raise rpe.ReferencePhaseError("No such option, use '1' or '2'.")

    # Save reference phase numpy arrays to disk.
    ref_phs_file = os.path.join(params[cf.TMPDIR], 'ref_phs.npy')
    if mpiops.rank == MASTER_PROCESS:
        collected_ref_phs = np.zeros(len(ifg_paths), dtype=np.float64)
        process_indices = mpiops.array_split(range(len(ifg_paths)))
        collected_ref_phs[process_indices] = ref_phs
        for r in range(1, mpiops.size):
            process_indices = mpiops.array_split(range(len(ifg_paths)), r)
            this_process_ref_phs = np.zeros(shape=len(process_indices),
                                            dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=r, tag=r)
            collected_ref_phs[process_indices] = this_process_ref_phs
        np.save(file=ref_phs_file, arr=collected_ref_phs)
    else:
        mpiops.comm.Send(ref_phs, dest=MASTER_PROCESS, tag=mpiops.rank)
    log.debug('Finished reference phase correction')

    # Preserve old return value so tests don't break.
    if isinstance(ifg_paths[0], Ifg):
        ifgs = ifg_paths
    else:
        ifgs = [Ifg(ifg_path) for ifg_path in ifg_paths]
    mpiops.comm.barrier()
    shared.save_numpy_phase(ifg_paths, params)
    return ref_phs, ifgs


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

    return process_ifgs(params)


def _stack_calc(params):
    """
    MPI wrapper for stacking calculation
    """
    tiles = params[cf.TILES]
    preread_ifgs = params[cf.PREREAD_IFGS]
    vcmt = params[cf.VCMT]
    ifg_paths = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating rate map from stacking')
    output_dir = params[cf.TMPDIR]
    for t in process_tiles:
        log.info('Stacking of tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
        mst_grid_n = np.load(os.path.join(output_dir, 'mst_mat_{}.npy'.format(t.index)))
        rate, error, samples = stack.stack_rate_array(ifg_parts, params, vcmt, mst_grid_n)
        # declare file names
        np.save(file=os.path.join(output_dir, 'stack_rate_{}.npy'.format(t.index)), arr=rate)
        np.save(file=os.path.join(output_dir, 'stack_error_{}.npy'.format(t.index)), arr=error)
        np.save(file=os.path.join(output_dir, 'stack_samples_{}.npy'.format(t.index)), arr=samples)
    mpiops.comm.barrier()
    log.debug("Finished stack rate calc!")


def _maxvar_vcm_calc(params):
    """
    MPI wrapper for maxvar and vcmt computation
    """
    preread_ifgs = params[cf.PREREAD_IFGS]
    ifg_paths = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    log.info('Calculating the temporal variance-covariance matrix')
    process_indices = mpiops.array_split(range(len(ifg_paths)))

    def _get_r_dist(ifg_path):
        """
        Get RDIst class object
        """
        ifg = Ifg(ifg_path)
        ifg.open()
        r_dist = vcm_module.RDist(ifg)()
        ifg.close()
        return r_dist

    r_dist = mpiops.run_once(_get_r_dist, ifg_paths[0])
    prcs_ifgs = mpiops.array_split(ifg_paths)
    process_maxvar = []
    for n, i in enumerate(prcs_ifgs):
        log.debug('Calculating maxvar for {} of process ifgs {} of total {}'.format(n+1, len(prcs_ifgs), len(ifg_paths)))
        process_maxvar.append(vcm_module.cvd(i, params, r_dist, calc_alpha=True, write_vals=True, save_acg=True)[0])
    if mpiops.rank == MASTER_PROCESS:
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        maxvar[process_indices] = process_maxvar
        for i in range(1, mpiops.size):  # pragma: no cover
            rank_indices = mpiops.array_split(range(len(ifg_paths)), i)
            this_process_ref_phs = np.empty(len(rank_indices), dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=i, tag=i)
            maxvar[rank_indices] = this_process_ref_phs
    else:  # pragma: no cover
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        mpiops.comm.Send(np.array(process_maxvar, dtype=np.float64), dest=MASTER_PROCESS, tag=mpiops.rank)

    mpiops.comm.barrier()
    maxvar = mpiops.comm.bcast(maxvar, root=0)
    vcmt = mpiops.run_once(vcm_module.get_vcmt, preread_ifgs, maxvar)
    log.debug("Finished maxvar and vcm calc!")
    params[cf.MAXVAR], params[cf.VCMT] = maxvar, vcmt
    return maxvar, vcmt


def _timeseries_calc(params):
    """
    MPI wrapper for time series calculation.
    """
    tiles = params[cf.TILES]
    preread_ifgs = params[cf.PREREAD_IFGS]
    vcmt = params[cf.VCMT]
    ifg_paths = [ifg_path.sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    if params[cf.TIME_SERIES_CAL] == 0:
        log.info('Time Series Calculation not required')
        return

    if params[cf.TIME_SERIES_METHOD] == 1:
        log.info('Calculating time series using Laplacian Smoothing method')
    elif params[cf.TIME_SERIES_METHOD] == 2:
        log.info('Calculating time series using SVD method')

    output_dir = params[cf.TMPDIR]
    total_tiles = len(tiles)
    process_tiles = mpiops.array_split(tiles)
    for t in process_tiles:
        log.debug("Calculating time series for tile "+str(t.index)+" out of "+str(total_tiles))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir, 'mst_mat_{}.npy'.format(t.index)))
        tsincr, tscuml, _ = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        np.save(file=os.path.join(output_dir, 'tscuml_{}.npy'.format(t.index)), arr=tscuml)
        # optional save of tsincr npy tiles
        if params["savetsincr"] == 1:
            np.save(file=os.path.join(output_dir, 'tsincr_{}.npy'.format(t.index)), arr=tsincr)
    mpiops.comm.barrier()
    log.debug("Finished timeseries calc!")


def _update_params_with_tiles(params: dict) -> None:
    ifg_path = params[cf.INTERFEROGRAM_FILES][0].sampled_path
    rows, cols = params["rows"], params["cols"]
    tiles = mpiops.run_once(get_tiles, ifg_path, rows, cols)
    # add tiles to params
    params['tiles'] = tiles


process_steps = {
    'orbfit': _orb_fit_calc,
    'refphase': _ref_phase_estimation,
    'mst': _mst_calc,
    'apscorrect': wrap_spatio_temporal_filter,
    'maxvar': _maxvar_vcm_calc,
    'timeseries': _timeseries_calc,
    'stack': _stack_calc
}


def process_ifgs(params: dict):
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

    if mpiops.size > 1:  # turn of multiprocessing during mpi jobs
        params[cf.PARALLEL] = False
    _update_params_with_tiles(params)
    _create_ifg_dict(params)
    refpx, refpy = _ref_pixel_calc(params)

    for step in process_steps:
        process_steps[step](params)
    log.info('PyRate workflow completed')
    return (refpx, refpy), params[cf.MAXVAR], params[cf.VCMT]
