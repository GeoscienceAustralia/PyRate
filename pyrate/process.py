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
import logging
import os
from os.path import join
import pickle as cp
from collections import OrderedDict
import numpy as np

from pyrate.core import (shared, algorithm, orbital, ref_phs_est as rpe, 
                         ifgconstants as ifc, mpiops, config as cf, 
                         timeseries, mst, covariance as vcm_module, 
                         stack, refpixel)
from pyrate.core.aps import wrap_spatio_temporal_filter
from pyrate.core.shared import Ifg, PrereadIfg, get_tiles, mpi_vs_multiprocess_logging
from pyrate.core.logger import pyratelogger as log

MASTER_PROCESS = 0


def _join_dicts(dicts):
    """
    Function to concatenate dictionaries
    """
    if dicts is None:  # pragma: no cover
        return
    assembled_dict = {k: v for D in dicts for k, v in D.items()}
    return assembled_dict


def _create_ifg_dict(dest_tifs, params, tiles):
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
    ifgs_dict = {}
    nifgs = len(dest_tifs)
    process_tifs = mpiops.array_split(dest_tifs)
    shared.save_numpy_phase(dest_tifs, tiles, params)
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

    preread_ifgs_file = join(params[cf.TMPDIR], 'preread_ifgs.pk')

    if mpiops.rank == MASTER_PROCESS:

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

    mpiops.comm.barrier()
    preread_ifgs = OrderedDict(sorted(cp.load(open(preread_ifgs_file, 'rb')).items()))
    log.debug('Finished converting phase_data to numpy in process {}'.format(mpiops.rank))
    return preread_ifgs


def _mst_calc(dest_tifs, params, tiles, preread_ifgs):
    """
    MPI wrapper function for MST calculation
    """
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


def _ref_pixel_calc(ifg_paths, params):
    """
    Wrapper for reference pixel calculation
    """
    refx = params[cf.REFX]
    refy = params[cf.REFY]

    ifg = Ifg(ifg_paths[0])
    ifg.open(readonly=True)

    if refx == -1 or refy == -1:

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

        refy, refx = refpixel_returned

        log.info('Selected reference pixel coordinate: ({}, {})'.format(refx, refy))
    else:
        log.info('Reusing reference pixel from config file: ({}, {})'.format(refx, refy))
    ifg.close()
    return refx, refy


def _orb_fit_calc(ifg_paths, params, preread_ifgs=None):
    """
    MPI wrapper for orbital fit correction
    """
    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required!')
        print('Orbital correction not required!')
        return
    log.info('Calculating orbital correction')

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
            orbital.remove_orbital_error(ifg_paths, params, preread_ifgs)
    mpiops.comm.barrier()
    log.debug('Finished Orbital error correction')


def _ref_phase_estimation(ifg_paths, params, refpx, refpy):
    """
    Wrapper for reference phase estimation.
    """
    log.info("Calculating reference phase")
    if len(ifg_paths) < 2:
        raise rpe.ReferencePhaseError(
            "At least two interferograms required for reference phase correction ({len_ifg_paths} "
            "provided).".format(len_ifg_paths=len(ifg_paths))
        )

    if mpiops.run_once(shared.check_correction_status, ifg_paths, ifc.PYRATE_REF_PHASE):
        log.debug('Finished reference phase estimation')
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
    log.debug('Finished reference phase estimation')

    # Preserve old return value so tests don't break.
    if isinstance(ifg_paths[0], Ifg):
        ifgs = ifg_paths
    else:
        ifgs = [Ifg(ifg_path) for ifg_path in ifg_paths]
    mpiops.comm.barrier()
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

    ifg_paths = []
    for ifg_path in params[cf.INTERFEROGRAM_FILES]:
        ifg_paths.append(ifg_path.sampled_path)

    rows, cols = params["rows"], params["cols"]

    return process_ifgs(ifg_paths, params, rows, cols)


def process_ifgs(ifg_paths, params, rows, cols):
    """
    Top level function to perform PyRate workflow on given interferograms

    :param list ifg_paths: List of interferogram paths
    :param dict params: Dictionary of configuration parameters
    :param int rows: Number of sub-tiles in y direction
    :param int cols: Number of sub-tiles in x direction

    :return: refpt: tuple of reference pixel x and y position
    :rtype: tuple
    :return: maxvar: array of maximum variance values of interferograms
    :rtype: ndarray
    :return: vcmt: Variance-covariance matrix array
    :rtype: ndarray
    """

    if mpiops.size > 1:  # turn of multiprocessing during mpi jobs
        params[cf.PARALLEL] = False

    tiles = mpiops.run_once(get_tiles, ifg_paths[0], rows, cols)

    preread_ifgs = _create_ifg_dict(ifg_paths, params=params, tiles=tiles)
    # _mst_calc(ifg_paths, params, tiles, preread_ifgs)

    refpx, refpy = _ref_pixel_calc(ifg_paths, params)

    log.debug("refpx, refpy: "+str(refpx) + " " + str(refpy))

    # remove non ifg keys
    _ = [preread_ifgs.pop(k) for k in ['gt', 'epochlist', 'md', 'wkt']]

    _orb_fit_calc(ifg_paths, params, preread_ifgs)

    _ref_phase_estimation(ifg_paths, params, refpx, refpy)

    _mst_calc(ifg_paths, params, tiles, preread_ifgs)

    # spatio-temporal aps filter
    wrap_spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs)

    maxvar, vcmt = _maxvar_vcm_calc(ifg_paths, params, preread_ifgs)
    # save phase data tiles as numpy array for timeseries and stackrate calc

    shared.save_numpy_phase(ifg_paths, tiles, params)

    _timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    _stack_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    log.info('PyRate workflow completed')
    return (refpx, refpy), maxvar, vcmt


def _stack_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    """
    MPI wrapper for stacking calculation
    """
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating rate map from stacking')
    output_dir = params[cf.TMPDIR]
    for t in process_tiles:
        log.info('Stacking of tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs, params) for p in ifg_paths]
        mst_grid_n = np.load(os.path.join(output_dir, 'mst_mat_{}.npy'.format(t.index)))
        rate, error, samples = stack.stack_rate(ifg_parts, params, vcmt, mst_grid_n)
        # declare file names
        np.save(file=os.path.join(output_dir, 'stack_rate_{}.npy'.format(t.index)), arr=rate)
        np.save(file=os.path.join(output_dir, 'stack_error_{}.npy'.format(t.index)), arr=error)
        np.save(file=os.path.join(output_dir, 'stack_samples_{}.npy'.format(t.index)), arr=samples)
    mpiops.comm.barrier()
    log.debug("Finished stack rate calc!")


def _maxvar_vcm_calc(ifg_paths, params, preread_ifgs):
    """
    MPI wrapper for maxvar and vcmt computation
    """
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
    return maxvar, vcmt


def _timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    """
    MPI wrapper for time series calculation.
    """
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
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, _ = res
        np.save(file=os.path.join(output_dir, 'tsincr_{}.npy'.format(t.index)), arr=tsincr)
        np.save(file=os.path.join(output_dir, 'tscuml_{}.npy'.format(t.index)), arr=tscum)
    mpiops.comm.barrier()
    log.debug("Finished timeseries calc!")
