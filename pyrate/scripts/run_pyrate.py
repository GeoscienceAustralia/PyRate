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

from pyrate import algorithm
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import linrate
from pyrate import mpiops
from pyrate import mst
from pyrate import orbital
from pyrate import ref_phs_est as rpe
from pyrate import refpixel
from pyrate import shared
from pyrate import timeseries
from pyrate import covariance as vcm_module
from pyrate.aps import _wrap_spatio_temporal_filter
from pyrate.config import ConfigException
from pyrate.shared import Ifg, PrereadIfg, get_tiles

MASTER_PROCESS = 0
log = logging.getLogger(__name__)


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
        ifgs_dict['epochlist'] = algorithm.get_epochs(ifgs_dict)[0]
        ifgs_dict['gt'] = gt
        ifgs_dict['md'] = md
        ifgs_dict['wkt'] = wkt
        # dump ifgs_dict file for later use
        cp.dump(ifgs_dict, open(preread_ifgs_file, 'wb'))

    mpiops.comm.barrier()
    preread_ifgs = OrderedDict(sorted(cp.load(open(preread_ifgs_file,
                                                   'rb')).items()))
    log.info('Finished converting phase_data to numpy '
             'in process {}'.format(mpiops.rank))
    return preread_ifgs


def _mst_calc(dest_tifs, params, tiles, preread_ifgs):
    """
    MPI wrapper function for MST calculation
    """
    process_tiles = mpiops.array_split(tiles)

    def _save_mst_tile(tile, i, preread_ifgs):
        """
        Convenient inner loop for mst tile saving
        """
        if params[cf.NETWORKX_OR_MATLAB_FLAG] == 1:
            log.info('Calculating minimum spanning tree matrix '
                     'using NetworkX method')
            mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        elif params[cf.NETWORKX_OR_MATLAB_FLAG] == 0:
            raise ConfigException('Matlab-style MST not supported')
        else:
            raise ConfigException('Only NetworkX MST is supported')
            # mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        # locally save the mst_mat
        mst_file_process_n = join(
            params[cf.TMPDIR], 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t in process_tiles:
        _save_mst_tile(t, t.index, preread_ifgs)
    log.info('finished mst calculation for process {}'.format(mpiops.rank))
    mpiops.comm.barrier()


def _ref_pixel_calc(ifg_paths, params):
    """
    Wrapper for reference pixel calculation
    """
    # unlikely, but possible the refpixel can be (0,0)
    # check if there is a pre-specified reference pixel coord
    refx = params[cf.REFX]
    ifg = Ifg(ifg_paths[0])
    ifg.open(readonly=True)
    if refx > ifg.ncols - 1:
        msg = ('Supplied reference pixel X coordinate is greater than '
               'the number of ifg columns: {}').format(refx)
        raise ValueError(msg)

    refy = params[cf.REFY]
    if refy > ifg.nrows - 1:
        msg = ('Supplied reference pixel Y coordinate is greater than '
               'the number of ifg rows: {}').format(refy)
        raise ValueError(msg)

    if refx <= 0 or refy <= 0:  # if either zero or negative
        log.info('Searching for best reference pixel location')

        half_patch_size, thresh, grid = refpixel.ref_pixel_setup(ifg_paths, 
                                                                 params)
        process_grid = mpiops.array_split(grid)
        refpixel.save_ref_pixel_blocks(process_grid, half_patch_size,
                                       ifg_paths, params)
        mean_sds = refpixel._ref_pixel_mpi(process_grid, half_patch_size,
                                           ifg_paths, thresh, params)
        mean_sds = mpiops.comm.gather(mean_sds, root=0)
        if mpiops.rank == MASTER_PROCESS:
            mean_sds = np.hstack(mean_sds)

        refy, refx = mpiops.run_once(refpixel.find_min_mean, mean_sds, grid)
        log.info('Selected reference pixel coordinate: '
                 '({}, {})'.format(refx, refy))
    else:  # pragma: no cover
        log.info('Reusing reference pixel from config file: '
                 '({}, {})'.format(refx, refy))
    ifg.close()
    return refx, refy


def _orb_fit_calc(ifg_paths, params, preread_ifgs=None):
    """
    MPI wrapper for orbital fit correction
    """
    log.info('Calculating orbfit correction')

    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required')
        return

    if preread_ifgs:  # don't check except for mpi tests
        # perform some general error/sanity checks
        log.info('Checking Orbital error correction status')
        if mpiops.run_once(shared.check_correction_status, preread_ifgs,
                           ifc.PYRATE_ORBITAL_ERROR):
            log.info('Finished Orbital error correction')
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
    log.info('Finished Orbital error correction')


def _ref_phase_estimation(ifg_paths, params, refpx, refpy, preread_ifgs=None):
    """
    Wrapper function for MPI-enabled reference phase estimation.
    Refer to documentation for ref_est_phs.estimate_ref_phase.
    """
    # perform some checks on existing ifgs
    #if preread_ifgs and mpiops.rank == MASTER_PROCESS:
    #TODO: implement MPI capability into ref_phs_est module and remove here
    if preread_ifgs:
        log.info('Checking reference phase estimation status')
        if mpiops.run_once(shared.check_correction_status, preread_ifgs,
                           ifc.PYRATE_REF_PHASE):
            log.info('Finished reference phase estimation')
            return  # return if True condition returned

    if params[cf.REF_EST_METHOD] == 1:
        # calculate phase sum for later use in ref phase method 1
        comp = _phase_sum(ifg_paths, params)
        log.info('Computing reference phase via method 1')
        process_ref_phs = _ref_phs_method1(ifg_paths, comp)
    elif params[cf.REF_EST_METHOD] == 2:
        log.info('Computing reference phase via method 2')
        process_ref_phs = _ref_phs_method2(ifg_paths, params, refpx, refpy)
    else:
        raise ConfigException('Ref phase estimation method must be 1 or 2')

    # Save reference phase numpy arrays to disk
    ref_phs_file = join(params[cf.TMPDIR], 'ref_phs.npy')
    if mpiops.rank == MASTER_PROCESS:
        ref_phs = np.zeros(len(ifg_paths), dtype=np.float64)
        process_indices = mpiops.array_split(range(len(ifg_paths)))
        ref_phs[process_indices] = process_ref_phs
        for r in range(1, mpiops.size):  # pragma: no cover
            process_indices = mpiops.array_split(range(len(ifg_paths)), r)
            this_process_ref_phs = np.zeros(shape=len(process_indices),
                                            dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=r, tag=r)
            ref_phs[process_indices] = this_process_ref_phs
        np.save(file=ref_phs_file, arr=ref_phs)
    else:  # pragma: no cover
        # send reference phase data to master process
        mpiops.comm.Send(process_ref_phs, dest=MASTER_PROCESS,
                         tag=mpiops.rank)
    log.info('Finished reference phase estimation')


def _phase_sum(ifg_paths, params):
    """
    Save phase data and phase sum used in the reference phase estimation
    """
    p_paths = mpiops.array_split(ifg_paths)
    ifg = Ifg(p_paths[0])
    ifg.open(readonly=True)
    shape = ifg.shape
    phs_sum = np.zeros(shape=shape, dtype=np.float64)
    ifg.close()

    for d in p_paths:
        ifg = Ifg(d)
        ifg.open()
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        phs_sum += ifg.phase_data
        ifg.close()

    if mpiops.rank == MASTER_PROCESS:
        phase_sum_all = phs_sum
        # loop is better for memory
        for i in range(1, mpiops.size):  # pragma: no cover
            phs_sum = np.zeros(shape=shape, dtype=np.float64)
            mpiops.comm.Recv(phs_sum, source=i, tag=i)
            phase_sum_all += phs_sum
        comp = np.isnan(phase_sum_all)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
    else:  # pragma: no cover
        comp = None
        mpiops.comm.Send(phs_sum, dest=0, tag=mpiops.rank)

    comp = mpiops.comm.bcast(comp, root=0)
    return comp


def _ref_phs_method2(ifg_paths, params, refpx, refpy):
    """
    MPI wrapper for reference phase computation using method 2.
    Refer to documentation for ref_est_phs.est_ref_phase_method2.
    """
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    process_ifg_paths = mpiops.array_split(ifg_paths)

    def _inner(ifg_path):
        """
        Convenient inner loop
        """
        ifg = Ifg(ifg_path)
        ifg.open(readonly=False)
        phase_data = ifg.phase_data
        ref_ph = rpe._est_ref_phs_method2(phase_data,
                                         half_chip_size,
                                         refpx, refpy, thresh)
        phase_data -= ref_ph
        md = ifg.meta_data
        md[ifc.PYRATE_REF_PHASE] = ifc.REF_PHASE_REMOVED
        ifg.write_modified_phase(data=phase_data)
        ifg.close()
        return ref_ph

    ref_phs = np.array([_inner(p) for p in process_ifg_paths])
    log.info('Ref phase computed in process {}'.format(mpiops.rank))
    return ref_phs


def _ref_phs_method1(ifg_paths, comp):
    """
    MPI wrapper for reference phase computation using method 1.
    Refer to documentation for ref_est_phs.est_ref_phase_method1.
    """
    def _inner(ifg_path):
        """
        Convenient inner loop
        """
        ifg = Ifg(ifg_path)
        ifg.open(readonly=False)
        phase_data = ifg.phase_data
        ref_phase = rpe._est_ref_phs_method1(phase_data, comp)
        phase_data -= ref_phase
        md = ifg.meta_data
        md[ifc.PYRATE_REF_PHASE] = ifc.REF_PHASE_REMOVED
        ifg.write_modified_phase(data=phase_data)
        ifg.close()
        return ref_phase
    this_process_ifgs = mpiops.array_split(ifg_paths)
    ref_phs = np.array([_inner(ifg) for ifg in this_process_ifgs])
    log.info('Ref phase computed in process {}'.format(mpiops.rank))
    return ref_phs


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


    # remove non ifg keys
    _ = [preread_ifgs.pop(k) for k in ['gt', 'epochlist', 'md', 'wkt']]

    _orb_fit_calc(ifg_paths, params, preread_ifgs)

    _ref_phase_estimation(ifg_paths, params, refpx, refpy, preread_ifgs)

    _mst_calc(ifg_paths, params, tiles, preread_ifgs)

    # spatio-temporal aps filter
    _wrap_spatio_temporal_filter(ifg_paths, params, tiles, preread_ifgs)

    maxvar, vcmt = _maxvar_vcm_calc(ifg_paths, params, preread_ifgs)

    # save phase data tiles as numpy array for timeseries and linrate calc
    shared.save_numpy_phase(ifg_paths, tiles, params)

    _timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    _linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    log.info('PyRate workflow completed')
    return (refpx, refpy), maxvar, vcmt


def _linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    """
    MPI wrapper for linrate calculation
    """
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating linear rate map')
    output_dir = params[cf.TMPDIR]
    for t in process_tiles:
        log.info('Calculating linear rate of tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_grid_n = np.load(os.path.join(output_dir,
                                          'mst_mat_{}.npy'.format(t.index)))
        rate, error, samples = linrate.linear_rate(ifg_parts, params,
                                                   vcmt, mst_grid_n)
        # declare file names
        np.save(file=os.path.join(output_dir,
                                  'linrate_{}.npy'.format(t.index)),
                arr=rate)
        np.save(file=os.path.join(output_dir,
                                  'linerror_{}.npy'.format(t.index)),
                arr=error)
        np.save(file=os.path.join(output_dir,
                                  'linsamples_{}.npy'.format(t.index)),
                arr=samples)
    mpiops.comm.barrier()


def _maxvar_vcm_calc(ifg_paths, params, preread_ifgs):
    """
    MPI wrapper for maxvar and vcmt computation
    """
    log.info('Calculating maxvar and vcm')
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
        log.info('Calculating maxvar for {} of process ifgs {} of '
                 'total {}'.format(n+1, len(prcs_ifgs), len(ifg_paths)))
        process_maxvar.append(vcm_module.cvd(i, params, r_dist,
                                             calc_alpha=True,
                                             write_vals=True,
                                             save_acg=True)[0])
    if mpiops.rank == MASTER_PROCESS:
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        maxvar[process_indices] = process_maxvar
        for i in range(1, mpiops.size):  # pragma: no cover
            rank_indices = mpiops.array_split(range(len(ifg_paths)), i)
            this_process_ref_phs = np.empty(len(rank_indices),
                                            dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=i, tag=i)
            maxvar[rank_indices] = this_process_ref_phs
    else:  # pragma: no cover
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        mpiops.comm.Send(np.array(process_maxvar, dtype=np.float64),
                         dest=MASTER_PROCESS, tag=mpiops.rank)

    maxvar = mpiops.comm.bcast(maxvar, root=0)
    vcmt = mpiops.run_once(vcm_module.get_vcmt, preread_ifgs, maxvar)
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
    process_tiles = mpiops.array_split(tiles)
    for t in process_tiles:
        log.info('Calculating time series for tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir, 'mst_mat_{}.npy'.format(t.index)))
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, _ = res
        np.save(file=os.path.join(output_dir, 'tsincr_{}.npy'.format(t.index)), arr=tsincr)
        np.save(file=os.path.join(output_dir, 'tscuml_{}.npy'.format(t.index)), arr=tscum)
    mpiops.comm.barrier()
