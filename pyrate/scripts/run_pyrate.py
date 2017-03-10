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
This Python module runs the main PyRate processing workflow
"""
from __future__ import print_function

import logging
import os
from os.path import join
import pickle as cp
import numpy as np

from pyrate import algorithm
from pyrate import config as cf
from pyrate.config import ConfigException
from pyrate import ifgconstants as ifc
from pyrate import linrate
from pyrate import mpiops
from pyrate import mst
from pyrate import orbital
from pyrate import ref_phs_est as rpe
from pyrate import refpixel
from pyrate import shared
from pyrate import timeseries
from pyrate import vcm as vcm_module
from pyrate.compat import PyAPS_INSTALLED
from pyrate.shared import Ifg, create_tiles, \
    PrereadIfg, prepare_ifg, save_numpy_phase, get_projection_info

if PyAPS_INSTALLED:  # pragma: no cover
    # from pyrate import aps
    from pyrate.aps import check_aps_ifgs, aps_delay_required

MASTER_PROCESS = 0
log = logging.getLogger(__name__)


def get_tiles(ifg_path, rows, cols):
    """
    Break up the ifgs into tiles based on user supplied rows and cols

    Parameters
    ----------
    ifg_path: str
        list of destination tifs
    rows: int
        number of rows to break each ifg into
    cols: int
        number of cols to break each ifg into

    Returns
    -------
    tiles: list
        list of shared.Tile instances
    """
    ifg = Ifg(ifg_path)
    ifg.open(readonly=True)
    tiles = create_tiles(ifg.shape, nrows=rows, ncols=cols)
    ifg.close()
    return tiles


def _join_dicts(dicts):
    """
    Parameters
    ----------
    dicts: list
        list of dicts to join

    Returns
    -------
    assembled_dict: dict
        dictionary after join
    """
    if dicts is None:  # pragma: no cover
        return
    assembled_dict = {k: v for D in dicts for k, v in D.items()}
    return assembled_dict


def create_ifg_dict(dest_tifs, params, tiles):
    """
    1. Convert ifg phase data into numpy binary files.
    2. Save the preread_ifgs dict with information about the ifgs that are
    later used for fast loading of Ifg files in IfgPart class

    Parameters
    ----------
    dest_tifs: list
        list of destination tifs
    params: dict
        config dict
    tiles: list
        list of all Tile instances

    Returns
    -------
    preread_ifgs: dict
        dict containing information regarding ifgs that are used downstream
    """
    ifgs_dict = {}
    process_tifs = mpiops.array_split(dest_tifs)
    save_numpy_phase(dest_tifs, tiles, params)
    for d in process_tifs:
        ifg = prepare_ifg(d, params)
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

    preread_ifgs_file = join(params[cf.OUT_DIR], 'preread_ifgs.pk')

    if mpiops.rank == MASTER_PROCESS:

        # add some extra information that's also useful later
        gt, md, wkt = get_projection_info(process_tifs[0])
        ifgs_dict['epochlist'] = algorithm.get_epochs(ifgs_dict)[0]
        ifgs_dict['gt'] = gt
        ifgs_dict['md'] = md
        ifgs_dict['wkt'] = wkt
        # dump ifgs_dict file for later use
        cp.dump(ifgs_dict, open(preread_ifgs_file, 'wb'))

    mpiops.comm.barrier()
    preread_ifgs = cp.load(open(preread_ifgs_file, 'rb'))
    log.info('finish converting phase_data to numpy '
             'in process {}'.format(mpiops.rank))
    return preread_ifgs


def mst_calc(dest_tifs, params, tiles, preread_ifgs):
    """
    MPI function that control each process during MPI run
    Reference phase computation using method 2
    Parameters
    ----------
    dest_tifs: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    tiles: list
        list of all tiles used during mpi processes
    preread_ifgs: dict
        dict containing ifg characteristics for efficient computing
    """
    process_tiles = mpiops.array_split(tiles)

    def save_mst_tile(tile, i, preread_ifgs):
        """ Convenient inner loop for mst tile saving"""
        if params[cf.NETWORKX_OR_MATLAB_FLAG] == 1:
            log.info('Calculating minimum spanning tree matrix '
                     'using NetworkX method')
            mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        elif params[cf.NETWORKX_OR_MATLAB_FLAG] == 0:
            raise ConfigException('Matlab mst not supported')
        else:
            raise ConfigException('Only NetworkX mst is supported')
            # mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        # locally save the mst_mat
        mst_file_process_n = join(
            params[cf.OUT_DIR], 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t in process_tiles:
        save_mst_tile(t, t.index, preread_ifgs)
    log.info('finished mst calculation for process {}'.format(mpiops.rank))
    mpiops.comm.barrier()


def ref_pixel_calc(ifg_paths, params):
    """
    Reference pixel calculation setup

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file

    Returns
    -------
    refx: float
        reference pixel x-coordinate
    refy: float
        reference pixel y-coordinate

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
        refy, refx = find_ref_pixel(ifg_paths, params)
        log.info('Selected reference pixel coordinate: '
                 '({}, {})'.format(refx, refy))
    else:  # pragma: no cover
        log.info('Reusing reference pixel from config file: '
                 '({}, {})'.format(refx, refy))
    ifg.close()
    return refx, refy


def find_ref_pixel(ifg_paths, params):
    """
    Find reference pixel using mpi
    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file

    Returns
    -------
    tuple of (refy, refx)
    """
    half_patch_size, thresh, grid = refpixel.ref_pixel_setup(ifg_paths, params)
    process_grid = mpiops.array_split(grid)
    save_ref_pixel_blocks(process_grid, half_patch_size, ifg_paths, params)
    mean_sds = refpixel.ref_pixel_mpi(process_grid, half_patch_size,
                                      ifg_paths, thresh, params)
    mean_sds = mpiops.comm.gather(mean_sds, root=0)
    if mpiops.rank == MASTER_PROCESS:
        mean_sds = np.hstack(mean_sds)
    return mpiops.run_once(refpixel.filter_means, mean_sds, grid)


def save_ref_pixel_blocks(grid, half_patch_size, ifg_paths, params):
    """
    Parameters
    ----------
    grid: list
        list of tuples (y, x) corresponding reference pixel grids
    half_patch_size: int
        patch size in pixels corresponding to ref pixel grids
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    """
    log.info('Saving ref pixel blocks')
    outdir = params[cf.OUT_DIR]
    for pth in ifg_paths:
        ifg = Ifg(pth)
        ifg.open(readonly=True)
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.convert_to_mm()
        for y, x in grid:
            data = ifg.phase_data[y - half_patch_size:y + half_patch_size + 1,
                                  x - half_patch_size:x + half_patch_size + 1]

            data_file = join(outdir, 'ref_phase_data_{b}_{y}_{x}.npy'.format(
                b=os.path.basename(pth).split('.')[0], y=y, x=x))
            np.save(file=data_file, arr=data)
        ifg.close()
    log.info('Saved ref pixel blocks')


def orb_fit_calc(ifg_paths, params, preread_ifgs=None):
    """
    Orbital fit correction

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    """
    log.info('Calculating orbfit correction')
    if params[cf.ORBITAL_FIT_METHOD] == 1:
        prcs_ifgs = mpiops.array_split(ifg_paths)
        orbital.remove_orbital_error(prcs_ifgs, params, preread_ifgs)
    else:
        # Here we do all the multilooking in one process, but in memory
        # can use multiple processes if we write data to disc during
        # remove_orbital_error step
        # A performance comparison should be performed be performed for saving
        # multilooked files on disc vs in memory single process multilooking
        if mpiops.rank == MASTER_PROCESS:
            orbital.remove_orbital_error(ifg_paths, params, preread_ifgs)
    mpiops.comm.barrier()
    log.info('Finished orbfit calculation in process {}'.format(mpiops.rank))


def ref_phase_estimation(ifg_paths, params, refpx, refpy):
    """
    Reference phase estimation

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    refpx: float
        reference pixel x-coordinate
    refpy: float
        reference pixel y-coordinate
    """

    log.info('Estimating and removing reference phase')
    if params[cf.REF_EST_METHOD] == 1:
        # calculate phase sum for later use in ref phase method 1
        comp = phase_sum(ifg_paths, params)
        process_ref_phs = ref_phs_method1(ifg_paths, comp)
    elif params[cf.REF_EST_METHOD] == 2:
        process_ref_phs = ref_phs_method2(ifg_paths, params, refpx, refpy)
    else:
        raise ConfigException('Ref phase estimation method must be 1 or 2')

    ref_phs_file = join(params[cf.OUT_DIR], 'ref_phs.npy')
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


def ref_phs_method2(ifg_paths, params, refpx, refpy):
    """
    Reference phase computation using method 2
    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    refpx: float
        reference pixel x-coordinate
    refpy: float
        reference pixel y-coordinate

    Returns
    -------
    ref_phs: ndarray
        array of reference phase of shape ifg.shape
    """
    half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
    chipsize = 2 * half_chip_size + 1
    thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
    process_ifg_paths = mpiops.array_split(ifg_paths)

    def _inner(ifg_path):
        ifg = Ifg(ifg_path)
        ifg.open(readonly=False)
        phase_data = ifg.phase_data
        ref_ph = rpe.est_ref_phs_method2(phase_data,
                                         half_chip_size,
                                         refpx, refpy, thresh)
        phase_data -= ref_ph
        md = ifg.meta_data
        md[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
        ifg.write_modified_phase(data=phase_data)
        ifg.close()
        return ref_ph

    ref_phs = np.array([_inner(p) for p in process_ifg_paths])
    log.info('Ref phase computed in process {}'.format(mpiops.rank))
    return ref_phs


def ref_phs_method1(ifg_paths, comp):
    """
    Reference phase computation using method 1

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    comp: ndarray
        array of phase sum of all ifgs of shape ifg.shape

    Returns
    -------
    ref_phs: ndarray
        array of reference phase of shape ifg.shape
    """

    def _inner(ifg_path):
        """ convenient inner loop """
        ifg = Ifg(ifg_path)
        ifg.open(readonly=False)
        phase_data = ifg.phase_data
        ref_phase = rpe.est_ref_phs_method1(phase_data, comp)
        phase_data -= ref_phase
        md = ifg.meta_data
        md[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
        ifg.write_modified_phase(data=phase_data)
        ifg.close()
        return ref_phase
    this_process_ifgs = mpiops.array_split(ifg_paths)
    ref_phs = np.array([_inner(ifg) for ifg in this_process_ifgs])
    log.info('Ref phase computed in process {}'.format(mpiops.rank))
    return ref_phs


def process_ifgs(ifg_paths, params, rows, cols):
    """
    Top level function to perform PyRate correction steps on given ifgs

    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    rows: int
        number of rows to break each ifg into
    cols: int
        number of cols to break each ifg into
    """
    if mpiops.size > 1:
        params[cf.PARALLEL] = False

    tiles = mpiops.run_once(get_tiles, ifg_paths[0], rows, cols)
    preread_ifgs = create_ifg_dict(ifg_paths,
                                   params=params,
                                   tiles=tiles)

    mst_calc(ifg_paths, params, tiles, preread_ifgs)

    # Estimate reference pixel location
    refpx, refpy = ref_pixel_calc(ifg_paths, params)

    # remove APS delay here, and write aps delay removed ifgs to disc
    # TODO: fix PyAPS integration
    if PyAPS_INSTALLED and \
            aps_delay_required(ifg_paths, params):  # pragma: no cover
        # ifgs = aps.remove_aps_delay(ifg_paths, params)
        log.info('Finished APS delay correction')
        # make sure aps correction flags are consistent
        if params[cf.APS_CORRECTION]:
            check_aps_ifgs(ifg_paths)

    # Estimate and remove orbit errors
    orb_fit_calc(ifg_paths, params, preread_ifgs)

    # calc and remove reference phase
    ref_phase_estimation(ifg_paths, params, refpx, refpy)

    maxvar, vcmt = maxvar_vcm_calc(ifg_paths, params, preread_ifgs)
    save_numpy_phase(ifg_paths, tiles, params)

    if params[cf.TIME_SERIES_CAL]:
        timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    # Calculate linear rate map
    linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    log.info('PyRate workflow completed')
    return (refpx, refpy), maxvar, vcmt


def linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    """
    mpi capable linrate calculation

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    vcmt: ndarrray
        vcmt array
    tiles: list
        list of all tiles used during mpi processes
    preread_ifgs: dict
        dict containing ifg characteristics for efficient computing
    """

    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating linear rate')
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        log.info('calculating lin rate of tile {}'.format(t.index))
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


def maxvar_vcm_calc(ifg_paths, params, preread_ifgs):
    """
    mpi capable maxvar and vcmt computation

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    preread_ifgs: dict
        dict containing ifg characteristics for efficient computing

    Returns
    -------
    maxvar: ndarray
        array of shape (nifgs, 1)
    vcmt: ndarray
        array of shape (nifgs, nifgs)
    """
    log.info('Calculating maxvar and vcm')
    process_indices = mpiops.array_split(range(len(ifg_paths)))
    prcs_ifgs = mpiops.array_split(ifg_paths)
    process_maxvar = []
    for n, i in enumerate(prcs_ifgs):
        log.info('Calculating maxvar for {} of process ifgs {} of '
                 'total {}'.format(n+1, len(prcs_ifgs), len(ifg_paths)))
        # TODO: cvd calculation is still pretty slow - revisit
        process_maxvar.append(vcm_module.cvd(i, params)[0])
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


def phase_sum(ifg_paths, params):
    """
    save phase data and phs_sum used in the reference phase estimation
    Parameters
    ----------
    ifg_paths: list:
        list of paths to ifgs
    params: dict
        config dict
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


def timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    """
    Time series calculation.

    Parameters
    ----------
    ifg_paths: list
        list of ifg paths
    params: dict
        parameters dict corresponding to config file
    vcmt: ndarrray
        vcmt array
    tiles: list
        list of all tiles used during mpi processes
    preread_ifgs: dict
        dict containing ifg characteristics for efficient computing

    """
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating time series')
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        log.info('Calculating time series for tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_tile = np.load(os.path.join(output_dir,
                                        'mst_mat_{}.npy'.format(t.index)))
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, _ = res
        np.save(file=os.path.join(output_dir, 'tsincr_{}.npy'.format(t.index)),
                arr=tsincr)
        np.save(file=os.path.join(output_dir, 'tscuml_{}.npy'.format(t.index)),
                arr=tscum)
    mpiops.comm.barrier()


def main(config_file, rows, cols):  # pragma: no cover
    """ linear rate and timeseries execution starts here """
    _, dest_paths, pars = cf.get_ifg_paths(config_file)
    process_ifgs(sorted(dest_paths), pars, rows, cols)
