"""
Main workflow script for PyRate
"""
from __future__ import print_function

import logging
import os
import pickle as cp
from os.path import join
import numpy as np

from pyrate import algorithm
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import linrate
from pyrate import mpiops
from pyrate import mst
from pyrate import orbital
from pyrate import prepifg
from pyrate import ref_phs_est as rpe
from pyrate import refpixel
from pyrate import shared
from pyrate import timeseries
from pyrate import vcm as vcm_module
from pyrate.compat import PyAPS_INSTALLED
from pyrate.shared import Ifg, write_output_geotiff, create_tiles, \
    PrereadIfg, prepare_ifg, save_numpy_phase, get_projection_info

if PyAPS_INSTALLED:
    from pyrate import aps
    from pyrate.aps import check_aps_ifgs, aps_delay_required

MASTER_PROCESS = 0
log = logging.getLogger(__name__)


def get_tiles(ifg_path, rows, cols):
    ifg = Ifg(ifg_path)
    ifg.open(readonly=True)
    tiles = create_tiles(ifg.shape, nrows=rows, ncols=cols)
    ifg.close()
    return tiles


def _join_dicts(dicts):
    if dicts is None:
        return
    d = {k: v for D in dicts for k, v in D.items()}
    return d


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
    ifgs_dict: dict
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
                                  ncols=ifg.ncols)
        ifg.close()
    ifgs_dict = _join_dicts(mpiops.comm.allgather(ifgs_dict))

    if mpiops.rank == MASTER_PROCESS:
        preread_ifgs = join(params[cf.OUT_DIR], 'preread_ifgs.pk')
        # add some extra information that's also useful later
        gt, md, wkt = get_projection_info(process_tifs[0])
        ifgs_dict['epochlist'] = algorithm.get_epochs(ifgs_dict)[0]
        ifgs_dict['gt'] = gt
        ifgs_dict['md'] = md
        ifgs_dict['wkt'] = wkt
        cp.dump(ifgs_dict, open(preread_ifgs, 'wb'))

    log.info('finish converting phase_data to numpy '
             'in process {}'.format(mpiops.rank))
    return ifgs_dict


def mst_calc(dest_tifs, params, tiles, preread_ifgs):
    """
    MPI function that control each process during MPI run
    """

    log.info('Calculating mst')
    log.info('Calculating minimum spanning tree matrix '
             'using NetworkX method')
    process_tiles = mpiops.array_split(tiles)

    def save_mst_tile(tile, i, preread_ifgs):
        if params[cf.NETWORKX_OR_MATLAB_FLAG]:
            mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        else:
            raise cf.ConfigException('Matlab mst not supported yet')
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

    # unlikely, but possible the refpixel can be (0,0)
    # check if there is a pre-specified reference pixel coord
    log.info('Starting reference pixel calculation')
    refx = params[cf.REFX]
    ifg = Ifg(ifg_paths[0])
    ifg.open(readonly=True)
    if refx > ifg.ncols - 1:
        raise ValueError("Invalid reference pixel X coordinate: %s" % refx)

    refy = params[cf.REFY]
    if refy > ifg.nrows - 1:
        raise ValueError("Invalid reference pixel Y coordinate: %s" % refy)

    if refx == 0 or refy == 0:  # matlab equivalent
        refy, refx = ref_pixel_mpi(ifg_paths, params)
        log.info('Found reference pixel coordinate: (%s, %s)' % (refx, refy))
    else:
        log.info('Reusing config file reference pixel (%s, %s)' % (refx, refy))
    ifg.close()
    return refx, refy


def ref_pixel_mpi(ifg_paths, params):
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
    log.info('Saving ref pixel blocks')
    outdir = params[cf.OUT_DIR]
    for p in ifg_paths:
        ifg = Ifg(p)
        ifg.open(readonly=True)
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.convert_to_mm()
        for y, x in grid:
            data = ifg.phase_data[y - half_patch_size:y + half_patch_size + 1,
                                  x - half_patch_size:x + half_patch_size + 1]

            data_file = join(outdir, 'ref_phase_data_{b}_{y}_{x}.npy'.format(
                b=os.path.basename(p).split('.')[0], y=y, x=x))
            np.save(file=data_file, arr=data)
        ifg.close()
    log.info('Saved ref pixel blocks')


def orb_fit_calc(ifg_paths, params):
    log.info('Calculating orbfit correction')
    if params[cf.ORBITAL_FIT_METHOD] != 1:
        raise cf.ConfigException('Only orbfit method 1 is supported')
    process_ifgs = mpiops.array_split(ifg_paths)
    mlooked = None
    # TODO: MPI orbfit method 2
    orbital.orbital_correction(process_ifgs, params, mlooked=mlooked)
    log.info('Finished orbfit calculation in process {}'.format(mpiops.rank))


def ref_phase_estimation(ifg_paths, params, refpx, refpy):
    # TODO: may benefit from tiling and using a median of median algorithms
    log.info('Estimating and removing reference phase')
    if params[cf.REF_EST_METHOD] == 1:
        # calculate phase sum for later use in ref phase method 1
        comp = phase_sum_mpi(ifg_paths, params)
        process_ref_phs = ref_phs_method1(ifg_paths, comp)
    elif params[cf.REF_EST_METHOD] == 2:
        process_ref_phs = ref_phs_method2(ifg_paths, params, refpx, refpy)
    else:
        raise cf.ConfigException('Ref phase estimation method must be 1 or 2')

    ref_phs_file = join(params[cf.OUT_DIR], 'ref_phs.npy')
    if mpiops.rank == MASTER_PROCESS:
        ref_phs = np.zeros(len(ifg_paths), dtype=np.float64)
        process_indices = mpiops.array_split(range(len(ifg_paths)))
        ref_phs[process_indices] = process_ref_phs
        for r in range(1, mpiops.size):
            process_indices = mpiops.array_split(range(len(ifg_paths)), r)
            this_process_ref_phs = np.zeros(shape=len(process_indices),
                                            dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=r, tag=r)
            ref_phs[process_indices] = this_process_ref_phs
        np.save(file=ref_phs_file, arr=ref_phs)
    else:
        # send reference phase data to master process
        mpiops.comm.Send(process_ref_phs, dest=MASTER_PROCESS,
                         tag=mpiops.rank)


def ref_phs_method2(ifg_paths, params, refpx, refpy):
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

    def _inner(ifg_path):
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
    ifgs: sequence of paths to interferrograms
    (NB: changes are saved into ifgs)
    params: dictionary of configuration parameters
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
    if PyAPS_INSTALLED and aps_delay_required(ifg_paths, params):
        ifgs = aps.remove_aps_delay(ifg_paths, params)
        log.info('Finished APS delay correction')
        # make sure aps correction flags are consistent
        if params[cf.APS_CORRECTION]:
            check_aps_ifgs(ifg_paths)

    # Estimate and remove orbit errors
    orb_fit_calc(ifg_paths, params)

    # calc and remove reference phase
    ref_phase_estimation(ifg_paths, params, refpx, refpy)

    maxvar, vcmt = maxvar_vcm_mpi(ifg_paths, params, preread_ifgs)
    save_numpy_phase(ifg_paths, tiles, params)

    if params[cf.TIME_SERIES_CAL]:
        timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    # Calculate linear rate map
    linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs)

    log.info('PyRate workflow completed')
    return (refpx, refpy), maxvar, vcmt


def linrate_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):

    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating linear rate')
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        i = t.index
        log.info('calculating lin rate of tile {}'.format(i))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_n = os.path.join(output_dir, 'mst_mat_{}.npy'.format(i))
        mst_grid_n = np.load(mst_n)
        res = linrate.linear_rate(ifg_parts, params, vcmt, mst_grid_n)

        for r in res:
            if r is None:
                raise ValueError('TODO: bad value')
        rate, error, samples = res
        # declare file names
        rate_file = os.path.join(output_dir, 'linrate_{}.npy'.format(i))
        error_file = os.path.join(output_dir, 'linerror_{}.npy'.format(i))
        samples_file = os.path.join(output_dir, 'linsamples_{}.npy'.format(i))
        np.save(file=rate_file, arr=rate)
        np.save(file=error_file, arr=error)
        np.save(file=samples_file, arr=samples)
    mpiops.comm.barrier()


def maxvar_vcm_mpi(ifg_paths, params, preread_ifgs):
    log.info('Calculating maxvar and vcm')
    process_indices = mpiops.array_split(range(len(ifg_paths)))
    process_ifgs = mpiops.array_split(ifg_paths)
    process_maxvar = []
    for n, i in enumerate(process_ifgs):
        log.info('Calculating maxvar for {} of process ifgs {} of '
                 'total {}'.format(n+1, len(process_ifgs), len(ifg_paths)))
        # TODO: cvd calculation is still pretty slow - revisit
        process_maxvar.append(vcm_module.cvd(i, params)[0])
    if mpiops.rank == MASTER_PROCESS:
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        maxvar[process_indices] = process_maxvar
        for i in range(1, mpiops.size):
            rank_indices = mpiops.array_split(range(len(ifg_paths)), i)
            this_process_ref_phs = np.empty(len(rank_indices),
                                            dtype=np.float64)
            mpiops.comm.Recv(this_process_ref_phs, source=i, tag=i)
            maxvar[rank_indices] = this_process_ref_phs
    else:
        maxvar = np.empty(len(ifg_paths), dtype=np.float64)
        mpiops.comm.Send(np.array(process_maxvar, dtype=np.float64),
                         dest=MASTER_PROCESS, tag=mpiops.rank)

    maxvar = mpiops.comm.bcast(maxvar, root=0)
    vcmt = mpiops.run_once(vcm_module.get_vcmt, preread_ifgs, maxvar)
    return maxvar, vcmt


def phase_sum_mpi(ifg_paths, params):
    """
    save phase data and phase_sum used in the reference phase estimation
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
    phase_sum = np.zeros(shape=shape, dtype=np.float64)
    ifg.close()

    for d in p_paths:
        ifg = Ifg(d)
        ifg.open()
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        phase_sum += ifg.phase_data
        ifg.close()

    if mpiops.rank == MASTER_PROCESS:
        phase_sum_all = phase_sum
        for i in range(1, mpiops.size):  # loop is better for memory
            phase_sum = np.zeros(shape=shape, dtype=np.float64)
            mpiops.comm.Recv(phase_sum, source=i, tag=i)
            phase_sum_all += phase_sum
        comp = np.isnan(phase_sum_all)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
    else:
        comp = None
        mpiops.comm.Send(phase_sum, dest=0, tag=mpiops.rank)

    comp = mpiops.comm.bcast(comp, root=0)
    return comp


def timeseries_calc(ifg_paths, params, vcmt, tiles, preread_ifgs):
    process_tiles = mpiops.array_split(tiles)
    log.info('Calculating time series')
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        i = t.index
        log.info('Calculating time series for tile {}'.format(t.index))
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_file_process_n = os.path.join(output_dir,
                                          'mst_mat_{}.npy'.format(i))
        mst_tile = np.load(mst_file_process_n)
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, _ = res
        tsincr_file = os.path.join(output_dir, 'tsincr_{}.npy'.format(i))
        tscum_file = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))
        np.save(file=tsincr_file, arr=tsincr)
        np.save(file=tscum_file, arr=tscum)
    mpiops.comm.barrier()


def remove_orbital_error(ifgs, params):
    log.info('Calculating orbital error correction')

    if not params[cf.ORBITAL_FIT]:
        log.info('Orbital correction not required')
        return

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_ORBITAL_ERROR) for i in ifgs]

    if all(flags):
        log.info('Skipped orbital correction, ifgs already corrected')
        return
    else:
        check_orbital_ifgs(ifgs, flags)

    mlooked = None

    if (params[cf.ORBITAL_FIT_LOOKS_X] > 1 or
            params[cf.ORBITAL_FIT_LOOKS_Y] > 1):
        # resampling here to use all prior corrections to orig data
        # can't do multiprocessing without writing to disc, but can do MPI
        # due to swig pickling issue. So multiprocesing is not implemented
        mlooked_dataset = prepifg.prepare_ifgs(
            [i.data_path for i in ifgs],
            crop_opt=prepifg.ALREADY_SAME_SIZE,
            xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
            ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
            thresh=params[cf.NO_DATA_AVERAGING_THRESHOLD],
            write_to_disc=False)
        mlooked = [Ifg(m[1]) for m in mlooked_dataset]

        for m in mlooked:
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]

    orbital.orbital_correction(ifgs, params, mlooked=mlooked)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == ifc.ORB_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        log.debug('Detected mix of corrected and uncorrected '
                  'orbital error in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '{}: prior orbital error correction detected'.format(i)
            else:
                msg = '{}: no orbital correction detected'.format(i)
            log.debug(msg)
            raise orbital.OrbitalError(msg)


def calculate_linear_rate(ifgs, params, vcmt, mst_mat=None):
    log.info('Calculating linear rate')
    res = linrate.linear_rate(ifgs, params, vcmt, mst_mat)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res
    write_linrate_tifs(ifgs, params, res)
    log.info('Linear rate calculated')
    return rate, error, samples


def write_linrate_tifs(ifgs, params, res):
    log.info('Writing linrate results')
    rate, error, samples = res
    gt, md, wkt = get_projection_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]
    dest = join(params[cf.OUT_DIR], "linrate.tif")
    md[ifc.MASTER_DATE] = epochlist.dates
    md[ifc.PRTYPE] = 'linrate'
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "linerror.tif")
    md[ifc.PRTYPE] = 'linerror'
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "linsamples.tif")
    md[ifc.PRTYPE] = 'linsamples'
    write_output_geotiff(md, gt, wkt, samples, dest, np.nan)
    write_linrate_numpy_files(error, rate, samples, params)


def write_linrate_numpy_files(error, rate, samples, params):
    rate_file = join(params[cf.OUT_DIR], 'rate.npy')
    error_file = join(params[cf.OUT_DIR], 'error.npy')
    samples_file = join(params[cf.OUT_DIR], 'samples.npy')
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)


def warp_required(xlooks, ylooks, crop):
    """
    Returns True if params show rasters need to be cropped and/or resized.
    """

    if xlooks > 1 or ylooks > 1:
        return True

    if crop is None:
        return False

    return True


def main(config_file, rows, cols):
    """ linear rate and timeseries execution starts here """
    _, dest_paths, pars = cf.get_ifg_paths(config_file)
    process_ifgs(sorted(dest_paths), pars, rows, cols)
