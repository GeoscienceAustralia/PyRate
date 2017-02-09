"""
Main workflow script for PyRate
"""
from __future__ import print_function
import logging
import os
from os.path import join
import numpy as np
from osgeo import gdal
import pickle as cp

from pyrate import config as cf
from pyrate import linrate
from pyrate import mst
from pyrate import orbital
from pyrate import prepifg
from pyrate import refpixel
from pyrate import timeseries
from pyrate import shared
from pyrate import algorithm
from pyrate import ifgconstants as ifc
from pyrate import matlab_mst as matlab_mst
from pyrate import ref_phs_est as rpe
from pyrate import vcm as vcm_module
from pyrate.shared import Ifg, write_output_geotiff, \
    pre_prepare_ifgs, create_tiles, PrereadIfg, prepare_ifg, save_numpy_phase
from pyrate.compat import PyAPS_INSTALLED
from pyrate import mpiops
if PyAPS_INSTALLED:
    from pyrate import remove_aps_delay as aps

PYRATEPATH = cf.PYRATEPATH
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
        gt, md, wkt = get_projection_info(process_tifs[0], params)
        ifgs_dict['epochlist'] = algorithm.get_epochs(ifgs_dict)
        ifgs_dict['gt'] = gt
        ifgs_dict['md'] = md
        ifgs_dict['wkt'] = wkt
        cp.dump(ifgs_dict, open(preread_ifgs, 'wb'))

    log.info('finish converting phase_data to numpy '
             'in process {}'.format(mpiops.rank))
    return ifgs_dict


def mpi_mst_calc(dest_tifs, params, tiles, preread_ifgs):
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


def temp_mst_grid_reconstruct(tiles, ifg_paths, params):
    ifg = Ifg(ifg_paths[0])
    ifg.open(readonly=True)

    mst_grid = np.empty(shape=(len(ifg_paths), ifg.shape[0], ifg.shape[1]),
                        dtype=bool)
    ifg.close()
    for t in tiles:
        mst_f = join(params[cf.OUT_DIR], 'mst_mat_{}.npy'.format(t.index))
        mst_grid[:, t.top_left_y:t.bottom_right_y,
                 t.top_left_x:t.bottom_right_x] = np.load(mst_f)
    return mst_grid


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
        refx, refy = ref_pixel_mpi(ifg_paths, params)
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
    refx, refy = mpiops.run_once(refpixel.filter_means, mean_sds, grid)
    return refx, refy


def save_ref_pixel_blocks(grid, half_patch_size, ifg_paths, params):
    outdir = params[cf.OUT_DIR]
    for p in ifg_paths:
        for y, x in grid:
            ifg = Ifg(p)
            ifg.open(readonly=True)
            ifg.nodata_value = params[cf.NO_DATA_VALUE]
            ifg.convert_to_nans()
            ifg.convert_to_mm()
            data = ifg.phase_data[y - half_patch_size:y + half_patch_size + 1,
                                  x - half_patch_size:x + half_patch_size + 1]

            data_file = join(outdir, 'ref_phase_data_{b}_{y}_{x}.npy'.format(
                                         b=os.path.basename(p).split('.')[0],
                                         y=y, x=x))
            np.save(file=data_file, arr=data)


def orb_fit_calc(ifg_paths, params):
    log.info('Calculating orbfit correction')
    if params[cf.ORBITAL_FIT_METHOD] != 1:
        raise cf.ConfigException('Only orbfit method 1 is supported')
    process_ifgs = mpiops.array_split(ifg_paths)
    mlooked = None
    # TODO: MPI orbfit method 2
    orbital.orbital_correction(process_ifgs, params, mlooked=mlooked)
    log.info('Finished orbfit calculation in process {}'.format(mpiops.rank))


def ref_phase_estimation_mpi(ifg_paths, params, refpx, refpy):
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
        ref_ph = rpe.est_ref_phase_method2_multi(phase_data,
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
        ref_phase = rpe.est_ref_phase_method1_multi(phase_data, comp)
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

    mpi_mst_calc(ifg_paths, params, tiles, preread_ifgs)
    mst_grid = temp_mst_grid_reconstruct(tiles, ifg_paths, params)

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
    ref_phase_estimation_mpi(ifg_paths, params, refpx, refpy)

    maxvar, vcmt = maxvar_vcm_mpi(ifg_paths, params, preread_ifgs)
    ifgs = pre_prepare_ifgs(ifg_paths, params)
    save_numpy_phase(ifg_paths, tiles, params)

    if params[cf.TIME_SERIES_CAL] != 0:
        time_series_mpi(ifg_paths, params, vcmt, tiles, preread_ifgs)
        # compute_time_series(ifgs, mst_grid, params, vcmt)

    # Calculate linear rate map
    rate, error, samples = calculate_linear_rate(ifgs, params, vcmt, mst_grid)
    linrate_mpi(ifg_paths, params, vcmt, tiles, preread_ifgs)

    log.info('PyRate workflow completed')
    return mst_grid, (refpx, refpy), maxvar, vcmt, rate, error, samples


def linrate_mpi(ifg_paths, params, vcmt, tiles, preread_ifgs):

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
            this_process_ref_phs = np.empty(len(rank_indices), dtype=np.float64)
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


def aps_delay_required(ifgs, params):
    log.info('Removing APS delay')

    if not params[cf.APS_CORRECTION]:
        log.info('APS delay removal not required')
        return False

    # perform some general error/sanity checks
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]

    if all(flags):
        log.info('Skipped APS delay removal, ifgs are already aps corrected')
        return False
    else:
        check_aps_ifgs(ifgs)

    return True


def check_aps_ifgs(ifgs):
    flags = [i.dataset.GetMetadataItem(ifc.PYRATE_APS_ERROR) for i in ifgs]
    count = sum([f == aps.APS_STATUS for f in flags])
    if (count < len(flags)) and (count > 0):
        log.debug('Detected mix of corrected and uncorrected '
                      'APS delay in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '%s: prior APS delay correction detected'
            else:
                msg = '%s: no APS delay correction detected'
            logging.debug(msg % i.data_path)
        raise aps.APSException('Mixed APS removal status in ifg list')


def mst_calculation(ifg_paths_or_instance, params):
    if isinstance(ifg_paths_or_instance, list):
        ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)
        log.info('Calculating minimum spanning tree matrix using '
                 'NetworkX')

        mst_grid = mst.mst_parallel(ifgs, params)
    else:
        nan_conversion = params[cf.NAN_CONVERSION]
        assert isinstance(ifg_paths_or_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_paths_or_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_mm()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_paths_or_instance,
                               nodata_value=params[cf.NO_DATA_VALUE],
                               nan_conversion=nan_conversion)
        log.info('Calculating minimum spanning tree matrix '
                 'using Matlab-algorithm method')
        mst_grid = matlab_mst.matlab_mst_boolean_array(ifg_instance_updated)

    # write mst output to a file
    mst_mat_binary_file = join(
        PYRATEPATH, params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    for i in ifgs:
        i.close()
    return mst_grid


def calculate_time_series(ifgs, params, vcmt, mst):
    res = timeseries.time_series(ifgs, params, vcmt, mst)
    for r in res:
        if len(r.shape) != 3:
            logging.error('TODO: time series result shape is incorrect')
            raise timeseries.TimeSeriesError

    log.info('Time series calculated')
    tsincr, tscum, tsvel = res
    return tsincr, tscum, tsvel


def time_series_mpi(ifg_paths, params, vcmt, tiles, preread_ifgs):
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
        tsincr, tscum, tsvel = res
        tsincr_file = os.path.join(output_dir, 'tsincr_{}.npy'.format(i))
        tscum_file = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))
        np.save(file=tsincr_file, arr=tsincr)
        np.save(file=tscum_file, arr=tscum)
    mpiops.comm.barrier()


def compute_time_series(ifgs, mst_grid, params, vcmt):
    log.info('Calculating time series')
    # Calculate time series
    tsincr, tscum, tsvel = calculate_time_series(
        ifgs, params, vcmt=vcmt, mst=mst_grid)

    # tsvel_file = join(params[cf.OUT_DIR], 'tsvel.npy')
    tsincr_file = join(params[cf.OUT_DIR], 'tsincr.npy')
    tscum_file = join(params[cf.OUT_DIR], 'tscum.npy')
    np.save(file=tsincr_file, arr=tsincr)
    np.save(file=tscum_file, arr=tscum)
    # np.save(file=tsvel_file, arr=tsvel)

    # TODO: write tests for these functions
    write_timeseries_geotiff(ifgs, params, tsincr, pr_type='tsincr')
    write_timeseries_geotiff(ifgs, params, tscum, pr_type='tscuml')
    # write_timeseries_geotiff(ifgs, params, tsvel, pr_type='tsvel')
    return tsincr, tscum, tsvel


def get_projection_info(ifg_path, params):
    p = join(params[cf.OUT_DIR], ifg_path)
    ds = gdal.Open(p)
    md = ds.GetMetadata()  # get metadata for writing on output tifs
    gt = ds.GetGeoTransform()  # get geographical bounds of data
    wkt = ds.GetProjection()  # get projection of data
    ds = None  # close dataset
    return gt, md, wkt


def write_timeseries_geotiff(ifgs, params, tsincr, pr_type):
    # setup metadata for writing into result files
    gt, md, wkt = get_projection_info(ifgs, params)
    epochlist = algorithm.get_epochs(ifgs)

    for i in range(tsincr.shape[2]):
        md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
        md['PR_SEQ_POS'] = i  # sequence position

        data = tsincr[:, :, i]
        dest = join(PYRATEPATH, params[cf.OUT_DIR], pr_type + "_" +
                    str(epochlist.dates[i + 1]) + ".tif")
        md[ifc.PRTYPE] = pr_type
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)


# This is not used anywhere now, but may be useful
def insert_time_series_interpolation(ifg_instance_updated, params):

    edges = matlab_mst.get_sub_structure(ifg_instance_updated,
                                         np.zeros(len(ifg_instance_updated.id),
                                                  dtype=bool))

    _, _, ntrees = matlab_mst.matlab_mst_kruskal(edges, ntrees=True)
    # if ntrees=1, no interpolation; otherwise interpolate
    if ntrees > 1:
        params[cf.TIME_SERIES_INTERP] = 1
    else:
        params[cf.TIME_SERIES_INTERP] = 0

    return params


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

    orbital.orbital_correction(ifgs,
                               params,
                               mlooked=mlooked)


def check_orbital_ifgs(ifgs, flags):
    count = sum([f == ifc.ORB_REMOVED for f in flags])
    if (count < len(flags)) and (count > 0):
        logging.debug('Detected mix of corrected and uncorrected '
                      'orbital error in ifgs')

        for i, flag in zip(ifgs, flags):
            if flag:
                msg = '{}: prior orbital error correction detected'.format(i)
            else:
                msg = '{}: no orbital correction detected'.format(i)
            logging.debug(msg % i.data_path)

        raise orbital.OrbitalError(msg)


def calculate_linear_rate(ifgs, params, vcmt, mst=None):
    log.info('Calculating linear rate')
    res = linrate.linear_rate(ifgs, params, vcmt, mst)
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
    gt, md, wkt = get_projection_info(ifgs[0].data_path, params)
    epochlist = algorithm.get_epochs(ifgs)
    dest = join(PYRATEPATH, params[cf.OUT_DIR], "linrate.tif")
    md[ifc.MASTER_DATE] = epochlist.dates
    md[ifc.PRTYPE] = 'linrate'
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = join(PYRATEPATH, params[cf.OUT_DIR], "linerror.tif")
    md[ifc.PRTYPE] = 'linerror'
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)
    dest = join(PYRATEPATH, params[cf.OUT_DIR], "linsamples.tif")
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


# general function template
#
# add check for pre-existing metadata flag / skip if required
# perform calculation
# optionally save modified data to disk if required
# optionally save correction component to disk (more useful for debugging)
# set flag in dataset for correction
# write to log file


def warp_required(xlooks, ylooks, crop):
    """
    Returns True if params show rasters need to be cropped and/or resized.
    """

    if xlooks > 1 or ylooks > 1:
        return True

    if crop is None:
        return False

    return True


def working_ifg_paths(src_paths, xlooks, ylooks, cropping):
    """
    Filter. Returns paths to ifgs to process (eg. checks for mlooked tifs)
    """
    if warp_required(xlooks, ylooks, cropping):
        mlooked_unw = [cf.mlooked_path(p, xlooks, crop_out=cropping)
                       for p in src_paths]
        mlooked_paths = [os.path.splitext(m)[0]+'.tif' for m in mlooked_unw]

        if not all([os.path.exists(p) for p in mlooked_paths]):
            msg = 'Multilooked ifgs do not exist ' \
                  '(execute "run_prepifg.py" first)'
            raise IOError(msg)

        log.info('Using mlooked interferograms...')
        return mlooked_paths
    return src_paths  # multi looking not specified, work with base ifgs


def dest_ifg_paths(ifg_paths, outdir):
    """
    Returns paths to out/dest ifgs.
    """

    bases = [os.path.basename(p) for p in ifg_paths]
    return [join(outdir, p) for p in bases]


def main(config_file, rows, cols):
    base_unw_paths, dest_paths, pars = cf.get_ifg_paths(config_file)
    process_ifgs(sorted(dest_paths), pars, rows, cols)
