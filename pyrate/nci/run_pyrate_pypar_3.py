import datetime
import glob
import os
import sys
import gc
from operator import itemgetter
import psutil
import numpy as np
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import linrate
from pyrate import mst
from pyrate import orbital
from pyrate import reference_phase_estimation as rpe
from pyrate import refpixel
from pyrate import remove_aps_delay as aps
from pyrate import shared
from pyrate import timeseries
from pyrate import vcm as vcm_module
from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate.scripts.run_pyrate import write_msg
from pyrate.shared import get_tmpdir
from collections import namedtuple
import cPickle as cp
from osgeo import gdal
gdal.SetCacheMax(64)

TMPDIR = get_tmpdir()

__author__ = 'sudipta'

# Constants
MASTER_PROCESS = 0
data_path = 'DATAPATH'

PrereadIfg = namedtuple('PrereadIfg', 'path nan_fraction master slave time_span')


def main(params, config_file=sys.argv[1]):

    # setup paths
    xlks, ylks, crop = run_pyrate.transform_params(params)
    base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])
    dest_tifs = run_pyrate.get_dest_paths(base_unw_paths, crop, params, xlks)

    # Setting up parallelisation
    parallel = Parallel(True)
    rank = parallel.rank

    # calculate process information
    ifg_shape, process_tiles, process_indices, tiles = \
        get_process_tiles(dest_tifs, parallel, params)

    output_dir = params[cf.OUT_DIR]
    preread_ifgs = os.path.join(output_dir, 'preread_ifgs.pk')

    parallel.barrier()

    # all processes need access to maxvar, and vcmt
    maxvar, vcmt = maxvar_vcm_mpi(rank, dest_tifs, parallel, params)

    parallel.barrier()

    if rank == MASTER_PROCESS:
        for d in dest_tifs:
            save_latest_phase(d, output_dir, tiles)

    # linrate mpi computation
    linrate_mpi(rank, dest_tifs, parallel, params, vcmt,
                process_tiles, process_indices, ifg_shape, preread_ifgs)
    # TODO: write linrate aggregation function from linrate tiles

    parallel.barrier()
    # time series mpi computation
    if params[cf.TIME_SERIES_CAL]:
        time_series_mpi(dest_tifs, params, vcmt, process_tiles,
                        process_indices, preread_ifgs)
        parallel.barrier()
        write_time_series_geotiff_mpi(dest_tifs, params, tiles, parallel, rank)

    parallel.finalize()


def save_latest_phase(d, output_dir, tiles):
    ifg = shared.Ifg(d)
    ifg.open()
    ifg.nodata_value = 0
    phase_data = ifg.phase_data
    for t in tiles:
        p_data = phase_data[t.top_left_y:t.bottom_right_y,
                 t.top_left_x:t.bottom_right_x]
        phase_file = 'phase_data_{}_{}.npy'.format(
            os.path.basename(d).split('.')[0], t.index)

        np.save(file=os.path.join(output_dir, phase_file),
                arr=p_data)
    return ifg


def write_time_series_geotiff_mpi(dest_tifs, params, tiles, parallel, MPI_id):
    ifgs = shared.prepare_ifgs_without_phase(dest_tifs, params)
    epochlist, gt, md, wkt = run_pyrate.setup_metadata(ifgs, params)

    # load the first tsincr file to determine the number of time series tifs
    tsincr_file = os.path.join(TMPDIR, 'tsincr_0.npy')
    tsincr = np.load(file=tsincr_file)

    no_ts_tifs = tsincr.shape[2]
    process_tifs = parallel.calc_indices(no_ts_tifs)

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop is helps us over come the memory limit
    print 'process {} will write {} ts tifs of total {}'.format(
        MPI_id, len(process_tifs), no_ts_tifs)

    for i in process_tifs:
        tsincr_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        tscum_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        for n, t in enumerate(tiles):
            tsincr_file = os.path.join(TMPDIR, 'tsincr_{}.npy'.format(n))
            tscum_file = os.path.join(TMPDIR, 'tscuml_{}.npy'.format(n))
            tsincr = np.load(file=tsincr_file)
            tscum = np.load(file=tscum_file)

            md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
            md['PR_SEQ_POS'] = i  # sequence position
            tsincr_g[t.top_left_y:t.bottom_right_y,
                     t.top_left_x:t.bottom_right_x] = tsincr[:, :, i]
            dest = os.path.join(params[cf.OUT_DIR],
                'tsincr' + "_" + str(epochlist.dates[i + 1]) + ".tif")
            md[ifc.PRTYPE] = 'tsincr'
            shared.write_output_geotiff(md, gt, wkt, tsincr_g, dest, np.nan)

            tscum_g[t.top_left_y:t.bottom_right_y,
                t.top_left_x:t.bottom_right_x] = tscum[:, :, i]
            dest = os.path.join(params[cf.OUT_DIR],
                'tscuml' + "_" + str(epochlist.dates[i + 1]) + ".tif")
            md[ifc.PRTYPE] = 'tscuml'
            shared.write_output_geotiff(md, gt, wkt, tscum_g, dest, np.nan)


def linrate_mpi(MPI_myID, ifg_paths, parallel, params, vcmt,
                process_tiles, process_indices, ifg_shape, preread_ifgs):
    write_msg('Calculating linear rate')
    for i, t in zip(process_indices, process_tiles):
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_n = os.path.join(TMPDIR, 'mst_mat_{}.npy'.format(i))
        mst_grid_n = np.load(mst_n)
        res = linrate.linear_rate(ifg_parts, params, vcmt, mst_grid_n)

        for r in res:
            if r is None:
                raise ValueError('TODO: bad value')
        rate, error, samples = res
        # declare file names
        rate_file = os.path.join(params[cf.OUT_DIR], 'rate_{}.npy'.format(i))
        error_file = os.path.join(params[cf.OUT_DIR], 'error_{}.npy'.format(i))
        samples_file = os.path.join(params[cf.OUT_DIR],
                                    'samples_{}.npy'.format(i))

        np.save(file=rate_file, arr=rate)
        np.save(file=error_file, arr=rate)
        np.save(file=samples_file, arr=samples)


def time_series_mpi(ifg_paths, params, vcmt, process_tiles, process_indices,
                    preread_ifgs):
    write_msg('Calculating time series')  # this should be logged

    for t in process_tiles:
        i = t.index
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_file_process_n = os.path.join(TMPDIR, 'mst_mat_{}.npy'.format(i))
        mst_tile = np.load(mst_file_process_n)
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, tsvel = res
        tsincr_file = os.path.join(TMPDIR, 'tsincr_{}.npy'.format(i))
        tscum_file = os.path.join(TMPDIR, 'tscuml_{}.npy'.format(i))
        np.save(file=tsincr_file, arr=tsincr)
        np.save(file=tscum_file, arr=tscum)


def ref_phase_estimation_mpi(MPI_myID, ifg_paths, parallel, params,
                             refpx, refpy):
    write_msg('Finding and removing reference phase')
    num_processors = parallel.size
    # ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifg_paths)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifg_paths) for p in process_indices]
    process_ref_phs = np.zeros(len(process_ifgs))
    output_dir = params[cf.OUT_DIR]

    if params[cf.REF_EST_METHOD] == 1:
        # TODO: revisit as this will likely hit memory limit in NCI
        for n, p in enumerate(process_ifgs):
            process_ref_phs[n] = ref_phase_method1_dummy(p, output_dir)

    elif params[cf.REF_EST_METHOD] == 2:
        half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
        chipsize = 2 * half_chip_size + 1
        thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
        for n, ifg in enumerate(process_ifgs):
            ref_phs = rpe.est_ref_phase_method2_multi(ifg.phase_data,
                                                      half_chip_size,
                                                      refpx, refpy, thresh)
            ifg.phase_data -= ref_phs
            ifg.meta_data[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
            ifg.write_modified_phase()
            process_ref_phs[n] = ref_phs
            ifg.close()
    else:
        raise cf.ConfigException('Ref phase estimation method must be 1 or 2')

    ref_phs_file = os.path.join(params[cf.OUT_DIR], 'ref_phs.npy')

    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        ref_phs = np.zeros(no_ifgs)
        ref_phs[process_indices] = process_ref_phs

        for i in range(1, num_processors):
            process_ref_phs = \
                parallel.receive(source=i, tag=-1,
                                 return_status=False)
            ref_phs[all_indices[i]] = process_ref_phs
        np.save(file=ref_phs_file, arr=ref_phs)
    else:
        # send reference phase data to master process
        parallel.send(process_ref_phs, destination=MASTER_PROCESS,
                      tag=MPI_myID)


def ref_phase_method1_dummy(ifg_path, output_dir):
    comp_file = os.path.join(output_dir, 'comp.npy')

    comp = np.load(comp_file)
    process = psutil.Process(os.getpid())
    print process.memory_info()
    print psutil.virtual_memory()
    numpy_file = os.path.join(
        output_dir, os.path.basename(ifg_path).split('.')[0] + '.npy')
    phase_data = np.load(numpy_file)
    ref_phs = rpe.est_ref_phase_method1_multi(phase_data, comp)
    phase_data -= ref_phs
    ifg = shared.Ifg(ifg_path)
    ifg.open()
    md = ifg.meta_data
    md[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
    ifg.write_modified_phase(data=phase_data)
    ifg.close()
    return ref_phs


def maxvar_vcm_mpi(rank, ifg_paths, parallel, params):
    num_processors = parallel.size
    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_maxvar = [vcm_module.cvd(i, params)[0] for i in process_ifgs]
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    vcmt_file = os.path.join(params[cf.OUT_DIR], 'vcmt.npy')
    if rank == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        maxvar = np.empty(len(ifgs))
        maxvar[all_indices[MASTER_PROCESS]] = process_maxvar

        for i in range(1, num_processors):
            process_maxvar = parallel.receive(source=i, tag=-1,
                                              return_status=False)
            maxvar[all_indices[i]] = process_maxvar
        np.save(file=maxvar_file, arr=maxvar)
    else:
        parallel.send(process_maxvar, destination=MASTER_PROCESS, tag=rank)
    parallel.barrier()
    maxvar = np.load(maxvar_file)
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)
    if rank == MASTER_PROCESS:
        np.save(file=vcmt_file, arr=vcmt)

    return maxvar, vcmt


def orb_fit_calc_mpi(MPI_myID, ifg_paths, num_processors, parallel, params):
    print 'calculating orbfit using MPI id:', MPI_myID
    if params[cf.ORBITAL_FIT_METHOD] != 1:
        raise cf.ConfigException('For now orbfit method must be 1')

    # ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifg_paths)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifg_paths) for p in process_indices]

    mlooked = None
    # difficult to enable MPI on orbfit method 2
    # so just do orbfit method 1 for now
    # TODO: MPI orbfit method 2
    orbital.orbital_correction(process_ifgs,
                               params,
                               mlooked=mlooked)
    # set orbfit tags after orbital error correction
    for i in process_ifgs:
        ifg = shared.Ifg(i)
        ifg.open()
        ifg.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR, ifc.ORB_REMOVED)
        ifg.write_modified_phase()
        ifg.close()
        # implement mpi logging


def ref_pixel_calc_mpi(MPI_myID, ifg_paths, num_processors, parallel, params):
    half_patch_size, thresh, grid = refpixel.ref_pixel_setup(ifg_paths, params)

    save_ref_pixel_blocks(grid, half_patch_size, ifg_paths, parallel, params)
    parallel.barrier()
    no_steps = len(grid)
    process_indices = parallel.calc_indices(no_steps)
    process_grid = [itemgetter(p)(grid) for p in process_indices]
    print 'Processor {mpi_id} has {processes} ' \
          'tiles out of {num_files}'.format(mpi_id=MPI_myID,
                                            processes=len(process_indices),
                                            num_files=no_steps)
    mean_sds = refpixel.ref_pixel_mpi(process_grid, half_patch_size,
                                      ifg_paths, thresh, params)
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_steps)
        mean_sds_final = np.empty(shape=no_steps)
        mean_sds_final[all_indices[MASTER_PROCESS]] = mean_sds
        for i in range(1, num_processors):
            process_mean_sds = parallel.receive(source=i, tag=-1,
                                                return_status=False)
            mean_sds_final[all_indices[i]] = process_mean_sds

        refx, refy = refpixel.filter_means(mean_sds_final, grid)
        print 'finished calculating ref pixel'
        return refx, refy
    else:
        parallel.send(mean_sds, destination=MASTER_PROCESS, tag=MPI_myID)
        print 'sent ref pixel to master'


def save_ref_pixel_blocks(grid, half_patch_size, ifg_paths, parallel, params):
    no_ifgs = len(ifg_paths)
    process_path_indices = parallel.calc_indices(no_ifgs)
    process_ifg_paths = [itemgetter(p)(ifg_paths) for p in process_path_indices]
    outdir = params[cf.OUT_DIR]
    for p in process_ifg_paths:
        for y, x in grid:
            ifg = shared.Ifg(p)
            ifg.open(readonly=True)
            ifg.nodata_value = params[cf.NO_DATA_VALUE]
            ifg.convert_to_nans()
            ifg.convert_to_mm()
            data = ifg.phase_data[y - half_patch_size:y + half_patch_size + 1,
                   x - half_patch_size:x + half_patch_size + 1]

            data_file = os.path.join(outdir,
                                     'ref_phase_data_{b}_{y}_{x}.npy'.format(
                                         b=os.path.basename(p).split('.')[0],
                                         y=y, x=x)
                                     )
            np.save(file=data_file, arr=data)


def mpi_mst_calc(dest_tifs, process_tiles, process_indices, preread_ifgs):
    """
    MPI function that control each process during MPI run
    :param MPI_myID:
    :param dest_tifs: paths of cropped amd resampled geotiffs
    :param mpi_log_filename:
    :param num_processors:
    :param parallel: MPI Parallel class instance
    :param params: config params dictionary
    :param mst_file: mst file (2d numpy array) save to disc
    :return:
    """
    write_msg('Calculating mst')

    write_msg('Calculating minimum spanning tree matrix '
                         'using NetworkX method')

    def save_mst_tile(tile, i, preread_ifgs):
        mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        # locally save the mst_mat
        mst_file_process_n = os.path.join(
            TMPDIR, 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t, p_ind in zip(process_tiles, process_indices):
        save_mst_tile(t, p_ind, preread_ifgs)


def get_process_tiles(dest_tifs, parallel, params):
    ifg = shared.pre_prepare_ifgs([dest_tifs[0]], params)[0]
    # TODO: changing the nrows and ncols here will break tests, fix this
    tiles = shared.create_tiles(ifg.shape,
                                n_ifgs=len(dest_tifs), nrows=3, ncols=4)
    no_tiles = len(tiles)
    process_indices = parallel.calc_indices(no_tiles)
    process_tiles = [itemgetter(p)(tiles) for p in process_indices]
    return ifg.shape, process_tiles, process_indices, tiles


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


if __name__ == '__main__':
    # read in the config file, and params for the simulation
    params = run_pyrate.get_ifg_paths()[2]
    main(params)
