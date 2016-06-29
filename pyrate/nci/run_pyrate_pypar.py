import datetime
import glob
import os
import sys
from operator import itemgetter

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
from pyrate.algorithm import get_epochs
from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate.scripts.run_pyrate import write_msg
from pyrate.shared import get_tmpdir

TMPDIR = get_tmpdir()

__author__ = 'sudipta'

# Constants
MASTER_PROCESS = 0


def main(params=None):
    # Setting up parallelisation
    parallel = Parallel(True)
    MPI_myID = parallel.rank
    num_processors = parallel.size
    ### Master Process ###
    if MPI_myID == MASTER_PROCESS:
        print "Master process found {} worker processors".format(num_processors)

    # Read config file, dest_tifs are input files to run_pyrate
    if params:
        xlks, ylks, crop = run_pyrate.transform_params(params)
        base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])
        dest_tifs = \
            run_pyrate.get_dest_paths(base_unw_paths, crop, params, xlks)
    else:
        _, dest_tifs, params = run_pyrate.get_ifg_paths()

    output_dir = params[cf.OUT_DIR]
    mpi_log_filename = os.path.join(output_dir, "mpi_run_pyrate.log")

    ### Master Process ###
    if MPI_myID == MASTER_PROCESS:
        output_log_file = open(mpi_log_filename, "w")
        config_filepath = sys.argv[1]
        configfile = open(config_filepath)
        output_log_file.write("Starting Simulation at: "
                              + datetime.datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"))
        output_log_file.write("Master process found " +
                              str(num_processors) +
                              " worker processors.\n")
        output_log_file.write("\n")
        output_log_file.write("\nConfig Settings: start\n")
        lines = configfile.read()
        for line in lines:
            output_log_file.write(line)
        output_log_file.write("\nConfig Settings: end\n")

        output_log_file.write("\n Input files for run_pyrate are:\n")
        for b in dest_tifs:
            output_log_file.write(b + "\n")

        output_log_file.close()

    parallel.barrier()

    # calculate process information
    ifg_shape, process_tiles, process_indices = \
        get_process_tiles(dest_tifs, parallel, params)

    print 'Processor {} has {} tiles'.format(MPI_myID, len(process_tiles))
    # Calc mst using MPI
    if MPI_myID == MASTER_PROCESS:
        mpi_mst_calc(dest_tifs, process_tiles, process_indices)
    else:
        mpi_mst_calc(dest_tifs, process_tiles, process_indices)

    if MPI_myID == MASTER_PROCESS:
        output_log_file = open(mpi_log_filename, "a")
        output_log_file.write("\n\n Mst caclulation finished\n")
        output_log_file.close()

    parallel.barrier()

    # Calc ref_pixel using MPI
    ref_pixel_file = os.path.join(params[cf.OUT_DIR], 'ref_pixel.npy')
    if MPI_myID == MASTER_PROCESS:
        refpx, refpy = ref_pixel_calc_mpi(MPI_myID, dest_tifs,
                                          num_processors, parallel, params)
        np.save(file=ref_pixel_file, arr=[refpx, refpy])
    else:
        ref_pixel_calc_mpi(MPI_myID, dest_tifs,
                           num_processors, parallel, params)

    parallel.barrier()
    # refpixel read in each process
    refpx, refpy = np.load(ref_pixel_file)
    print 'Found reference pixel', refpx, refpy
    parallel.barrier()

    # remove APS delay here
    if params[cf.APS_CORRECTION]:
        ifgs = shared.pre_prepare_ifgs(dest_tifs, params)
        if run_pyrate.aps_delay_required(ifgs, params):
            no_ifgs = len(ifgs)
            process_indices = parallel.calc_indices(no_ifgs)
            # process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
            ifgs = aps.remove_aps_delay(ifgs, params, process_indices)

        for i in ifgs:
            i.close()

    parallel.barrier()
    # required as all processes need orbital corrected ifgs
    orb_fit_calc_mpi(MPI_myID, dest_tifs,
                     num_processors, parallel, params)
    parallel.barrier()

    # estimate and remove reference phase
    ref_phase_estimation_mpi(MPI_myID, dest_tifs, parallel,
                             params, refpx, refpy)
    parallel.barrier()

    # all processes need access to maxvar, and vcmt
    maxvar, vcmt = maxvar_vcm_mpi(MPI_myID, dest_tifs, parallel, params)

    vcmt_file = os.path.join(params[cf.OUT_DIR], 'vcmt.npy')
    if MPI_myID == MASTER_PROCESS:
        np.save(file=vcmt_file, arr=vcmt)

    parallel.barrier()

    # linrate mpi computation
    linrate_mpi(MPI_myID, dest_tifs, parallel, params, vcmt,
                process_tiles, process_indices, ifg_shape)

    # time series mpi computation
    if params[cf.TIME_SERIES_CAL]:
        time_series_mpi(MPI_myID, dest_tifs, parallel, params, vcmt,
                        process_tiles, process_indices, ifg_shape)

    parallel.finalize()


def linrate_mpi(MPI_myID, ifg_paths, parallel, params, vcmt,
                process_tiles, process_indices, ifg_shape):
    write_msg('Calculating linear rate')
    num_processors = parallel.size
    # initialize empty arrays
    rows, cols = ifg_shape
    rate = np.zeros([rows, cols], dtype=np.float32)
    # may not need error and samples? save memory?
    error = np.zeros([rows, cols], dtype=np.float32)
    samples = np.zeros([rows, cols], dtype=np.float32)

    for i, t in zip(process_indices, process_tiles):
        r_start, c_start = t.top_left
        r_end, c_end = t.bottom_right
        ifg_parts = [shared.IfgPart(p, t) for p in ifg_paths]
        mst_file_process_n = os.path.join(
            TMPDIR, 'mst_mat_{}.npy'.format(i))
        mst_grid_ifg_parts = np.load(mst_file_process_n)
        # mst_grid_ifg_parts = mst_grid[:, r_start: r_end, c_start: c_end]
        res = linrate.linear_rate(ifg_parts, params, vcmt, mst_grid_ifg_parts)

        for r in res:
            if r is None:
                raise ValueError('TODO: bad value')
        rate[r_start:r_end, c_start:c_end] = res[0]
        error[r_start:r_end, c_start:c_end] = res[1]
        samples[r_start:r_end, c_start:c_end] = res[2]

    process_res = (rate, error, samples)

    if MPI_myID == MASTER_PROCESS:
        for i in range(1, num_processors):
            process_res = parallel.receive(source=i, tag=-1,
                                           return_status=False)
            rate += process_res[0]
            error += process_res[1]
            samples += process_res[2]

        # declare file names
        rate_file = os.path.join(params[cf.OUT_DIR], 'rate.npy')
        error_file = os.path.join(params[cf.OUT_DIR], 'error.npy')
        samples_file = os.path.join(params[cf.OUT_DIR], 'samples.npy')

        np.save(file=rate_file, arr=rate)
        np.save(file=error_file, arr=error)
        np.save(file=samples_file, arr=samples)
    else:
        parallel.send(process_res, destination=MASTER_PROCESS, tag=MPI_myID)


def time_series_mpi(MPI_myID, ifg_paths, parallel, params, vcmt,
                    process_tiles, process_indices, ifg_shape):
    write_msg('Calculating time series')

    num_processors = parallel.size

    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    epochlist = get_epochs(ifgs)
    nepoch = len(epochlist.dates)  # epoch number
    nvelpar = nepoch - 1  # velocity parameters number

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    tsvel = np.zeros(shape=ifg_shape + (nvelpar, ), dtype=np.float32)
    tsincr = np.zeros(shape=ifg_shape + (nvelpar, ), dtype=np.float32)
    tscum = np.zeros(shape=ifg_shape + (nvelpar, ), dtype=np.float32)

    for i, t in zip(process_indices, process_tiles):
        r_start, c_start = t.top_left
        r_end, c_end = t.bottom_right
        ifg_parts = [shared.IfgPart(p, t) for p in ifg_paths]
        mst_file_process_n = os.path.join(
            TMPDIR, 'mst_mat_{}.npy'.format(i))
        mst_grid_ifg_parts = np.load(mst_file_process_n)
        # mst_grid_ifg_parts = mst_grid[:, r_start: r_end, c_start: c_end]
        res = timeseries.time_series(ifg_parts, params,
                                     vcmt, mst_grid_ifg_parts)
        tsincr[r_start:r_end, c_start:c_end, :] = res[0]
        tscum[r_start:r_end, c_start:c_end, :] = res[1]
        tsvel[r_start:r_end, c_start:c_end, :] = res[2]

    process_res = (tsincr, tscum, tsvel)

    if MPI_myID == MASTER_PROCESS:

        # putting it together
        for i in range(1, num_processors):
            process_res = parallel.receive(source=i, tag=-1, return_status=False)
            tsincr += process_res[0]
            tscum += process_res[1]
            tsvel += process_res[2]
        tsvel_file = os.path.join(params[cf.OUT_DIR], 'tsvel.npy')
        tsincr_file = os.path.join(params[cf.OUT_DIR], 'tsincr.npy')
        tscum_file = os.path.join(params[cf.OUT_DIR], 'tscum.npy')
        np.save(file=tscum_file, arr=tscum)
        np.save(file=tsincr_file, arr=tsincr)
        np.save(file=tsvel_file, arr=tsvel)
    else:
        parallel.send(process_res, destination=MASTER_PROCESS,
                      tag=MPI_myID)


def ref_phase_estimation_mpi(MPI_myID, ifg_paths, parallel, params,
                             refpx, refpy):
    write_msg('Finding and removing reference phase')
    num_processors = parallel.size
    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_ref_phs = np.zeros(len(process_ifgs))

    if params[cf.REF_EST_METHOD] == 1:
        ifg_phase_data_sum = np.zeros(ifgs[0].shape, dtype=np.float64)
        # TODO: revisit as this will likely hit memory limit in NCI
        for ifg in ifgs:
            ifg_phase_data_sum += ifg.phase_data

        comp = np.isnan(ifg_phase_data_sum)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
        for n, ifg in enumerate(process_ifgs):
            # shared.nan_and_mm_convert(ifg, params)
            ref_phs = rpe.est_ref_phase_method1_multi(ifg.phase_data, comp)
            ifg.phase_data -= ref_phs
            ifg.meta_data[ifc.REF_PHASE] = ifc.REF_PHASE_REMOVED
            ifg.write_modified_phase()
            process_ref_phs[n] = ref_phs
            ifg.close()

    elif params[cf.REF_EST_METHOD] == 2:
        half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
        chipsize = 2 * half_chip_size + 1
        thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
        for n, ifg in enumerate(process_ifgs):
            # shared.nan_and_mm_convert(ifg, params)
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
        ref_phs = np.zeros(len(ifgs))
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


def maxvar_vcm_mpi(MPI_myID, ifg_paths, parallel, params):
    num_processors = parallel.size
    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_maxvar = [vcm_module.cvd(i, params)[0] for i in process_ifgs]
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        maxvar = np.empty(len(ifgs))
        maxvar[all_indices[MASTER_PROCESS]] = process_maxvar

        for i in range(1, num_processors):
            process_maxvar = parallel.receive(source=i, tag=-1,
                                              return_status=False)
            maxvar[all_indices[i]] = process_maxvar
        np.save(file=maxvar_file, arr=maxvar)
    else:
        parallel.send(process_maxvar, destination=MASTER_PROCESS, tag=MPI_myID)
    parallel.barrier()
    maxvar = np.load(maxvar_file)
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)

    return maxvar, vcmt


def orb_fit_calc_mpi(MPI_myID, ifg_paths, num_processors, parallel, params):
    print 'calculating orbfit using MPI id:', MPI_myID
    if params[cf.ORBITAL_FIT_METHOD] != 1:
        raise cf.ConfigException('For now orbfit method must be 1')

    ifgs = shared.prepare_ifgs_without_phase(ifg_paths, params)
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]

    mlooked = None
    # difficult to enable MPI on orbfit method 2
    # so just do orbfit method 1 for now
    # TODO: MPI orbfit method 2
    orbital.orbital_correction(process_ifgs,
                               params,
                               mlooked=mlooked)
    # set orbfit tags after orbital error correction
    for i in process_ifgs:
        i.write_modified_phase()
        # implement mpi logging


def ref_pixel_calc_mpi(MPI_myID, ifg_paths, num_processors, parallel, params):
    half_patch_size, thresh, grid = refpixel.ref_pixel_setup(ifg_paths, params)
    no_steps = len(grid)
    process_indices = parallel.calc_indices(no_steps)
    process_grid = [itemgetter(p)(grid) for p in process_indices]
    print 'Processor {mpi_id} has {processes} ' \
          'tiles out of {num_files}'.format(mpi_id=MPI_myID,
                                            processes=len(process_indices),
                                            num_files=no_steps)
    mean_sds = refpixel.ref_pixel_mpi(process_grid,
                                      half_patch_size, ifg_paths, thresh, params)
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


def mpi_mst_calc(dest_tifs, process_tiles, process_indices):
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

    def save_mst_tile(tile, i):
        mst_tile = mst.mst_multiprocessing(tile, dest_tifs)
        # locally save the mst_mat
        mst_file_process_n = os.path.join(
            TMPDIR, 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t, p_ind in zip(process_tiles, process_indices):
        save_mst_tile(t, p_ind)


def get_process_tiles(dest_tifs, parallel, params):
    ifg = shared.pre_prepare_ifgs([dest_tifs[0]], params)[0]
    tiles = shared.create_tiles(ifg.shape,
                                n_ifgs=len(dest_tifs), nrows=3, ncols=4)
    no_tiles = len(tiles)
    process_indices = parallel.calc_indices(no_tiles)
    process_tiles = [itemgetter(p)(tiles) for p in process_indices]
    return ifg.shape, process_tiles, process_indices


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


if __name__ == '__main__':
    main()
