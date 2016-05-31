import sys
import os
import datetime
from operator import itemgetter
import glob
import numpy as np
from numpy import isnan, nan

from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate import config as cf
from pyrate.scripts import run_prepifg
from pyrate import gamma
from pyrate.shared import Ifg
from pyrate import prepifg
from pyrate import mst
from pyrate import refpixel
from pyrate import orbital
from pyrate import ifgconstants as ifc
from pyrate import vcm as vcm_module
from pyrate import reference_phase_estimation as rpe
from pyrate import timeseries
from pyrate import linrate
from pyrate import remove_aps_delay as aps
from pyrate import shared
import gdal
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

    # Read config file, cropped_and_sampled_tifs are input files to run_pyrate
    if params:
        xlks, ylks, crop = run_pyrate.transform_params(params)
        base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])
        cropped_and_sampled_tifs = \
            run_pyrate.get_dest_paths(base_unw_paths, crop, params, xlks)
    else:
        _, cropped_and_sampled_tifs, params = run_pyrate.get_ifg_paths()

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
        for b in cropped_and_sampled_tifs:
            output_log_file.write(b + "\n")

        output_log_file.close()

    parallel.barrier()

    # Calc mst using MPI
    mst_mat_binary_file = os.path.join(params[cf.OUT_DIR], 'mst_mat.npy')
    run_pyrate.write_msg('Calculating mst')
    if MPI_myID == MASTER_PROCESS:
        _, ifgs = mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs,
                               mpi_log_filename,
                               num_processors, parallel, params,
                               mst_mat_binary_file)
    else:
        _, ifgs = mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs,
                               mpi_log_filename,
                               num_processors, parallel, params,
                               mst_mat_binary_file)
    print 'finished mst computation'

    parallel.barrier()

    mst_grid = np.load(file=mst_mat_binary_file)

    # Calc ref_pixel using MPI
    ref_pixel_file = os.path.join(params[cf.OUT_DIR], 'ref_pixel.npy')

    if MPI_myID == MASTER_PROCESS:
        refpx, refpy = ref_pixel_calc_mpi(MPI_myID, ifgs,
                                          num_processors, parallel, params)
        np.save(file=ref_pixel_file, arr=[refpx, refpy])
    else:
        ref_pixel_calc_mpi(MPI_myID, ifgs, num_processors, parallel, params)

    parallel.barrier()
    # refpixel read in each process
    refpx, refpy = np.load(ref_pixel_file)

    parallel.barrier()

    # remove APS delay here
    if run_pyrate.aps_delay_required(ifgs, params):
        no_ifgs = len(ifgs)
        process_indices = parallel.calc_indices(no_ifgs)
        # process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
        ifgs = aps.remove_aps_delay(ifgs, params, process_indices)

    parallel.barrier()

    parallel.barrier()
    # every proces writing ifgs - BAD, TODO: Remove this step
    # required as all processes need orbital corrected ifgs
    run_pyrate.remove_orbital_error(ifgs, params)
    # orb_fit_calc_mpi(MPI_myID, ifgs, num_processors, parallel, params)

    # estimate and remove reference phase, ifgs are modifed in place
    ref_phase_estimation_mpi(MPI_myID, ifgs, num_processors, parallel, params,
                             refpx, refpy)

    # all processes need access to maxvar
    maxvar = maxvar_mpi(MPI_myID, ifgs, num_processors, parallel, params)
    # get vcmt is very fast and not mpi'ed at the moment
    # TODO: revisit this if performance is low in NCI
    # vcmt extremely fast already and no need to parallize
    # offers easy parallisation if required
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)  # all processes

    vcmt_file = os.path.join(params[cf.OUT_DIR], 'vcmt.npy')
    if MPI_myID == MASTER_PROCESS:
        np.save(file=vcmt_file, arr=vcmt)

    parallel.barrier()

    # # time series mpi computation
    time_series_mpi(MPI_myID, ifgs, mst_grid, num_processors, parallel, params,
                    vcmt)

    parallel.barrier()  # may not need this

    # linrate mpi computation
    linrate_mpi(MPI_myID, ifgs, mst_grid, num_processors, parallel, params,
                vcmt)

    parallel.finalize()


def linrate_mpi(MPI_myID, ifgs, mst_grid, num_processors, parallel, params,
                vcmt):
    run_pyrate.write_msg('Calculating linear rate')
    MAXSIG, NSIG, PTHRESH, ncols, error, mst, obs, _, processes, rate, \
    nrows, samples, span = linrate.linrate_setup(ifgs, mst_grid, params)
    rows = range(nrows)
    process_indices = parallel.calc_indices(nrows)
    process_rows = [itemgetter(p)(rows) for p in process_indices]
    process_res = np.empty(shape=(len(process_rows), ncols, 3))
    for n, row in enumerate(process_rows):
        process_res[n, :, :] = linrate.linear_rate_by_rows(
            row, ncols, mst, NSIG, obs, PTHRESH, span, vcmt)
    rate_file = os.path.join(params[cf.OUT_DIR], 'rate.npy')
    error_file = os.path.join(params[cf.OUT_DIR], 'error.npy')
    samples_file = os.path.join(params[cf.OUT_DIR], 'samples.npy')
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(nrows)
        # rate, error, samples have already been zero initialized in setup
        rate[all_indices[MASTER_PROCESS]] = process_res[:, :, 0]
        error[all_indices[MASTER_PROCESS]] = process_res[:, :, 1]
        samples[all_indices[MASTER_PROCESS]] = process_res[:, :, 2]
        for i in range(1, num_processors):
            process_res = parallel.receive(source=i, tag=-1,
                                           return_status=False)
            rate[all_indices[i]] = process_res[:, :, 0]
            error[all_indices[i]] = process_res[:, :, 1]
            samples[all_indices[i]] = process_res[:, :, 2]
        mask = ~isnan(error)
        mask[mask] &= error[mask] > MAXSIG
        rate[mask] = nan
        error[mask] = nan
        np.save(file=rate_file, arr=rate)
        np.save(file=error_file, arr=error)
        np.save(file=samples_file, arr=samples)
    else:
        parallel.send(process_res, destination=MASTER_PROCESS,
                      tag=MPI_myID)


def time_series_mpi(MPI_myID, ifgs, mst_grid, num_processors, parallel, params,
                    vcmt):
    run_pyrate.write_msg('Calculating time series')
    B0, INTERP, PTHRESH, SMFACTOR, SMORDER, TSMETHOD, ifg_data, mst, \
    ncols, nrows, nvelpar, _, processes, span, tsvel_matrix = \
        timeseries.time_series_setup(ifgs=ifgs, mst=mst_grid, params=params)
    # mpi by the rows
    # still problem is ifg_data which contains data for all ifgs
    # TODO: must stop sending all phase_data to each mpi process
    rows = range(nrows)
    process_indices = parallel.calc_indices(nrows)
    process_rows = [itemgetter(p)(rows) for p in process_indices]
    process_tsvel_matrix = []
    for row in process_rows:
        process_tsvel_matrix.append(timeseries.time_series_by_rows(
            row, B0, SMFACTOR, SMORDER, ifg_data, mst, ncols,
            nvelpar, PTHRESH, vcmt, TSMETHOD, INTERP))
    tsvel_file = os.path.join(params[cf.OUT_DIR], 'tsvel.npy')
    tsincr_file = os.path.join(params[cf.OUT_DIR], 'tsincr.npy')
    tscum_file = os.path.join(params[cf.OUT_DIR], 'tscum.npy')
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(nrows)
        tsvel_matrix = np.empty(shape=(nrows, ncols, nvelpar), dtype=np.float32)
        tsvel_matrix[all_indices[MASTER_PROCESS]] = process_tsvel_matrix

        # putting it together
        for i in range(1, num_processors):
            process_tsvel_matrix = parallel.receive(source=i, tag=-1,
                                                    return_status=False)
            tsvel_matrix[all_indices[i]] = process_tsvel_matrix
        tsvel_matrix = np.where(tsvel_matrix == 0, np.nan, tsvel_matrix)
        tsincr = tsvel_matrix * span
        tscum = np.cumsum(tsincr, 2)
        np.save(file=tsvel_file, arr=tsvel_matrix)
        np.save(file=tscum_file, arr=tscum)
        np.save(file=tsincr_file, arr=tsincr)
    else:
        parallel.send(process_tsvel_matrix, destination=MASTER_PROCESS,
                      tag=MPI_myID)


def ref_phase_estimation_mpi(MPI_myID, ifgs, num_processors, parallel, params,
                             refpx, refpy):
    run_pyrate.write_msg('Finding and removing reference phase')
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_phase_data = [i.phase_data for i in process_ifgs]
    process_ref_phs = np.zeros(len(process_ifgs))
    if params[cf.REF_EST_METHOD] == 1:
        ifg_phase_data_sum = np.zeros(ifgs[0].shape, dtype=np.float64)
        # TODO: revisit as this will likely hit memory limit in NCI
        for ifg in ifgs:
            ifg_phase_data_sum += ifg.phase_data

        comp = np.isnan(ifg_phase_data_sum)  # this is the same as in Matlab
        comp = np.ravel(comp, order='F')  # this is the same as in Matlab
        for n, ifg in enumerate(process_ifgs):
            process_ref_phs[n] = \
                rpe.est_ref_phase_method1_multi(process_phase_data[n], comp)

    elif params[cf.REF_EST_METHOD] == 2:
        half_chip_size = int(np.floor(params[cf.REF_CHIP_SIZE] / 2.0))
        chipsize = 2 * half_chip_size + 1
        thresh = chipsize * chipsize * params[cf.REF_MIN_FRAC]
        for n, ifg in enumerate(process_ifgs):
            process_ref_phs[n] = \
                rpe.est_ref_phase_method2_multi(process_phase_data[n],
                                                half_chip_size,
                                                refpx, refpy, thresh)

    else:
        raise
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
    parallel.barrier()
    ref_phs = np.load(ref_phs_file)
    for n, ifg in enumerate(ifgs):
        ifg.phase_data -= ref_phs[n]


def maxvar_mpi(MPI_myID, ifgs, num_processors, parallel, params):
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_maxvar = [vcm_module.cvd(i)[0] for i in process_ifgs]
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
    return maxvar


def orb_fit_calc_mpi(MPI_myID, ifgs, num_processors, parallel, params):
    print 'calculating orbfit in parallel'
    # exts calculation is very fast not mpi-ed
    exts = prepifg.getAnalysisExtent(
        cropOpt=prepifg.ALREADY_SAME_SIZE,
        rasters=ifgs,
        xlooks=params[cf.ORBITAL_FIT_LOOKS_X],
        ylooks=params[cf.ORBITAL_FIT_LOOKS_Y],
        userExts=None)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    data_paths = [i.data_path for i in ifgs]
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_data_paths = [itemgetter(p)(data_paths) for p in process_indices]
    process_mlooked_dataset = [prepifg.prepare_ifg(
        d,
        params[cf.ORBITAL_FIT_LOOKS_X],
        params[cf.ORBITAL_FIT_LOOKS_Y],
        exts,
        thresh,
        prepifg.ALREADY_SAME_SIZE,
        False) for d in process_data_paths]
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        mlooked_dataset = list(np.empty(shape=no_ifgs))
        # process_indices and all_indices[0] will be same for master process
        for i, p in enumerate(process_indices):
            mlooked_dataset[p] = process_mlooked_dataset[i]

        # this loop does not work due to swig import pickling issue
        for i in range(1, num_processors):
            process_mlooked_dataset = parallel.receive(source=i, tag=-1,
                                                       return_status=False)
            for i, p in enumerate(all_indices[i]):
                mlooked_dataset[p] = process_mlooked_dataset[i]

        mlooked = [Ifg(m) for m in mlooked_dataset]

        for m, i in zip(mlooked, ifgs):
            m.initialize()
            m.nodata_value = params[cf.NO_DATA_VALUE]

        # difficult to enable MPI on orbfit method 2
        # so just do orbfit method 1 for now
        # TODO: MPI orbfit method 2
        orbital.orbital_correction(ifgs,
                                   params,
                                   mlooked=mlooked)
        # write data to disc after orbital error correction
        for i in ifgs:
            i.dataset.SetMetadataItem(ifc.PYRATE_ORBITAL_ERROR,
                                      run_pyrate.ORB_REMOVED)
            i.write_modified_phase()
    else:
        parallel.send(process_mlooked_dataset, destination=MASTER_PROCESS,
                      tag=MPI_myID)
        print 'sent mlooked phase data'


def ref_pixel_calc_mpi(MPI_myID, ifgs, num_processors, parallel, params):
    half_patch_size, _, thresh, grid = refpixel.ref_pixel_setup(ifgs, params)
    no_steps = len(grid)
    process_indices = parallel.calc_indices(no_steps)
    process_grid = [itemgetter(p)(grid) for p in process_indices]
    print 'Processor {mpi_id} has {processes} ' \
          'tiles out of {num_files}'.format(mpi_id=MPI_myID,
                                            processes=len(process_indices),
                                            num_files=no_steps)
    mean_sds = refpixel.ref_pixel_mpi(process_grid,
                                      half_patch_size, ifgs, thresh)
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


def mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs, mpi_log_filename,
                 num_processors, parallel, params,
                 mst_file):
    run_pyrate.write_msg('Calculating minimum spanning tree matrix '
                         'using NetworkX method')
    ifgs = run_pyrate.pre_prepare_ifgs(cropped_and_sampled_tifs,
                                       params)
    top_lefts, bottom_rights, no_tiles = shared.setup_tiles(
        ifgs[0].shape, processes=num_processors)
    # parallel.calc_indices(no_tiles)
    process_indices = parallel.calc_indices(no_tiles)
    process_top_lefts = [itemgetter(p)(top_lefts)
                         for p in process_indices]
    process_bottom_rights = [itemgetter(p)(bottom_rights)
                             for p in process_indices]
    print 'Processor {mpi_id} has {processes} ' \
          'tiles out of {num_files}'.format(mpi_id=MPI_myID,
                                            processes=len(process_indices),
                                            num_files=no_tiles)
    result_process = mst.mst_multiprocessing_map(
        process_top_lefts, process_bottom_rights,
        ifgs, ifgs[0].shape, no_ifgs=len(ifgs)
    )
    parallel.barrier()

    if MPI_myID == MASTER_PROCESS:
        result = result_process
        # combine the mst from the other processes
        for i in range(1, num_processors):
            result_remote_processes = \
                parallel.receive(source=i, tag=-1, return_status=False)
            result += result_remote_processes

        output_log_file = open(mpi_log_filename, "a")
        output_log_file.write("\n\n Mst caclulation finished\n")
        output_log_file.close()
        # write mst output to a file
        np.save(file=mst_file, arr=result)
        return result, ifgs
    else:
        # send the result arrays
        parallel.send(result_process, destination=MASTER_PROCESS, tag=MPI_myID)
        print "sent mst result from process", MPI_myID
        return result_process, ifgs


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


if __name__ == '__main__':
    main()
