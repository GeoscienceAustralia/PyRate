import pyrate.shared

__author__ = 'sudipta'
import sys
import os
import datetime
from operator import itemgetter
import glob
import numpy as np

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
                + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
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
    if MPI_myID == MASTER_PROCESS:
        mst_grid, ifgs = mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs, mpi_log_filename,
                 num_processors, parallel, params)
        # write mst output to a file
        np.save(file=mst_mat_binary_file, arr=mst_grid)
    else:
        _, ifgs = mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs, mpi_log_filename,
                 num_processors, parallel, params)

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

    # every proces writing ifgs - BAD, TODO: Remove this step
    # required as all processes need orbital corrected ifgs
    run_pyrate.remove_orbital_error(ifgs, params)
    # orb_fit_calc_mpi(MPI_myID, ifgs, num_processors, parallel, params)

    # estimate reference phase

    _, ifgs = rpe.estimate_ref_phase(ifgs, params, refpx, refpy)

    # all processes need access to maxvar
    maxvar = maxvar_mpi(MPI_myID, ifgs, num_processors, parallel, params)
    # get vcmt is very fast and not mpi'ed at the moment
    # TODO: revisit this if performance is low in NCI
    vcmt = vcm_module.get_vcmt(ifgs, maxvar)  # all processes

    parallel.finalize()


def maxvar_mpi(MPI_myID, ifgs, num_processors, parallel, params):
    no_ifgs = len(ifgs)
    process_indices = parallel.calc_indices(no_ifgs)
    process_ifgs = [itemgetter(p)(ifgs) for p in process_indices]
    process_maxvar = [vcm_module.cvd(i)[0] for i in process_ifgs]
    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    if MPI_myID == MASTER_PROCESS:
        all_indices = parallel.calc_all_indices(no_ifgs)
        maxvar = np.empty(len(ifgs))
        maxvar[process_indices] = process_maxvar

        for i in range(1, num_processors):
            process_maxvar = parallel.receive(source=i, tag=-1,
                                              return_status=False)
            maxvar[all_indices[i]] = process_maxvar
        np.save(file=maxvar_file, arr=maxvar)
        original_maxvar = [vcm_module.cvd(i)[0] for i in ifgs]
        np.testing.assert_array_equal(original_maxvar, maxvar)
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
        mean_sds_final[all_indices[0]] = mean_sds
        for i in range(1, num_processors):
            process_mean_sds = parallel.receive(source=i, tag=-1,
                                                return_status=False)
            mean_sds_final[all_indices[i]] = process_mean_sds

        refx, refy = refpixel.filter_means(mean_sds_final, grid)
        print 'finished calculating ref pixel'
        return refx, refy
    else:
        parallel.send(mean_sds, destination=MASTER_PROCESS, tag=MPI_myID)
        print 'sent ref pixel requirements to master'


def mpi_mst_calc(MPI_myID, cropped_and_sampled_tifs, mpi_log_filename,
                 num_processors, parallel, params):
    ifgs = run_pyrate.prepare_ifgs_for_networkx_mst(cropped_and_sampled_tifs,
                                                    params)
    top_lefts, bottom_rights, no_tiles = pyrate.shared.setup_tiles(
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
        return result, ifgs
    else:
        # send the result arrays
        parallel.send(result_process, destination=MASTER_PROCESS, tag=MPI_myID)
        print "sent result from process", MPI_myID
        return result_process, ifgs


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


if __name__ == '__main__':
    main()
