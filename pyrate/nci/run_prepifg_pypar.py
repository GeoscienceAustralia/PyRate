__author__ = 'sudipta'
import sys
import os
import datetime
from operator import itemgetter

from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate import config as cf
from pyrate.scripts import run_prepifg
from pyrate import gamma
from pyrate.shared import Ifg
from pyrate import prepifg

# Constants
MASTER_PROCESS = 0

###============================================================================
#    Usage:  Reads in the gamma config file and produces the tifs from gamma
#    processed unwrapped interferrograms (ifgs).
#    mpirun -np 4 python pyrate/nci/run_prepifg_pypar.py gamma_config_file.conf
#===============================================================================


def main():

    # Setting up parallelisation
    parallel = Parallel(True)
    MPI_myID = parallel.rank
    num_processors = parallel.size
    ### Master Process ###
    if MPI_myID == MASTER_PROCESS:
        print "Master process found {} worker processors".format(num_processors)

    # Read config file, dest_paths are final mlooked/sampled and cropped tifs
    base_ifg_paths, dest_paths, params = run_pyrate.get_ifg_paths()

    # logfile
    output_dir = params[cf.OUT_DIR]
    mpi_log_filename = os.path.join(output_dir, "mpi.log")

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

        output_log_file.write("\nInput unwrapped interferrograms:\n")
        for b in base_ifg_paths:
            output_log_file.write(b + "\n")

        output_log_file.write("\nOutput gamma processed interferrograms:\n")
        for b in dest_paths:
            output_log_file.write(b + "\n")

        output_log_file.close()

    num_files = len(dest_paths)
    parallel.calc_lo_hi(num_files)
    process_subset_indices = parallel.calc_indices(num_files)

    process_base_paths = [itemgetter(p)(base_ifg_paths)
                          for p in process_subset_indices]

    print 'Processor {mpi_id} has {processes} ' \
          'sites out of {num_files}'.format(mpi_id=MPI_myID,
                                            processes=len(process_base_paths),
                                            num_files=num_files)

    msg = "running gamma prepifg"
    print msg

    # location of geo_tif's
    dest_base_ifgs = [os.path.join(
        params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
                      for q in base_ifg_paths]

    # calculate geo_tifs using mpi
    [run_prepifg.gamma_multiprocessing(b, params) for b in process_base_paths]

    # need to come back to the main thread as all ifgs are needed for exts calc
    parallel.barrier()
    ifgs = [Ifg(p) for p in dest_base_ifgs]
    xlooks, ylooks, crop = run_pyrate.transform_params(params)
    exts = prepifg.getAnalysisExtent(crop, ifgs, xlooks, ylooks, userExts=None)
    thresh = params[cf.NO_DATA_AVERAGING_THRESHOLD]
    del ifgs  # to save memory
    # go mpi again for prep_ifg
    for p in process_subset_indices:
        data_path = itemgetter(p)(dest_base_ifgs)
        prepifg.prepare_ifg(data_path, xlooks, ylooks, exts, thresh, crop)
    parallel.finalize()


if __name__ == "__main__":
    main()
