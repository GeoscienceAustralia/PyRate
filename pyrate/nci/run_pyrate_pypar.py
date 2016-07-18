import cPickle as cp
import datetime
import glob
import os
from os.path import join
import sys
from collections import namedtuple
from operator import itemgetter
import numpy as np
from osgeo import gdal

from pyrate import config as cf
from pyrate import mst
from pyrate import refpixel
from pyrate import remove_aps_delay as aps
from pyrate import shared
from pyrate.nci import common_nci
from pyrate.nci.common_nci import save_latest_phase
from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate.scripts.run_pyrate import write_msg

gdal.SetCacheMax(64)
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
    num_processors = parallel.size

    # calculate process information
    ifg_shape, process_tiles, process_indices, tiles = \
        common_nci.get_process_tiles(dest_tifs, parallel, params)

    output_dir = params[cf.OUT_DIR]

    if rank == MASTER_PROCESS:
        print "Master process found {} worker processors".format(num_processors)

    preread_ifgs = convert_phase_to_numpy(dest_tifs, output_dir, tiles, parallel)

    parallel.barrier()
    mpi_log_filename = os.path.join(output_dir, "mpi_run_pyrate.log")

    if rank == MASTER_PROCESS:
        output_log_file = open(mpi_log_filename, "w")
        configfile = open(config_file)
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

    print 'Processor {} has {} tiles'.format(rank, len(process_tiles))
    # Calc mst using MPI
    if rank == MASTER_PROCESS:
        mpi_mst_calc(dest_tifs, process_tiles, process_indices,
                     preread_ifgs, params, parallel)
    else:
        mpi_mst_calc(dest_tifs, process_tiles, process_indices,
                     preread_ifgs, params, parallel)

    parallel.barrier()

    # Calc ref_pixel using MPI
    ref_pixel_file = os.path.join(params[cf.OUT_DIR], 'ref_pixel.npy')
    if rank == MASTER_PROCESS:
        refpx, refpy = ref_pixel_calc_mpi(rank, dest_tifs,
                                          num_processors, parallel, params)
        np.save(file=ref_pixel_file, arr=[refpx, refpy])
    else:
        ref_pixel_calc_mpi(rank, dest_tifs,
                           num_processors, parallel, params)

    parallel.barrier()

    # remove APS delay here
    if params[cf.APS_CORRECTION]:
        ifgs = shared.pre_prepare_ifgs(dest_tifs, params)
        if run_pyrate.aps_delay_required(ifgs, params):
            no_ifgs = len(ifgs)
            process_indices_aps = parallel.calc_indices(no_ifgs)
            ifgs = aps.remove_aps_delay(ifgs, params, process_indices_aps)

        for i in ifgs:
            i.close()

    parallel.finalize()


def convert_phase_to_numpy(dest_tifs, output_dir, tiles, parallel):
    """
    1. Convert ifg phase data into numpy files.
    2. Saves the preread_ifgs dict with information about the ifgs that are
    later used for fast loading of Ifg files in IfgPart class
    :param dest_tifs:
    :param output_dir:
    :param tiles:
    :param parallel:
    :return:
    """
    rank = parallel.rank
    num_processors = parallel.size

    preread_ifgs_dict = {}
    process_preread_ifgs = os.path.join(output_dir,
                                        'preread_ifgs_{}.pk'.format(rank))
    no_tifs = len(dest_tifs)
    process_indices = parallel.calc_indices(no_tifs)
    process_tifs = [itemgetter(p)(dest_tifs) for p in process_indices]
    for d in process_tifs:
        ifg = save_latest_phase(d, output_dir, tiles)
        nan_fraction = ifg.nan_fraction
        master = ifg.master
        slave = ifg.slave
        time_span = ifg.time_span

        preread_ifgs_dict[d] = {PrereadIfg(path=d,
                                           nan_fraction=nan_fraction,
                                           master=master,
                                           slave=slave,
                                           time_span=time_span)}
    cp.dump(preread_ifgs_dict, open(process_preread_ifgs, 'w'))
    parallel.barrier()
    preread_ifgs = os.path.join(output_dir, 'preread_ifgs.pk')
    if rank == MASTER_PROCESS:
        for r in xrange(1, num_processors):
            temp_file = join(output_dir, 'preread_ifgs_{}.pk'.format(r))
            temp_dict = cp.load(open(temp_file, 'r'))
            preread_ifgs_dict.update(temp_dict)
            os.remove(temp_file)
        cp.dump(preread_ifgs_dict, open(preread_ifgs, 'w'))
    print 'finish converting phase_data to numpy in process {}'.format(rank)
    return preread_ifgs


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


def mpi_mst_calc(dest_tifs, process_tiles, process_indices, preread_ifgs,
                 params, parallel):
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
    output_dir = params[cf.OUT_DIR]

    def save_mst_tile(tile, i, preread_ifgs):
        mst_tile = mst.mst_multiprocessing(tile, dest_tifs, preread_ifgs)
        # locally save the mst_mat
        mst_file_process_n = os.path.join(
            output_dir, 'mst_mat_{}.npy'.format(i))
        np.save(file=mst_file_process_n, arr=mst_tile)

    for t, p_ind in zip(process_tiles, process_indices):
        save_mst_tile(t, p_ind, preread_ifgs)
    print 'finished mst calculation for process {}'.format(parallel.rank)


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


if __name__ == '__main__':
    # read in the config file, and params for the simulation
    params = run_pyrate.get_ifg_paths()[2]
    main(params)
