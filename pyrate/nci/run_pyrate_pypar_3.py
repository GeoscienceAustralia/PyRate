import os
import sys
import numpy as np
from collections import namedtuple
from operator import itemgetter
from osgeo import gdal

from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import linrate
from pyrate import shared
from pyrate import timeseries
from pyrate.nci.parallel import Parallel
from pyrate.scripts import run_pyrate
from pyrate.scripts.run_pyrate import write_msg
from pyrate.nci import common_nci
gdal.SetCacheMax(64)

__author__ = 'sudipta'

# Constants
MASTER_PROCESS = 0
data_path = 'DATAPATH'

PrereadIfg = namedtuple('PrereadIfg', 'path nan_fraction master slave time_span')


def main(params):

    # setup paths
    xlks, ylks, crop = run_pyrate.transform_params(params)
    base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])
    dest_tifs = run_pyrate.get_dest_paths(base_unw_paths, crop, params, xlks)

    # Setting up parallelisation
    parallel = Parallel(True)
    rank = parallel.rank

    # calculate process information
    ifg_shape, process_tiles, process_indices, tiles = \
        common_nci.get_process_tiles(dest_tifs, parallel, params)

    preread_ifgs = os.path.join(params[cf.OUT_DIR], 'preread_ifgs.pk')

    maxvar_file = os.path.join(params[cf.OUT_DIR], 'maxvar.npy')
    vcmt_file = os.path.join(params[cf.OUT_DIR], 'vcmt.npy')

    maxvar, vcmt = np.load(maxvar_file), np.load(vcmt_file)
    output_dir = params[cf.OUT_DIR]

    # save latest phase data for use in linrate and mpi
    no_tifs = len(dest_tifs)
    process_indices = parallel.calc_indices(no_tifs)
    process_tifs = [itemgetter(p)(dest_tifs) for p in process_indices]
    for d in process_tifs:
        common_nci.save_latest_phase(d, output_dir, tiles)

    parallel.barrier()
    print 'linrate computation started'
    # linrate mpi computation
    linrate_mpi(dest_tifs, params, vcmt, process_tiles, preread_ifgs)

    # time series mpi computation
    if params[cf.TIME_SERIES_CAL]:
        time_series_mpi(dest_tifs, params, vcmt, process_tiles, preread_ifgs)
        parallel.barrier()
        write_time_series_geotiff_mpi(dest_tifs, params, tiles, parallel, rank)

    parallel.barrier()
    # linrate aggregation
    # TODO: identify processes with less load, make them do the remaining tasks
    if rank == 1:
        save_linrate_mpi(dest_tifs, params, tiles, out_type='linrate')
    elif rank == 2:
        save_linrate_mpi(dest_tifs, params, tiles, out_type='linerror')
    elif rank == 3:
        save_linrate_mpi(dest_tifs, params, tiles, out_type='linsamples')
    else:
        pass

    # TODO: Clean up

    parallel.finalize()


def save_linrate_mpi(dest_tifs, params, tiles, out_type):
    # TODO: fix linrate test
    # TODO: write tests for this function
    print 'saving linrate output type', out_type
    ifgs = shared.prepare_ifgs_without_phase(dest_tifs, params)
    epochlist, gt, md, wkt = run_pyrate.setup_metadata(ifgs, params)
    dest = os.path.join(params[cf.OUT_DIR], out_type + ".tif")
    md[ifc.MASTER_DATE] = epochlist.dates
    md[ifc.PRTYPE] = out_type
    output_dir = params[cf.OUT_DIR]
    rate = np.zeros(shape=ifgs[0].shape, dtype=np.float32)
    for t in tiles:
        rate_file = os.path.join(output_dir, out_type + '_{}.npy'.format(t.index))
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y,
            t.top_left_x:t.bottom_right_x] = rate_tile
    shared.write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    npy_rate_file = os.path.join(params[cf.OUT_DIR], out_type + '.npy')
    np.save(file=npy_rate_file, arr=rate)


def write_time_series_geotiff_mpi(dest_tifs, params, tiles, parallel, MPI_id):
    # TODO: This is currently pretty slow - investigate
    ifgs = shared.prepare_ifgs_without_phase(dest_tifs, params)
    epochlist, gt, md, wkt = run_pyrate.setup_metadata(ifgs, params)
    output_dir = params[cf.OUT_DIR]
    # load the first tsincr file to determine the number of time series tifs
    tsincr_file = os.path.join(output_dir, 'tsincr_0.npy')
    tsincr = np.load(file=tsincr_file)

    no_ts_tifs = tsincr.shape[2]
    # we create 2 x no_ts_tigs as we are splitting tsincr and tscuml
    # to all processes. This saves a lot of time.
    process_tifs = parallel.calc_indices(no_ts_tifs * 2)

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop is helps us over come the memory limit
    print 'process {} will write {} ts (incr/cuml) tifs of total {}'.format(
        MPI_id, len(process_tifs), no_ts_tifs * 2)
    for i in process_tifs:
        tscum_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        if i < no_ts_tifs:
            for n, t in enumerate(tiles):
                tscum_file = os.path.join(output_dir, 'tscuml_{}.npy'.format(n))
                tscum = np.load(file=tscum_file)

                md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
                md['PR_SEQ_POS'] = i  # sequence position
                tscum_g[t.top_left_y:t.bottom_right_y,
                    t.top_left_x:t.bottom_right_x] = tscum[:, :, i]
                dest = os.path.join(params[cf.OUT_DIR],
                    'tscuml' + "_" + str(epochlist.dates[i + 1]) + ".tif")
                md[ifc.PRTYPE] = 'tscuml'
                shared.write_output_geotiff(md, gt, wkt, tscum_g, dest, np.nan)
        else:
            tsincr_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
            i %= no_ts_tifs
            for n, t in enumerate(tiles):
                tsincr_file = os.path.join(output_dir, 'tsincr_{}.npy'.format(n))
                tsincr = np.load(file=tsincr_file)

                md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
                md['PR_SEQ_POS'] = i  # sequence position
                tsincr_g[t.top_left_y:t.bottom_right_y,
                t.top_left_x:t.bottom_right_x] = tsincr[:, :, i]
                dest = os.path.join(params[cf.OUT_DIR],
                                    'tsincr' + "_" + str(
                                        epochlist.dates[i + 1]) + ".tif")
                md[ifc.PRTYPE] = 'tsincr'
                shared.write_output_geotiff(md, gt, wkt, tsincr_g, dest,
                                            np.nan)
    print 'process {} finished writing {} ts (incr/cuml) tifs ' \
          'of total {}'.format(MPI_id, len(process_tifs), no_ts_tifs * 2)


def linrate_mpi(ifg_paths, params, vcmt, process_tiles, preread_ifgs):
    write_msg('Calculating linear rate')
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        i = t.index
        print 'calculating lin rate of tile {}'.format(i)
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


def time_series_mpi(ifg_paths, params, vcmt, process_tiles,
                    preread_ifgs):
    write_msg('Calculating time series')  # this should be logged
    output_dir = params[cf.OUT_DIR]
    for t in process_tiles:
        i = t.index
        print 'calculating time series for tile', t.index
        ifg_parts = [shared.IfgPart(p, t, preread_ifgs) for p in ifg_paths]
        mst_file_process_n = os.path.join(output_dir, 'mst_mat_{}.npy'.format(i))
        mst_tile = np.load(mst_file_process_n)
        res = timeseries.time_series(ifg_parts, params, vcmt, mst_tile)
        tsincr, tscum, tsvel = res
        tsincr_file = os.path.join(output_dir, 'tsincr_{}.npy'.format(i))
        tscum_file = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))
        np.save(file=tsincr_file, arr=tsincr)
        np.save(file=tscum_file, arr=tscum)


if __name__ == '__main__':
    # read in the config file, and params for the simulation
    paras = run_pyrate.get_ifg_paths()[2]
    main(paras)
