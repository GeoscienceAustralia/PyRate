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
    parallel.finalize()


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
