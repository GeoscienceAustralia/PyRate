import os
from os.path import join
import logging
import numpy as np
from osgeo import gdal
import pickle as cp

from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import shared
from pyrate.scripts import run_pyrate
from pyrate import mpiops
from pyrate.shared import PrereadIfg
gdal.SetCacheMax(64)

log = logging.getLogger(__name__)

# Constants
MASTER_PROCESS = 0


def main(config_file, rows, cols):
    # setup paths
    base_unw_paths, dest_paths, params= cf.get_ifg_paths(config_file)
    postprocess_linrate(rows, cols, params)
    if params[cf.TIME_SERIES_CAL]:
        postprocess_timeseries(rows, cols, params)


def postprocess_linrate(rows, cols, params):
    # setup paths
    xlks, ylks, crop = cf.transform_params(params)
    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])
    dest_tifs = cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    # load previously saved prepread_ifgs dict
    preread_ifgs_file = join(params[cf.OUT_DIR], 'preread_ifgs.pk')
    ifgs = cp.load(open(preread_ifgs_file, 'rb'))
    tiles = run_pyrate.get_tiles(dest_tifs[0], rows, cols)

    # linrate aggregation
    if mpiops.size >= 3:
        [save_linrate(ifgs, params, tiles, out_type=t)
         for i, t in enumerate(['linrate', 'linerror', 'linsamples'])
         if i == mpiops.rank]
    else:
        if mpiops.rank == MASTER_PROCESS:
            [save_linrate(ifgs, params, tiles, out_type=t)
             for t in ['linrate', 'linerror', 'linsamples']]


def save_linrate(ifgs_dict, params, tiles, out_type):
    log.info('Stating postprocessing {}'.format(out_type))
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, PrereadIfg)]
    dest = os.path.join(params[cf.OUT_DIR], out_type + ".tif")
    md[ifc.MASTER_DATE] = epochlist.dates
    md[ifc.PRTYPE] = out_type
    output_dir = params[cf.OUT_DIR]
    rate = np.zeros(shape=ifgs[0].shape, dtype=np.float32)
    for t in tiles:
        rate_file = os.path.join(output_dir, out_type +
                                 '_{}.npy'.format(t.index))
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y,
             t.top_left_x:t.bottom_right_x] = rate_tile
    shared.write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    npy_rate_file = os.path.join(params[cf.OUT_DIR], out_type + '.npy')
    np.save(file=npy_rate_file, arr=rate)
    log.info('Finished postprocessing {}'.format(out_type))


def postprocess_timeseries(rows, cols, params):
    xlks, ylks, crop = cf.transform_params(params)
    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])
    dest_tifs = cf.get_dest_paths(base_unw_paths, crop, params, xlks)
    output_dir = params[cf.OUT_DIR]

    # load previously saved prepread_ifgs dict
    preread_ifgs_file = join(output_dir, 'preread_ifgs.pk')
    ifgs = cp.load(open(preread_ifgs_file, 'rb'))

    # metadata and projections
    gt, md, wkt = ifgs['gt'], ifgs['md'], ifgs['wkt']
    epochlist = ifgs['epochlist']
    ifgs = [v for v in ifgs.values() if isinstance(v, PrereadIfg)]

    tiles = run_pyrate.get_tiles(dest_tifs[0], rows, cols)

    # load the first tsincr file to determine the number of time series tifs
    tsincr_file = os.path.join(output_dir, 'tsincr_0.npy')
    tsincr = np.load(file=tsincr_file)

    no_ts_tifs = tsincr.shape[2]
    # we create 2 x no_ts_tifs as we are splitting tsincr and tscuml
    # to all processes.
    process_tifs = mpiops.array_split(range(2 * no_ts_tifs))

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop helps us overcome the memory limit
    log.info('process {} will write {} ts (incr/cuml) tifs '
             'of total {}'.format(mpiops.rank, len(process_tifs),
                                  no_ts_tifs * 2))
    for i in process_tifs:
        tscum_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        if i < no_ts_tifs:
            for n, t in enumerate(tiles):
                tscum_file = os.path.join(output_dir,
                                          'tscuml_{}.npy'.format(n))
                tscum = np.load(file=tscum_file)

                md[ifc.MASTER_DATE] = epochlist.dates[i + 1]
                md['PR_SEQ_POS'] = i  # sequence position
                tscum_g[t.top_left_y:t.bottom_right_y,
                    t.top_left_x:t.bottom_right_x] = tscum[:, :, i]
                dest = os.path.join(params[cf.OUT_DIR],
                                    'tscuml' + "_" +
                                    str(epochlist.dates[i + 1]) + ".tif")
                md[ifc.PRTYPE] = 'tscuml'
                shared.write_output_geotiff(md, gt, wkt, tscum_g, dest, np.nan)
        else:
            tsincr_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
            i %= no_ts_tifs
            for n, t in enumerate(tiles):
                tsincr_file = os.path.join(output_dir,
                                           'tsincr_{}.npy'.format(n))
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
    log.info('process {} finished writing {} ts (incr/cuml) tifs '
             'of total {}'.format(mpiops.rank,
                                  len(process_tifs), no_ts_tifs * 2))
