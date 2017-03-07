#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module does post-processing steps to assemble the linear
rate and time series outputs and save as geotiff files
"""
import os
from os.path import join
import logging
import pickle as cp
import numpy as np
from osgeo import gdal

from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate import shared
from pyrate.scripts import run_pyrate
from pyrate import mpiops
from pyrate.shared import PrereadIfg
gdal.SetCacheMax(64)

# pylint: disable=logging-format-interpolation
log = logging.getLogger(__name__)

# Constants
MASTER_PROCESS = 0


def main(config_file, rows, cols):
    """ postprocessing main func"""
    # setup paths
    _, _, params = cf.get_ifg_paths(config_file)
    postprocess_linrate(rows, cols, params)
    if params[cf.TIME_SERIES_CAL]:
        postprocess_timeseries(rows, cols, params)


def postprocess_linrate(rows, cols, params):
    """ Postprocess linear rate """
    # pylint: disable=expression-not-assigned
    # setup paths
    xlks, _, crop = cf.transform_params(params)
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
    """ Save linear rate outputs"""
    log.info('Starting PyRate postprocessing {}'.format(out_type))
    gt, md, wkt = ifgs_dict['gt'], ifgs_dict['md'], ifgs_dict['wkt']
    epochlist = ifgs_dict['epochlist']
    ifgs = [v for v in ifgs_dict.values() if isinstance(v, PrereadIfg)]
    dest = os.path.join(params[cf.OUT_DIR], out_type + ".tif")
    md[ifc.EPOCH_DATE] = epochlist.dates
    if out_type == 'linrate':
        md[ifc.DATA_TYPE] = ifc.LINRATE
    elif out_type == 'linerror':
        md[ifc.DATA_TYPE] = ifc.LINERROR
    else:
        md[ifc.DATA_TYPE] = ifc.LINSAMP

    rate = np.zeros(shape=ifgs[0].shape, dtype=np.float32)
    for t in tiles:
        rate_file = os.path.join(params[cf.OUT_DIR], out_type +
                                 '_{}.npy'.format(t.index))
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y,
             t.top_left_x:t.bottom_right_x] = rate_tile
    shared.write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    npy_rate_file = os.path.join(params[cf.OUT_DIR], out_type + '.npy')
    np.save(file=npy_rate_file, arr=rate)
    log.info('Finished PyRate postprocessing {}'.format(out_type))


def postprocess_timeseries(rows, cols, params):
    """ Postprocess time series output """
    # pylint: disable=too-many-locals
    xlks, _, crop = cf.transform_params(params)
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

    # pylint: disable=no-member
    no_ts_tifs = tsincr.shape[2]
    # we create 2 x no_ts_tifs as we are splitting tsincr and tscuml
    # to all processes.
    process_tifs = mpiops.array_split(range(2 * no_ts_tifs))

    # depending on nvelpar, this will not fit in memory
    # e.g. nvelpar=100, nrows=10000, ncols=10000, 32bit floats need 40GB memory
    # 32 * 100 * 10000 * 10000 / 8 bytes = 4e10 bytes = 40 GB
    # the double for loop helps us overcome the memory limit
    log.info('process {} will write {} ts (incr/cuml) tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))
    for i in process_tifs:
        tscum_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
        if i < no_ts_tifs:
            for n, t in enumerate(tiles):
                tscum_file = os.path.join(output_dir,
                                          'tscuml_{}.npy'.format(n))
                tscum = np.load(file=tscum_file)

                md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
                # sequence position; first time slice is #0
                md['SEQUENCE_POSITION'] = i+1
                tscum_g[t.top_left_y:t.bottom_right_y,
                        t.top_left_x:t.bottom_right_x] = tscum[:, :, i]
                dest = os.path.join(params[cf.OUT_DIR],
                                    'tscuml' + "_" +
                                    str(epochlist.dates[i + 1]) + ".tif")
                md[ifc.DATA_TYPE] = ifc.CUML
                shared.write_output_geotiff(md, gt, wkt, tscum_g, dest, np.nan)
        else:
            tsincr_g = np.empty(shape=ifgs[0].shape, dtype=np.float32)
            i %= no_ts_tifs
            for n, t in enumerate(tiles):
                tsincr_file = os.path.join(output_dir,
                                           'tsincr_{}.npy'.format(n))
                tsincr = np.load(file=tsincr_file)

                md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
                # sequence position; first time slice is #0
                md['SEQUENCE_POSITION'] = i+1
                tsincr_g[t.top_left_y:t.bottom_right_y,
                         t.top_left_x:t.bottom_right_x] = tsincr[:, :, i]
                dest = os.path.join(params[cf.OUT_DIR],
                                    'tsincr' + "_" + str(
                                        epochlist.dates[i + 1]) + ".tif")
                md[ifc.DATA_TYPE] = ifc.INCR
                shared.write_output_geotiff(md, gt, wkt, tsincr_g, dest,
                                            np.nan)
    log.info('process {} finished writing {} ts (incr/cuml) tifs of '
             'total {}'.format(mpiops.rank, len(process_tifs), no_ts_tifs * 2))
