from __future__ import print_function
import logging
import numpy as np
from operator import itemgetter
import glob
import os
from os.path import basename, join

from pyrate import shared
from pyrate import config as cf
from pyrate import mpiops


log = logging.getLogger(__name__)


def get_process_tiles(dest_tifs, parallel, params):
    ifg = shared.pre_prepare_ifgs([dest_tifs[0]], params)[0]
    # TODO: changing the nrows and ncols here will break tests, fix this
    tiles = shared.create_tiles(ifg.shape, nrows=3, ncols=4)
    no_tiles = len(tiles)
    process_indices = parallel.calc_indices(no_tiles)
    process_tiles = [itemgetter(p)(tiles) for p in process_indices]
    print('finished processing tile information')
    return ifg.shape, process_tiles, process_indices, tiles


def clean_up_old_files():
    files = glob.glob(join('out', '*.tif'))
    for f in files:
        os.remove(f)
        log.info('removed {}'.format(f))


def prepare_ifg(ifg_path, params):
    """
    Parameters
    ----------
    ifg_path: str
        ifg path
    tiles: list
        list of Tile instances
    params: dict
        config dict
    Returns
    -------
    ifg: ifg class instance
    """
    ifg = shared.Ifg(ifg_path)
    ifg.open()
    shared.nan_and_mm_convert(ifg, params)
    return ifg


def save_numpy_phase(ifg_paths, tiles, params):
    """
    :param ifg_paths:
    :param params:
    :param tiles:
    :return:
    """
    process_ifgs = mpiops.array_split(ifg_paths)
    outdir = params[cf.OUT_DIR]
    for ifg_path in process_ifgs:
        ifg = shared.Ifg(ifg_path)
        ifg.open()
        phase_data = ifg.phase_data
        bname = basename(ifg_path).split('.')[0]
        for t in tiles:
            p_data = phase_data[
                     t.top_left_y:t.bottom_right_y,
                     t.top_left_x:t.bottom_right_x
                     ]
            phase_file = 'phase_data_{}_{}.npy'.format(bname, t.index)
            np.save(file=join(outdir, phase_file),
                    arr=p_data)
        ifg.close()
    mpiops.comm.barrier()
