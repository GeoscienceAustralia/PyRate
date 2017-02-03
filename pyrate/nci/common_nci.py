from __future__ import print_function
import logging
import numpy as np
from operator import itemgetter
import glob
import os
from os.path import basename, join

from pyrate import shared
from pyrate import config as cf


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


def save_latest_phase(ifg_path, tiles, params):
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
    phase_data = ifg.phase_data

    # TODO: remove saving whole ifg phase data when possible
    numpy_file = os.path.join(
        params[cf.OUT_DIR], os.path.basename(ifg_path).split('.')[0] + '.npy')
    np.save(file=numpy_file, arr=phase_data)
    for t in tiles:
        p_data = phase_data[
                 t.top_left_y:t.bottom_right_y,
                 t.top_left_x:t.bottom_right_x
                 ]
        phase_file = 'phase_data_{}_{}.npy'.format(
            basename(ifg_path).split('.')[0], t.index)

        np.save(file=join(params[cf.OUT_DIR], phase_file),
                arr=p_data)
    return ifg