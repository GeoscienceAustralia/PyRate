from pyrate import shared
from operator import itemgetter
import glob
import os


def get_process_tiles(dest_tifs, parallel, params):
    ifg = shared.pre_prepare_ifgs([dest_tifs[0]], params)[0]
    # TODO: changing the nrows and ncols here will break tests, fix this
    tiles = shared.create_tiles(ifg.shape,
                                n_ifgs=len(dest_tifs), nrows=3, ncols=4)
    no_tiles = len(tiles)
    process_indices = parallel.calc_indices(no_tiles)
    process_tiles = [itemgetter(p)(tiles) for p in process_indices]
    return ifg.shape, process_tiles, process_indices, tiles


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f
