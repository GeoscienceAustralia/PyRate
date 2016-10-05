import numpy as np
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
    print 'finished processing tile information'
    return ifg.shape, process_tiles, process_indices, tiles


def clean_up_old_files():
    files = glob.glob(os.path.join('out', '*.tif'))
    for f in files:
        os.remove(f)
        print 'removed', f


def save_latest_phase(d, output_dir, tiles):
    ifg = shared.Ifg(d)
    ifg.open()
    ifg.nodata_value = 0
    phase_data = ifg.phase_data
    for t in tiles:
        p_data = phase_data[t.top_left_y:t.bottom_right_y,
                 t.top_left_x:t.bottom_right_x]
        phase_file = 'phase_data_{}_{}.npy'.format(
            os.path.basename(d).split('.')[0], t.index)

        np.save(file=os.path.join(output_dir, phase_file),
                arr=p_data)
    return ifg