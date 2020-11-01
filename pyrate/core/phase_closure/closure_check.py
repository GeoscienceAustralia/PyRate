from collections import namedtuple
from pathlib import Path
from typing import List, Dict
import numpy as np
from pyrate.core.shared import Ifg, dem_or_ifg
from pyrate.core.phase_closure.mst_closure import find_signed_closed_loops
from pyrate.core.phase_closure.sum_closure import sum_phase_values_for_each_loop
from pyrate.core.logger import pyratelogger as log

THRESHOLD_TO_REMOVE_PIXEL = 0.25
LARGE_DEVIATION_THRESHOLD_FOR_PIXEL = 3.14152/4  # pi
THRESHOLD_TO_REMOVE_IFG = 0.1  # ifgs with more than this fraction of pixels with error will be dropped
PHASE_UNWRAP_ERROR_THRESHOLD = 5  # pixel with phase unwrap error in more than this many pixels will be flagged
MAX_LOOP_LENGTH = 5  # loops upto this many edges are considered for closure checks


def detect_ps_with_unwrapping_errors(check_ps, ifg_num):
    nrows, ncols, n_ifgs = check_ps.shape
    ps_unwrap_error = np.zeros(shape=(nrows, ncols), dtype=np.int16)
    for i in range(n_ifgs):
        ps_idx = check_ps[:, :, i] == ifg_num[i]  # ifg_num is zero for dropped IFGs
        ps_unwrap_error[ps_idx] += 1  # number of IFGs with unwrapping errors per PS
    # PS pixels with unwrapping errors in one or more SBAS IFGs will be marked.
    # mark_ix = ps_unwrap_error > 0  # don't need to output this

    # keep_ix = ~ (ps_unwrap_error >= PHASE_UNWRAP_ERROR_THRESHOLD)
    # log.info(f'Of {nrows * ncols} pixels, {np.sum(~keep_ix)} '
    #          f'have phase unwrapping error in {PHASE_UNWRAP_ERROR_THRESHOLD} or more pixels')
    # can move mark_ix an keep_ix in
    return ps_unwrap_error


# TODO: automated drop ifg logic

new_test_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*_unw.tif')
ifg_files = [f.as_posix() for f in new_test_files]
signed_loops = find_signed_closed_loops(ifg_files=ifg_files)
retained_loops = [sl for sl in signed_loops if len(sl) <= MAX_LOOP_LENGTH]
closure, check_ps, ifg_num = sum_phase_values_for_each_loop(ifg_files, retained_loops, LARGE_DEVIATION_THRESHOLD_FOR_PIXEL)
ps_unwrap_error = detect_ps_with_unwrapping_errors(check_ps, ifg_num)

