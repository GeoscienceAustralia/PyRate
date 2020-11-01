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
THRESHOLD_TO_REMOVE_IFG = 0.2  # ifgs with more than this fraction of pixels with error will be dropped
LOOP_COUNT_FOR_THRESHOLD_TO_REMOVE_IFG = 2  # pixel with phase unwrap error in at least this many loops
PHASE_UNWRAP_ERROR_THRESHOLD = 5  # pixel with phase unwrap error in more than this many ifgs will be flagged
MAX_LOOP_LENGTH = 5  # loops upto this many edges are considered for closure checks


def detect_ps_with_unwrapping_errors(check_ps, num_occurences_each_ifg):
    nrows, ncols, n_ifgs = check_ps.shape
    ps_unwrap_error = np.zeros(shape=(nrows, ncols), dtype=np.int16)
    for i in range(n_ifgs):
        ps_idx = check_ps[:, :, i] == num_occurences_each_ifg[i]
        ps_unwrap_error[ps_idx] += 1  # number of IFGs with unwrapping errors per PS
    # PS pixels with unwrapping errors in one or more SBAS IFGs will be marked.
    # mark_ix = ps_unwrap_error > 0  # don't need to output this

    # keep_ix = ~ (ps_unwrap_error >= PHASE_UNWRAP_ERROR_THRESHOLD)
    # log.info(f'Of {nrows * ncols} pixels, {np.sum(~keep_ix)} '
    #          f'have phase unwrapping error in {PHASE_UNWRAP_ERROR_THRESHOLD} or more pixels')
    # can move mark_ix an keep_ix in wrapper if at all required
    return ps_unwrap_error


def drop_ifgs_exceeding_threshold(orig_ifg_files, check_ps, num_occurences_each_ifg):
    orig_ifg_files.sort()
    nrows, ncols, n_ifgs = check_ps.shape
    selected_ifg_files = []
    for i, ifg_file in enumerate(orig_ifg_files):
        loop_count_of_this_ifg = num_occurences_each_ifg[i]
        if loop_count_of_this_ifg:  # if the ifg participated in at least one loop
            ifg_remove_threshold_breached = np.absolute(np.nansum(check_ps[:, :, i]))/loop_count_of_this_ifg/nrows/ncols > THRESHOLD_TO_REMOVE_IFG
            print(Path(ifg_file).stem, ifg_remove_threshold_breached, num_occurences_each_ifg[i])
            if not (
                    (num_occurences_each_ifg[i] > LOOP_COUNT_FOR_THRESHOLD_TO_REMOVE_IFG)  # min loops
                    and
                    ifg_remove_threshold_breached  # and breached threshold
            ):
                selected_ifg_files.append(ifg_file)
        else:
            selected_ifg_files.append(ifg_file)

    return selected_ifg_files


def closure_check_wrapper():
    ifg_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*_unw.tif')
    ifg_files = [f.as_posix() for f in ifg_files]
    for _ in range(1):  # run closure loop twice
        print('len(ifg_files):', len(ifg_files))
        ifg_files = wrap_closure_check(ifg_files)


def wrap_closure_check(ifg_files):
    signed_loops = find_signed_closed_loops(ifg_files=ifg_files)
    print('len(signed_loops): ', len(signed_loops))
    retained_loops = [sl for sl in signed_loops if len(sl) <= MAX_LOOP_LENGTH]
    print('len(retained_loops): ', len(retained_loops))
    closure, check_ps, num_occurences_each_ifg = sum_phase_values_for_each_loop(
        ifg_files, retained_loops, LARGE_DEVIATION_THRESHOLD_FOR_PIXEL
    )
    print(np.sum(closure), np.sum(check_ps), np.sum(num_occurences_each_ifg))
    # ps_unwrap_error = detect_ps_with_unwrapping_errors(check_ps, num_occurences_each_ifg)
    selcted_ifg_files = drop_ifgs_exceeding_threshold(ifg_files, check_ps, num_occurences_each_ifg)
    return selcted_ifg_files


closure_check_wrapper()
