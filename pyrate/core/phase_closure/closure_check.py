from collections import namedtuple
from pathlib import Path
from typing import List, Dict
import numpy as np
from pyrate.core.shared import Ifg, dem_or_ifg
from pyrate.core.phase_closure.mst_closure import find_signed_closed_loops
from pyrate.core.phase_closure.sum_closure import sum_phase_values_for_each_loop
from pyrate.core.logger import pyratelogger as log

LARGE_DEVIATION_THRESHOLD_FOR_PIXEL = np.pi/2  # pi
THRESHOLD_TO_REMOVE_IFG = 0.1  # ifgs with more than this fraction of pixels with error will be dropped
LOOP_COUNT_FOR_THRESHOLD_TO_REMOVE_IFG = 2  # pixel with phase unwrap error in at least this many loops
PHASE_UNWRAP_ERROR_THRESHOLD = 5  # pixel with phase unwrap error in more than this many ifgs will be flagged
MAX_LOOP_LENGTH = 4  # loops upto this many edges are considered for closure checks


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
    """
    We demand two thresholds breaches for an ifg to be dropped.
    1. The first one is the basic ifg loop participation count check.
    2. The second threshold check is a weighted average check of pixels breached taking all loops into account.
        (a) check_ps contains unwrapping error count for each pixel for each ifg seen in any loop
        (b) sum(check_ps[:, :, i]) is pixel total count with unwrapping error for i-th ifg over all loops
        (c) divide by loop_count_of_this_ifg and num of cells (nrows x ncols) for a weighted measure of threshold
    """
    orig_ifg_files.sort()
    nrows, ncols, n_ifgs = check_ps.shape
    selected_ifg_files = []
    for i, ifg_file in enumerate(orig_ifg_files):
        loop_count_of_this_ifg = num_occurences_each_ifg[i]
        if loop_count_of_this_ifg:  # if the ifg participated in at least one loop
            ifg_remove_threshold_breached = np.sum(check_ps[:, :, i])/loop_count_of_this_ifg/nrows/ncols > THRESHOLD_TO_REMOVE_IFG
            if not (
                    (num_occurences_each_ifg[i] > LOOP_COUNT_FOR_THRESHOLD_TO_REMOVE_IFG)  # min loops count # check 1
                    and
                    ifg_remove_threshold_breached  # and breached threshold
            ):
                selected_ifg_files.append(ifg_file)
        else:
            selected_ifg_files.append(ifg_file)

    return selected_ifg_files


def closure_check_wrapper():
    from pyrate.core.phase_closure.plot_closure import plot_closure
    from pyrate.constants import PYRATEPATH
    GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')
    ifg_files = GEOTIFF.glob('*_unw.tif')
    ifg_files = [f.as_posix() for f in ifg_files]
    while True:  # iterate till ifgs/loops are stable
        print('len(ifg_files):', len(ifg_files))
        new_ifg_files, closure = wrap_closure_check(ifg_files)
        plot_closure(closure=closure)
        if len(ifg_files) == len(new_ifg_files):
            break
        else:
            ifg_files = new_ifg_files  # exit condition could be some other check like number_of_loops


def wrap_closure_check(ifg_files):
    signed_loops = find_signed_closed_loops(ifg_files=ifg_files)
    retained_loops = [sl for sl in signed_loops if len(sl) <= MAX_LOOP_LENGTH]
    closure, check_ps, num_occurences_each_ifg = sum_phase_values_for_each_loop(
        ifg_files, retained_loops, LARGE_DEVIATION_THRESHOLD_FOR_PIXEL
    )
    # ps_unwrap_error = detect_ps_with_unwrapping_errors(check_ps, num_occurences_each_ifg)
    selcted_ifg_files = drop_ifgs_exceeding_threshold(ifg_files, check_ps, num_occurences_each_ifg)
    return selcted_ifg_files, closure


closure_check_wrapper()
