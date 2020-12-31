from collections import defaultdict
from typing import List
import numpy as np
from pyrate.core.shared import dem_or_ifg
from pyrate.core import config as cf
from pyrate.core.phase_closure.mst_closure import find_signed_closed_loops, sort_loops_based_on_weights_and_date, \
    WeightedLoop, Edge
from pyrate.core.phase_closure.sum_closure import sum_phase_values_for_each_loop
from pyrate.core.phase_closure.plot_closure import plot_closure
from pyrate.core.logger import pyratelogger as log


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


def drop_ifgs_if_not_part_of_any_loop(ifg_files: List[str], loops: List[WeightedLoop]) -> List[str]:
    loop_ifgs = set()
    for weighted_loop in loops:
        for edge in weighted_loop.loop:
            loop_ifgs.add(Edge(edge.first, edge.second))

    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    selected_ifg_files = []
    for i, f in zip(ifgs, ifg_files):
        if Edge(i.first, i.second) in loop_ifgs:
            selected_ifg_files.append(f)
    if len(ifg_files) != len(selected_ifg_files):
        log.info(f'Only {len(selected_ifg_files)} of the original {len(ifg_files)} '
                 f'participate in one or more loops, and selected for further pyrate analysis')
    return selected_ifg_files


def drop_ifgs_exceeding_threshold(orig_ifg_files, check_ps, num_occurences_each_ifg, params):
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
            ifg_remove_threshold_breached = \
                np.sum(check_ps[:, :, i])/loop_count_of_this_ifg/nrows/ncols > params[cf.THRESHOLD_TO_REMOVE_IFG]
            if not (
                # min loops count # check 1
                (num_occurences_each_ifg[i] > params[cf.LOOP_COUNT_FOR_THRESHOLD_TO_REMOVE_IFG])
                and
                ifg_remove_threshold_breached  # and breached threshold
            ):
                selected_ifg_files.append(ifg_file)
        else:
            selected_ifg_files.append(ifg_file)

    return selected_ifg_files


def filter_to_closure_checked_ifgs(params, interactive_plot=True):
    ifg_files = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    log.info(f"Performing closure check on original set of {len(ifg_files)} ifgs")

    while True:  # iterate till ifgs/loops are stable
        rets = wrap_closure_check(ifg_files, params)
        if rets is None:
            return
        new_ifg_files, closure, loops = rets
        if interactive_plot:
            plot_closure(closure=closure, loops=loops, params=params, thr=params[cf.LARGE_DEV_THR])
        if len(ifg_files) == len(new_ifg_files):
            break
        else:
            ifg_files = new_ifg_files  # exit condition could be some other check like number_of_loops

    log.info(f"After closure check {len(ifg_files)} ifgs are retained")
    return ifg_files


def discard_loops_containing_max_ifg_count(loops: List[WeightedLoop], params) -> List[WeightedLoop]:
    # available_edges = setup_edges(ifg_files)
    selected_loops = []
    ifg_counter = defaultdict(int)
    for l in loops:
        edge_apperances = np.array([ifg_counter[e] for e in l.edges])
        if not np.all(edge_apperances > params[cf.MAX_LOOP_COUNT_FOR_EACH_IFGS]):
            selected_loops.append(l)
            for e in l.edges:
                ifg_counter[e] += 1
        else:
            log.debug(f"Loop {l.loop} is ignored due to all it's ifgs already seen "
                     f"{params[cf.MAX_LOOP_COUNT_FOR_EACH_IFGS]} times or more")
    return selected_loops


def wrap_closure_check(ifg_files, params):
    signed_loops = find_signed_closed_loops(ifg_files=ifg_files)
    sorted_signed_loops = sort_loops_based_on_weights_and_date(signed_loops)
    retained_loops_meeting_max_loop_criretia = [sl for sl in sorted_signed_loops
                                                if len(sl) <= params[cf.MAX_LOOP_LENGTH]]
    msg = f"After applying MAX_LOOP_LENGTH={params[cf.MAX_LOOP_LENGTH]} criteria, " \
          f"{len(retained_loops_meeting_max_loop_criretia)} loops are retained"

    if len(retained_loops_meeting_max_loop_criretia) < 1:
        return None
    else:
        log.info(msg)

    retained_loops = discard_loops_containing_max_ifg_count(retained_loops_meeting_max_loop_criretia, params)
    ifgs_with_loops = drop_ifgs_if_not_part_of_any_loop(ifg_files, retained_loops)

    msg = f"After applying MAX_LOOP_COUNT_FOR_EACH_IFGS={params[cf.MAX_LOOP_COUNT_FOR_EACH_IFGS]} criteria, " \
             f"{len(retained_loops)} loops are retained"
    if len(retained_loops) < 1:
        return None
    else:
        log.info(msg)

    closure, check_ps, num_occurences_each_ifg = sum_phase_values_for_each_loop(
        ifgs_with_loops, retained_loops, params[cf.LARGE_DEV_THR],
        params[cf.SUBTRACT_MEDIAN_IN_CLOSURE_CHECK]
    )
    # ps_unwrap_error = detect_ps_with_unwrapping_errors(check_ps, num_occurences_each_ifg)
    selcted_ifg_files = drop_ifgs_exceeding_threshold(ifgs_with_loops, check_ps, num_occurences_each_ifg, params)
    return selcted_ifg_files, closure, retained_loops

