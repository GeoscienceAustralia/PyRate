#   This Python module is part of the PyRate software package.
#
#   Copyright 2021 Geoscience Australia
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

from collections import defaultdict
from typing import List, Dict, Tuple
import numpy as np
from pyrate.core import config as cf, mpiops
from pyrate.core.phase_closure.mst_closure import sort_loops_based_on_weights_and_date, WeightedLoop, Edge
from pyrate.configuration import Configuration
from pyrate.core.phase_closure.sum_closure import sum_phase_closures
from pyrate.core.phase_closure.plot_closure import plot_closure
from pyrate.core.shared import Ifg
from pyrate.core.logger import pyratelogger as log


def detect_ps_with_unwrapping_errors(check_ps: np.ndarray, num_occurences_each_ifg: np.ndarray, params: dict) \
        -> np.ndarray:
    """
    find where in the phase data exceed the PHS_UNW_ERR_THR, and assign nans to those pixels in all ifgs
    :param check_ps: unwrapping issues at pixels in all loops
    :param num_occurences_each_ifg:  frequency of ifgs appearing in all loops
    :param params: params dict
    :return: ps_unwrap_error: number of ifgs with unwrapping errors at each pixel
    """
    nrows, ncols, n_ifgs = check_ps.shape
    ps_unwrap_error = np.zeros(shape=(nrows, ncols), dtype=np.int16)
    for i in range(n_ifgs):
        ps_idx = check_ps[:, :, i] == num_occurences_each_ifg[i]
        ps_unwrap_error[ps_idx] += 1  # number of IFGs with unwrapping errors per PS

    # PS pixels with unwrapping errors in one or more SBAS IFGs will be marked.
    # mark_ix = ps_unwrap_error > 0  # don't need to output this

    nan_index = ps_unwrap_error >= params[cf.PHS_UNW_ERR_THR]

    log.info("Updating phase data of retained ifgs")

    for i, m_p in enumerate(params[cf.INTERFEROGRAM_FILES]):
        ifg = Ifg(m_p.tmp_sampled_path)
        ifg.open()
        ifg.nodata_value = params[cf.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.phase_data[nan_index] = np.nan
        ifg.write_modified_phase()

    log.info(f"Updated phase data of {i+1} retained ifgs after phase closure")

    # log.info(f'Of {nrows * ncols} pixels, {np.sum(~keep_ix)} '
    #          f'have phase unwrapping error in {PHS_UNW_ERR_THR} or more pixels')
    # can move mark_ix an keep_ix in wrapper if at all required
    return ps_unwrap_error


def __drop_ifgs_if_not_part_of_any_loop(ifg_files: List[str], loops: List[WeightedLoop], params: dict) -> List[str]:
    """
    Check if an ifg is part of any of the loops, otherwise drop it from the list of interferograms for further pyrate
    processing.
    """
    loop_ifgs = set()
    for weighted_loop in loops:
        for edge in weighted_loop.loop:
            loop_ifgs.add(Edge(edge.first, edge.second))

    ifgs = [Ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = params[cf.NO_DATA_VALUE]
    selected_ifg_files = []
    for i, f in zip(ifgs, ifg_files):
        if Edge(i.first, i.second) in loop_ifgs:
            selected_ifg_files.append(f)
    if len(ifg_files) != len(selected_ifg_files):
        log.info(f'Only {len(selected_ifg_files)} of the original {len(ifg_files)} '
                 f'participate in one or more loops, and selected for further pyrate analysis')
    return selected_ifg_files


def __drop_ifgs_exceeding_threshold(orig_ifg_files: List[str], check_ps, num_occurences_each_ifg, params):
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
                np.sum(check_ps[:, :, i]) / loop_count_of_this_ifg / nrows / ncols > params[cf.AVG_IFG_ERR_THR]
            if not (
                    # min loops count # check 1
                    (num_occurences_each_ifg[i] > params[cf.LOOPS_THR_IFG])
                    and
                    ifg_remove_threshold_breached  # and breached threshold
            ):
                selected_ifg_files.append(ifg_file)

    return selected_ifg_files


def filter_to_closure_checked_ifgs(config, interactive_plot=True) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """
    This function iterates to a stable list of interferogram files!
    :param config: Configuration class instance
    :param interactive_plot: bool, whether to plot sum closures of loops
    :return: stable list of ifg files, their check_ps, and number of occurrences of ifgs in loops
    """
    params = config.__dict__
    ifg_files = [ifg_path.tmp_sampled_path for ifg_path in params[cf.INTERFEROGRAM_FILES]]
    log.info(f"Performing closure check on original set of {len(ifg_files)} ifgs")

    while True:  # iterate till ifgs/loops are stable
        rets = wrap_closure_check(ifg_files, config)
        if rets is None:
            return
        new_ifg_files, closure, check_ps, num_occurences_each_ifg, loops = rets
        if interactive_plot:
            if mpiops.rank == 0:
                plot_closure(closure=closure, loops=loops, params=params, thr=params[cf.LARGE_DEV_THR])
        if len(ifg_files) == len(new_ifg_files):
            break
        else:
            ifg_files = new_ifg_files  # exit condition could be some other check like number_of_loops

    mpiops.comm.barrier()

    log.info(f"After closure check {len(ifg_files)} ifgs are retained")
    return ifg_files, check_ps, num_occurences_each_ifg


def discard_loops_containing_max_ifg_count(loops: List[WeightedLoop], params) -> List[WeightedLoop]:
    """
    This function will discard loops when each ifg participating in a loop has met the max loop count criteria.

    :param loops: list of loops
    :param params: params dict
    :return: selected loops satisfying MAX_LOOPS_IN_IFG criteria
    """
    selected_loops = []
    ifg_counter = defaultdict(int)
    for loop in loops:
        edge_appearances = np.array([ifg_counter[e] for e in loop.edges])
        if not np.all(edge_appearances > params[cf.MAX_LOOPS_IN_IFG]):
            selected_loops.append(loop)
            for e in loop.edges:
                ifg_counter[e] += 1
        else:
            log.debug(f"Loop {loop.loop} is ignored due to all it's ifgs already seen "
                      f"{params[cf.MAX_LOOPS_IN_IFG]} times or more")
    return selected_loops


def wrap_closure_check(ifg_files: List[str],  config: Configuration) -> \
        Tuple[List[str], np.ndarray, np.ndarray, np.ndarray, List[WeightedLoop]]:
    """
    This wrapper function returning the closure check outputs of a single run of closure check.

    :param ifg_files: list of ifg files
    :param config: Configuration class instance
    For return variables see docstring in `sum_phase_closures`.
    """
    params = config.__dict__
    ifg_files.sort()
    sorted_signed_loops = mpiops.run_once(sort_loops_based_on_weights_and_date, ifg_files)
    retained_loops_meeting_max_loop_criretia = [sl for sl in sorted_signed_loops
                                                if len(sl) <= params[cf.MAX_LOOP_LENGTH]]
    msg = f"After applying MAX_LOOP_LENGTH={params[cf.MAX_LOOP_LENGTH]} criteria, " \
          f"{len(retained_loops_meeting_max_loop_criretia)} loops are retained"

    if len(retained_loops_meeting_max_loop_criretia) < 1:
        return None
    else:
        log.info(msg)

    retained_loops = mpiops.run_once(discard_loops_containing_max_ifg_count,
                                     retained_loops_meeting_max_loop_criretia, params)
    ifgs_with_loops = mpiops.run_once(__drop_ifgs_if_not_part_of_any_loop, ifg_files, retained_loops, params)

    msg = f"After applying MAX_LOOPS_IN_IFG={params[cf.MAX_LOOPS_IN_IFG]} criteria, " \
          f"{len(retained_loops)} loops are retained"
    if len(retained_loops) < 1:
        return None
    else:
        log.info(msg)

    closure, check_ps, num_occurences_each_ifg = sum_phase_closures(ifgs_with_loops, retained_loops, params)

    if mpiops.rank == 0:
        closure_ins = config.closure()
        np.save(closure_ins.closure, closure)
        np.save(closure_ins.check_ps, check_ps)
        np.save(closure_ins.num_occurences_each_ifg, num_occurences_each_ifg)
        np.save(closure_ins.loops, retained_loops, allow_pickle=True)

    selected_ifg_files = mpiops.run_once(__drop_ifgs_exceeding_threshold,
                                         ifgs_with_loops, check_ps, num_occurences_each_ifg, params)
    return selected_ifg_files, closure, check_ps, num_occurences_each_ifg, retained_loops
