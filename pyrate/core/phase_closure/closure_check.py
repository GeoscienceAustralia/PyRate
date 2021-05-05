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
from typing import List, Dict, Tuple, Any
from nptyping import NDArray, UInt16, Float32
import numpy as np

import pyrate.constants as C
from pyrate.core import mpiops
from pyrate.core.phase_closure.mst_closure import sort_loops_based_on_weights_and_date, WeightedLoop, Edge
from pyrate.configuration import Configuration, MultiplePaths
from pyrate.core.phase_closure.sum_closure import sum_phase_closures
from pyrate.core.phase_closure.plot_closure import plot_closure
from pyrate.core.shared import Ifg
from pyrate.core.logger import pyratelogger as log


def mask_pixels_with_unwrapping_errors(ifgs_breach_count: NDArray[(Any, Any, Any), UInt16],
                                     num_occurrences_each_ifg: NDArray[(Any,), UInt16],
                                     params: dict) -> NDArray[(Any, Any), UInt16]:
    """
    Find pixels in the phase data exceeding the PHS_UNW_ERR_THR, and mask
    (assign nans) to those pixels in all ifgs
    :param ifgs_breach_count: unwrapping issues at pixels in all loops
    :param num_occurrences_each_ifg:  frequency of ifgs appearing in all loops
    :param params: params dict
    :return: pix_unwrap_error: number of ifgs with unwrapping errors at each pixel
    """
    nrows, ncols, n_ifgs = ifgs_breach_count.shape
    pix_unwrap_error = np.zeros(shape=(nrows, ncols), dtype=np.uint16)
    for i in range(n_ifgs):
        pix_idx = ifgs_breach_count[:, :, i] == num_occurrences_each_ifg[i]
        pix_unwrap_error[pix_idx] += 1  # number of IFGs with unwrapping errors per pixel

    # Pixels with unwrapping errors in one or more SBAS IFGs will be marked.
    # mark_ix = pix_unwrap_error > 0  # don't need to output this

    nan_index = pix_unwrap_error >= params[C.PHS_UNW_ERR_THR]

    log.debug("Updating phase data of retained ifgs")

    for i, m_p in enumerate(params[C.INTERFEROGRAM_FILES]):
        ifg = Ifg(m_p.tmp_sampled_path)
        ifg.open()
        ifg.nodata_value = params[C.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.phase_data[nan_index] = np.nan
        ifg.write_modified_phase()

    log.info(f"Updated phase data of {i + 1} retained ifgs after phase closure")

    # log.info(f'Of {nrows * ncols} pixels, {np.sum(~keep_ix)} '
    #          f'have phase unwrapping error in {PHS_UNW_ERR_THR} or more pixels')
    # can move mark_ix an keep_ix in wrapper if at all required
    return pix_unwrap_error


def __drop_ifgs_if_not_part_of_any_loop(ifg_files: List[str], loops: List[WeightedLoop], params: dict) -> List[str]:
    """
    Check if an ifg is part of any of the loops, otherwise drop it from the list of interferograms for further PyRate
    processing.
    """
    loop_ifgs = set()
    for weighted_loop in loops:
        for edge in weighted_loop.loop:
            loop_ifgs.add(Edge(edge.first, edge.second))

    ifgs = [Ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = params[C.NO_DATA_VALUE]
    selected_ifg_files = []
    for i, f in zip(ifgs, ifg_files):
        if Edge(i.first, i.second) in loop_ifgs:
            selected_ifg_files.append(f)
    if len(ifg_files) != len(selected_ifg_files):
        log.info(f'Only {len(selected_ifg_files)} (out of {len(ifg_files)}) ifgs participate in '
                 f'one or more closure loops, and are selected for further PyRate analysis')
    return selected_ifg_files


def __drop_ifgs_exceeding_threshold(orig_ifg_files: List[str], ifgs_breach_count, num_occurences_each_ifg, params):
    """
    We demand two thresholds breaches for an ifg to be dropped.
    1. The first one is the basic ifg loop participation count check.
    2. The second threshold check is a weighted average check of pixels breached taking all loops into account.
        (a) ifgs_breach_count contains unwrapping error count for each pixel for each ifg seen in all loops
        (b) sum(ifgs_breach_count[:, :, i]) is pixel total count with unwrapping error for i-th ifg over all loops
        (c) divide by loop_count_of_this_ifg and num of cells (nrows x ncols) for a weighted measure of threshold

    """
    orig_ifg_files.sort()
    nrows, ncols, n_ifgs = ifgs_breach_count.shape
    selected_ifg_files = []
    for i, ifg_file in enumerate(orig_ifg_files):
        loop_count_of_this_ifg = num_occurences_each_ifg[i]
        if loop_count_of_this_ifg:  # if the ifg participated in at least one loop
            ifg_remove_threshold_breached = \
                np.sum(ifgs_breach_count[:, :, i]) / loop_count_of_this_ifg / nrows / ncols > params[
                    C.AVG_IFG_ERR_THR]
            if not (
                    # min loops count # check 1
                    (num_occurences_each_ifg[i] > params[C.LOOPS_THR_IFG])
                    and
                    ifg_remove_threshold_breached  # and breached threshold
            ):
                selected_ifg_files.append(ifg_file)

    return selected_ifg_files


def iterative_closure_check(config, interactive_plot=True) -> \
        Tuple[List[str], NDArray[(Any, Any, Any), UInt16], NDArray[(Any,), UInt16]]:
    """
    This function iterates the closure check until a stable list of interferogram files is returned.
    :param config: Configuration class instance
    :param interactive_plot: bool, whether to plot sum closures of loops
    :return: stable list of ifg files, their ifgs_breach_count, and number of occurrences of ifgs in loops
    """
    params = config.__dict__
    ifg_files = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    i = 1  # iteration counter

    while True:  # iterate till ifgs/loops are stable
        log.info(f"Closure check iteration #{i}: working on {len(ifg_files)} ifgs")
        rets = __wrap_closure_check(config)
        if rets is None:
            return
        new_ifg_files, closure, ifgs_breach_count, num_occurences_each_ifg, loops = rets
        if interactive_plot:
            if mpiops.rank == 0:
                plot_closure(closure=closure, loops=loops, config=config, 
                                thr=params[C.LARGE_DEV_THR], iteration=i)
        if len(ifg_files) == len(new_ifg_files):
            break
        else:
            i += 1
            ifg_files = new_ifg_files  # exit condition could be some other check like number_of_loops

    mpiops.comm.barrier()

    log.info(f"Stable list of ifgs achieved after iteration #{i}. {len(ifg_files)} ifgs are retained")
    return ifg_files, ifgs_breach_count, num_occurences_each_ifg


def discard_loops_containing_max_ifg_count(loops: List[WeightedLoop], params) -> List[WeightedLoop]:
    """
    This function will discard loops when each ifg participating in a loop has met the max loop count criteria.

    :param loops: list of loops
    :param params: params dict
    :return: selected loops satisfying MAX_LOOPS_PER_IFG criteria
    """
    selected_loops = []
    ifg_counter = defaultdict(int)
    for loop in loops:
        edge_appearances = np.array([ifg_counter[e] for e in loop.edges])
        if not np.all(edge_appearances > params[C.MAX_LOOPS_PER_IFG]):
            selected_loops.append(loop)
            for e in loop.edges:
                ifg_counter[e] += 1
        else:
            log.debug(f"Loop {loop.loop} ignored: all constituent ifgs have been in a loop "
                      f"{params[C.MAX_LOOPS_PER_IFG]} times or more")
    return selected_loops


def __wrap_closure_check(config: Configuration) -> \
        Tuple[
            List[str],
            NDArray[(Any, Any), Float32],
            NDArray[(Any, Any, Any), UInt16],
            NDArray[(Any,), UInt16],
            List[WeightedLoop]]:
    """
    This wrapper function returns the closure check outputs for a single iteration of closure check.

    :param config: Configuration class instance
    For return variables see docstring in `sum_phase_closures`.
    """
    params = config.__dict__
    ifg_files = [ifg_path.tmp_sampled_path for ifg_path in params[C.INTERFEROGRAM_FILES]]
    ifg_files.sort()
    log.debug(f"The number of ifgs in the list is {len(ifg_files)}")
    sorted_signed_loops = mpiops.run_once(sort_loops_based_on_weights_and_date, params)
    log.info(f"Total number of selected closed loops with up to MAX_LOOP_LENGTH = "
             f"{params[C.MAX_LOOP_LENGTH]} edges is {len(sorted_signed_loops)}")

    if len(sorted_signed_loops) < 1:
        return None

    retained_loops = mpiops.run_once(discard_loops_containing_max_ifg_count,
                                     sorted_signed_loops, params)
    ifgs_with_loops = mpiops.run_once(__drop_ifgs_if_not_part_of_any_loop, ifg_files, retained_loops, params)

    msg = f"After applying MAX_LOOPS_PER_IFG = {params[C.MAX_LOOPS_PER_IFG]} criteria, " \
          f"{len(retained_loops)} loops are retained"
    if len(retained_loops) < 1:
        return None
    else:
        log.info(msg)

    closure, ifgs_breach_count, num_occurences_each_ifg = sum_phase_closures(ifgs_with_loops, retained_loops, params)

    if mpiops.rank == 0:
        closure_ins = config.closure()
        np.save(closure_ins.closure, closure)
        np.save(closure_ins.ifgs_breach_count, ifgs_breach_count)
        np.save(closure_ins.num_occurences_each_ifg, num_occurences_each_ifg)
        np.save(closure_ins.loops, retained_loops, allow_pickle=True)

    selected_ifg_files = mpiops.run_once(__drop_ifgs_exceeding_threshold,
                                         ifgs_with_loops, ifgs_breach_count, num_occurences_each_ifg, params)

    # update the ifg list in the parameters dictionary
    params[C.INTERFEROGRAM_FILES] = \
        mpiops.run_once(update_ifg_list, selected_ifg_files, params[C.INTERFEROGRAM_FILES])
    return selected_ifg_files, closure, ifgs_breach_count, num_occurences_each_ifg, retained_loops

# TODO: Consider whether this helper function is better homed in a generic module and used more widely
def update_ifg_list(ifg_files: List[str], multi_paths: List[MultiplePaths]) -> List[MultiplePaths]:
    """
    Function to extract full paths for a subsetted list of interferograms
    :param ifg_files: list of interferograms to subset
    :param multi_paths: list of full paths for list of original interferograms
    :return: filtered_multi_paths: list of full paths for the subset
    """
    filtered_multi_paths = []
    for m_p in multi_paths:
        if m_p.tmp_sampled_path in ifg_files:
            filtered_multi_paths.append(m_p)
    return filtered_multi_paths
