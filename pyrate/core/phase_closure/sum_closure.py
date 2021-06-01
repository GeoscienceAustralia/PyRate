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

import resource
from collections import namedtuple
from typing import List, Dict, Tuple, Any
from nptyping import NDArray, Float32, UInt16
import numpy as np

import pyrate.constants as C
from pyrate.core import mpiops
from pyrate.core.shared import Ifg, join_dicts
from pyrate.core.phase_closure.mst_closure import Edge, WeightedLoop
from pyrate.core.logger import pyratelogger as log

IndexedIfg = namedtuple('IndexedIfg', ['index', 'IfgPhase'])


class IfgPhase:
    """
    workaround class to only hold phase data for mpi SwigPyObject pickle error
    """

    def __init__(self, phase_data):
        self.phase_data = phase_data


def __create_ifg_edge_dict(ifg_files: List[str], params: dict) -> Dict[Edge, IndexedIfg]:
    """Returns a dictionary of indexed ifg 'edges'"""
    ifg_files.sort()
    ifgs = [Ifg(i) for i in ifg_files]

    def _func(ifg, index):
        ifg.open()
        ifg.nodata_value = params[C.NO_DATA_VALUE]
        ifg.convert_to_nans()
        ifg.convert_to_radians()
        idx_ifg = IndexedIfg(index, IfgPhase(ifg.phase_data))
        return idx_ifg

    process_ifgs = mpiops.array_split(list(enumerate(ifgs)))
    ret_combined = {}
    for idx, _ifg in process_ifgs:
        ret_combined[Edge(_ifg.first, _ifg.second)] = _func(_ifg, idx)
        _ifg.close()

    ret_combined = join_dicts(mpiops.comm.allgather(ret_combined))
    return ret_combined


def sum_phase_closures(ifg_files: List[str], loops: List[WeightedLoop], params: dict) -> \
        Tuple[NDArray[(Any, Any, Any), Float32], NDArray[(Any, Any, Any), UInt16], NDArray[(Any, ), UInt16]]:
    """
    Compute the closure sum for each pixel in each loop, and count the number of times a pixel
    contributes to a failed closure loop (where the summed closure is above/below the 
    LARGE_DEV_THR threshold).
    :param ifg_files: list of ifg files
    :param loops: list of loops
    :param params: params dict
    :return: Tuple of closure, ifgs_breach_count, num_occurrences_each_ifg
        closure: summed closure for each loop.
        ifgs_breach_count: shape=(ifg.shape, n_ifgs) number of times a pixel in an ifg fails the closure
            check (i.e., has unwrapping error) in all loops under investigation.
        num_occurrences_each_ifg: frequency of ifg appearance in all loops.
    """
    edge_to_indexed_ifgs = __create_ifg_edge_dict(ifg_files, params)
    ifgs = [v.IfgPhase for v in edge_to_indexed_ifgs.values()]
    n_ifgs = len(ifgs)
    closure_dict = {}

    if params[C.PARALLEL]:
        # rets = Parallel(n_jobs=params[cf.PROCESSES], verbose=joblib_log_level(cf.LOG_LEVEL))(
        #     delayed(__compute_ifgs_breach_count)(ifg0, n_ifgs, weighted_loop, edge_to_indexed_ifgs, params)
        #     for weighted_loop in loops
        # )
        # for k, r in enumerate(rets):
        #     closure_dict[k], ifgs_breach_count_dict[k] = r
        # TODO: enable multiprocessing - needs pickle error fix
        ifgs_breach_count = np.zeros(shape=(ifgs[0].phase_data.shape + (n_ifgs,)), dtype=np.uint16)
        for k, weighted_loop in enumerate(loops):
            closure_dict[k], ifgs_breach_count_l = __compute_ifgs_breach_count(weighted_loop, edge_to_indexed_ifgs,
                                                                               params)
            ifgs_breach_count += ifgs_breach_count_l
    else:
        loops_with_index = list(enumerate(loops))
        process_loops = mpiops.array_split(loops_with_index)
        ifgs_breach_count = np.zeros(shape=(ifgs[0].phase_data.shape + (n_ifgs,)), dtype=np.uint16)
        for k, weighted_loop in process_loops:
            closure_dict[k], ifgs_breach_count_l = __compute_ifgs_breach_count(weighted_loop, edge_to_indexed_ifgs,
                                                                               params)
            ifgs_breach_count += ifgs_breach_count_l  # process
        closure_dict = join_dicts(mpiops.comm.gather(closure_dict, root=0))

        total_gb = mpiops.comm.allreduce(ifgs_breach_count.nbytes / 1e9, op=mpiops.MPI.SUM)
        log.info("Memory usage due to ifgs_breach_count {:2.4f}GB of data".format(total_gb))
        ifgs_breach_count = mpiops.comm.reduce(ifgs_breach_count, op=mpiops.sum0_op, root=0)  # global

    closure, num_occurrences_each_ifg = None, None
    if mpiops.rank == 0:
        num_occurrences_each_ifg = _find_num_occurrences_each_ifg(loops, edge_to_indexed_ifgs, n_ifgs)
        closure = np.dstack([v for k, v in sorted(closure_dict.items(), key=lambda x: x[0])])

    return closure, ifgs_breach_count, num_occurrences_each_ifg


def _find_num_occurrences_each_ifg(loops: List[WeightedLoop],
                                   edge_to_indexed_ifgs: Dict[Edge, IndexedIfg],
                                   n_ifgs: int) -> NDArray[(Any, ), UInt16]:
    """find how many times each ifg appears in total in all loops"""
    num_occurrences_each_ifg = np.zeros(shape=n_ifgs, dtype=np.uint16)
    for weighted_loop in loops:
        for signed_edge in weighted_loop.loop:
            indexed_ifg = edge_to_indexed_ifgs[signed_edge.edge]
            ifg_index = indexed_ifg.index
            num_occurrences_each_ifg[ifg_index] += 1
    return num_occurrences_each_ifg


def __compute_ifgs_breach_count(weighted_loop: WeightedLoop,
                                edge_to_indexed_ifgs: Dict[Edge, IndexedIfg], params: dict) \
        -> Tuple[NDArray[(Any, Any), Float32], NDArray[(Any, Any, Any), UInt16]]:
    """Compute summed `closure` of each loop, and compute `ifgs_breach_count` for each pixel."""
    n_ifgs = len(edge_to_indexed_ifgs)
    indexed_ifg = list(edge_to_indexed_ifgs.values())[0]
    ifg = indexed_ifg.IfgPhase
    large_dev_thr = params[C.LARGE_DEV_THR] * np.pi
    use_median = params[C.SUBTRACT_MEDIAN]

    closure = np.zeros(shape=ifg.phase_data.shape, dtype=np.float32)
    # initiate variable for check of unwrapping issues at the same pixels in all loops
    ifgs_breach_count = np.zeros(shape=(ifg.phase_data.shape + (n_ifgs,)), dtype=np.uint16)

    for signed_edge in weighted_loop.loop:
        indexed_ifg = edge_to_indexed_ifgs[signed_edge.edge]
        ifg = indexed_ifg.IfgPhase
        closure += signed_edge.sign * ifg.phase_data
    if use_median:
        closure -= np.nanmedian(closure)  # may be able to drop median

    # this will deal with nans in `closure`, i.e., nans are not selected in indices_breaching_threshold
    indices_breaching_threshold = np.absolute(closure) > large_dev_thr

    for signed_edge in weighted_loop.loop:
        ifg_index = edge_to_indexed_ifgs[signed_edge.edge].index
        # the variable 'ifgs_breach_count' is increased by 1 for that pixel
        # make sure we are not incrementing the nan positions in the closure
        # as we don't know the phase of these pixels and also they were converted to zero before large_dev_thr check
        # Therefore, we leave them out of ifgs_breach_count, i.e., we don't increment their ifgs_breach_count values
        ifgs_breach_count[indices_breaching_threshold, ifg_index] += 1
    return closure, ifgs_breach_count
