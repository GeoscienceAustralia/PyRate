from collections import namedtuple
from typing import List, Dict, Tuple
from joblib import Parallel, delayed
import numpy as np
from pyrate.core import config as cf, mpiops
from pyrate.core.shared import Ifg, dem_or_ifg, join_dicts, joblib_log_level
from pyrate.core.phase_closure.mst_closure import Edge, SignedEdge, WeightedLoop
IndexedIfg = namedtuple('IndexedIfg', ['index', 'Ifg'])


def create_ifg_edge_dict(ifg_files: List[str]) -> Dict[Edge, IndexedIfg]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
        i.convert_to_nans()
    return {Edge(ifg.first, ifg.second): IndexedIfg(index, ifg) for index, ifg in enumerate(ifgs)}


def sum_phase_values_for_each_loop(ifg_files: List[str], loops: List[WeightedLoop], params: dict) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
    edge_to_indexed_ifgs = create_ifg_edge_dict(ifg_files)

    ifgs = [v.Ifg for v in edge_to_indexed_ifgs.values()]
    n_ifgs = len(ifgs)
    ifg0 = ifgs[0]

    closure_dict = {}
    check_ps_dict = {}
    if params[cf.PARALLEL]:
        # rets = Parallel(n_jobs=params[cf.PROCESSES], verbose=joblib_log_level(cf.LOG_LEVEL))(
        #     delayed(__compute_check_ps)(ifg0, n_ifgs, weighted_loop, edge_to_indexed_ifgs, params)
        #     for weighted_loop in loops
        # )
        # for k, r in enumerate(rets):
        #     closure_dict[k], check_ps_dict[k] = r
        # TODO: enable multiprocessing - needs pickle error fix
        for k, weighted_loop in enumerate(loops):
            closure_dict[k], check_ps_dict[k] = __compute_check_ps(
                ifg0, n_ifgs, weighted_loop, edge_to_indexed_ifgs, params
            )
    else:
        loops_with_index = list(enumerate(loops))
        process_loops = mpiops.array_split(loops_with_index)
        for k, weighted_loop in process_loops:
            closure_dict[k], check_ps_dict[k] = __compute_check_ps(
                ifg0, n_ifgs, weighted_loop, edge_to_indexed_ifgs, params
            )
        closure_dict = join_dicts(mpiops.comm.gather(closure_dict, root=0))
        check_ps_dict = join_dicts(mpiops.comm.gather(check_ps_dict, root=0))
    if mpiops.rank == 0:
        num_occurences_each_ifg = _find_num_occurences_each_ifg(loops, edge_to_indexed_ifgs, n_ifgs)
        closure = np.dstack(tuple(closure_dict.values()))
        check_ps = np.sum(np.stack(tuple(check_ps_dict.values()), axis=3), axis=3)
    else:
        closure, check_ps, num_occurences_each_ifg = None, None, None

    return closure, check_ps, num_occurences_each_ifg


def _find_num_occurences_each_ifg(loops, edge_to_indexed_ifgs, n_ifgs):
    """find how many times each ifg appears in total in all loops"""
    num_occurences_each_ifg = np.zeros(shape=n_ifgs, dtype=np.uint16)
    for weighted_loop in loops:
        for signed_edge in weighted_loop.loop:
            indexed_ifg = edge_to_indexed_ifgs[signed_edge.edge]
            ifg_index = indexed_ifg.index
            num_occurences_each_ifg[ifg_index] += 1
    return num_occurences_each_ifg


def __compute_check_ps(ifg: Ifg, n_ifgs: int, weighted_loop: WeightedLoop,
                       edge_to_indexed_ifgs: Dict[Edge, IndexedIfg], params: dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    find sum `closure` of each loop, and compute `check_ps` for each pixel.
    PS: Persistent Scatterer
    """
    large_dev_thr = params[cf.LARGE_DEV_THR] * cf
    use_median = params[cf.SUBTRACT_MEDIAN_IN_CLOSURE_CHECK]
    closure = np.zeros(shape=ifg.phase_data.shape, dtype=np.float32)
    # initiate variable for check of unwrapping issues at the same pixels in all loops
    check_ps = np.zeros(shape=(ifg.phase_data.shape + (n_ifgs,)), dtype=np.uint16)

    for signed_edge in weighted_loop.loop:
        indexed_ifg = edge_to_indexed_ifgs[signed_edge.edge]
        ifg = indexed_ifg.Ifg
        closure += signed_edge.sign * ifg.phase_data
    if use_median:
        closure -= np.nanmedian(closure)  # may be able to drop median
    # handle nans elegantly
    nan_indices = np.isnan(closure)
    closure[nan_indices] = 0  # values with nans can't be large_dev_thr checked
    indices_breaching_threshold = np.absolute(closure) > large_dev_thr
    closure[nan_indices] = np.nan  # set them to nan again  - this is useful when we plot
    for signed_edge in weighted_loop.loop:
        ifg_index = edge_to_indexed_ifgs[signed_edge.edge].index
        #  the variable check_ps is increased by 1 for that pixel
        # make sure we are not incrementing the nan positions in the closure
        # as we don't know the PS of these pixels and also they were converted to zero before large_dev_thr check
        # Therefore, we leave them out of check_ps, i.e., we don't increment their check_ps values
        check_ps[np.logical_and(indices_breaching_threshold, ~nan_indices), ifg_index] += 1
    return closure, check_ps

