from collections import namedtuple
from typing import List, Dict, Tuple
import numpy as np
from pyrate.core.shared import Ifg, dem_or_ifg
from pyrate.core.phase_closure.mst_closure import Edge, SignedEdge
IndexedIfg = namedtuple('IndexedIfg', ['index', 'Ifg'])


def create_ifg_edge_dict(ifg_files) -> Dict[Edge, IndexedIfg]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
        i.convert_to_nans()
    return {Edge(ifg.first, ifg.second): IndexedIfg(index, ifg) for index, ifg in enumerate(ifgs)}


def sum_phase_values_for_each_loop(ifg_files: List[str], loops: List[List[SignedEdge]], threshold: float):
    edge_to_indexed_ifgs = create_ifg_edge_dict(ifg_files)
    ifgs = [v.Ifg for v in edge_to_indexed_ifgs.values()]
    n_ifgs = len(ifgs)
    n_loops = len(loops)
    closure = np.zeros(shape=(ifgs[0].phase_data.shape + (n_loops,)), dtype=np.float32)

    num_occurences_each_ifg = np.zeros(shape=n_ifgs, dtype=np.int16)

    # initiate variable for check of unwrapping issues at the same pixels in all loops
    check_ps = np.zeros(shape=(ifgs[0].phase_data.shape + (n_ifgs,)))
    for k, loop in enumerate(loops):
        for signed_edge in loop:
            indexed_ifg = edge_to_indexed_ifgs[signed_edge.edge]
            ifg = indexed_ifg.Ifg
            ifg_index = indexed_ifg.index
            closure[:, :, k] += signed_edge.sign * ifg.phase_data
            num_occurences_each_ifg[ifg_index] += 1

        # closure[:, :, k] -= np.nanmedian(closure[:, :, k])  # may be able to drop median
        # handle nans elegantly
        nan_indices = np.isnan(closure[:, :, k])
        closure[:, :, k][nan_indices] = 0  # values with nans can't be threshold checked

        indices_breaching_threshold = np.absolute(closure[:, :, k]) > threshold
        closure[:, :, k][nan_indices] = np.nan  # set them to nan again  - this is useful when we plot
        for signed_edge in loop:
            ifg_index = edge_to_indexed_ifgs[signed_edge.edge].index
            #  the variable check_ps is increased by 1 for that pixel
            # make sure we are not incrementing the nan positions in the closure
            # as we don't know the PS of these pixels and also they were converted to zero before threshold check
            # Therefore, we leave them out of check_ps, i.e., we don't increment their check_ps values
            check_ps[np.logical_and(indices_breaching_threshold, ~nan_indices), ifg_index] += 1

    return closure, check_ps, num_occurences_each_ifg
