from collections import namedtuple
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np
from pyrate.core.shared import Ifg, dem_or_ifg
from pyrate.core.phase_closure.mst_closure import Edge, WeightedEdge, SignedEdge, find_signed_closed_loops
from pyrate.core.logger import pyratelogger as log
IndexedIfg = namedtuple('IndexedIfg', ['index', 'Ifg'])


def create_ifg_edge_dict(ifg_files) -> Dict[Edge, IndexedIfg]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    return {Edge(ifg.first, ifg.second): IndexedIfg(index, ifg) for index, ifg in enumerate(ifgs)}


def sum_phase_values_for_each_loop(ifg_files: List[str], loops: List[List[SignedEdge]], threshold: float):
    edge_to_indexed_ifgs = create_ifg_edge_dict(ifg_files)
    ifgs = [v.Ifg for v in edge_to_indexed_ifgs.values()]
    n_ifgs = len(ifgs)
    n_loops = len(loops)
    closure = np.zeros(shape=(ifgs[0].phase_data.shape + (n_loops,)))

    num_occurences_each_ifg = np.zeros(shape=n_ifgs)

    # initiate variable for check of unwrapping issues at the same pixels in all loops
    check_ps = np.zeros(shape=(ifgs[0].phase_data.shape + (n_ifgs,)))
    for k, loop in enumerate(loops):
        for signed_edge in loop:
            ifg = edge_to_indexed_ifgs[signed_edge.edge].Ifg
            ifg_index = edge_to_indexed_ifgs[signed_edge.edge].index
            closure[:, :, k] += signed_edge.sign * ifg.phase_data
            num_occurences_each_ifg[ifg_index] += 1

        closure[:, :, k] -= np.median(closure[:, :, k])

        indices_breaching_threshold = np.absolute(closure[:, :, k]) > threshold
        for signed_edge in loop:
            ifg_index = edge_to_indexed_ifgs[signed_edge.edge].index
            #  the variable check_ps is increased by 1 for that pixel
            check_ps[indices_breaching_threshold, ifg_index] += 1

    return closure, check_ps, num_occurences_each_ifg
