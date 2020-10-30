from collections import namedtuple
from pathlib import Path
from typing import List, Dict
import numpy as np
from pyrate.core.shared import Ifg, dem_or_ifg
from pyrate.core.phase_closure.mst_closure import Edge, WeightedEdge, SignedEdge, mst_closure_wrapped

THRESHOLD_TO_REMOVE_PIXEL = 0.25
LARGE_DEVIATION_THRESHOLD_FOR_PIXEL = 3.14152  # pi
THRESHOLD_TO_REMOVE_IFG = 0.5
MAX_LOOP_LENGTH = 5  # loops upto this many edges are considered for closure checks
IndexedIfg = namedtuple('IndexedIfg', ['index', 'Ifg'])


def create_ifg_edge_dict(ifg_files) -> Dict[Edge, IndexedIfg]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    return {Edge(ifg.first, ifg.second): IndexedIfg(index, ifg) for index, ifg in enumerate(ifgs)}


def sum_phase_values_for_each_loop(ifg_files: List[str], loops: List[List[SignedEdge]], threshold: float = None):
    edge_to_indexed_ifgs = create_ifg_edge_dict(ifg_files)
    ifgs = [v.Ifg for v in edge_to_indexed_ifgs.values()]
    n_ifgs = len(ifgs)
    n_loops = len(loops)
    closure = np.zeros(shape=(ifgs[0].phase_data.shape + (n_loops,)))

    # number of occurences per IFG
    ifg_num = np.zeros(shape=n_ifgs)

    # initiate variable for check of unwrapping issues at the same pixels in all loops
    check_ps = np.zeros(shape=(ifgs[0].phase_data.shape + (n_ifgs,)))
    for k, loop in enumerate(loops):
        for signed_edge in loop:
            ifg = edge_to_indexed_ifgs[signed_edge.edge].Ifg
            index = edge_to_indexed_ifgs[signed_edge.edge].index  # iline
            closure[:, :, k] += signed_edge.sign * ifg.phase_data
            ifg_num[index] += 1

        closure[:, :, k] -= np.nanmedian(closure[:, :, k])

        if threshold:
            indices_breaching_threshold = np.absolute(closure[:, :, k] > threshold)
            for signed_edge in loop:
                index = edge_to_indexed_ifgs[signed_edge.edge].index  # iline
                #  the variable check_ps is increased by 1 for that pixel
                check_ps[indices_breaching_threshold, index] += 1

    return closure, check_ps


new_test_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*_unw.tif')
ifg_files = [f.as_posix() for f in new_test_files]
signed_loops = mst_closure_wrapped(ifg_files=ifg_files)
retained_loops = [sl for sl in signed_loops if len(sl) <= MAX_LOOP_LENGTH]
closure, check_ps = sum_phase_values_for_each_loop(ifg_files, retained_loops, LARGE_DEVIATION_THRESHOLD_FOR_PIXEL)
import IPython; IPython.embed()

