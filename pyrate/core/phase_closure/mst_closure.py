from collections import namedtuple
from pathlib import Path
from typing import List, Dict, Union
from datetime import datetime, date, timedelta
import numpy as np
import networkx as nx
from pyrate.core.shared import dem_or_ifg

Edge = namedtuple('Edge', ['first', 'second'])
WeightedEdge = namedtuple('Edge', ['first', 'second', 'weight'])


def find_closed_loops(weighted_edges: List[WeightedEdge]) -> List[date]:
    g = nx.Graph()
    g.add_weighted_edges_from(weighted_edges)
    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]  # discard edges
    # TODO: also discard loops when the loop members are the same?
    return simple_cycles


def filter_loops_to_increasing_sequence_loops(loops):
    filtered_loops = set()
    for i, l in enumerate(loops):
        gt = [l[ii+1] > l[ii] for ii in range(len(l)-1)]
        if np.sum(~np.array(gt)) > 1:  # more than one sequence
            continue  # don't take this one
        elif np.sum(~np.array(gt)) == 1:  # can be turned into increasing sequence
            # take this but change to increasing sequence
            index_of_false = gt.index(False)
            new_l = l[index_of_false+1:] + l[0:index_of_false + 1]
            if all([new_l[ii+1] > new_l[ii] for ii in range(len(new_l)-1)]):
                filtered_loops.add(tuple(new_l))
        else:
            if all(gt):  # increasing sequence returned by networkx
                filtered_loops.add(tuple(l))
    return filtered_loops


def associate_ifgs_with_loops(loops, available_edges):
    loops_to_ifgs = {}

    for i, l in enumerate(loops):
        loop_ifgs = []
        for ii, ll in enumerate(l[:-1]):
            ifg = Edge(ll, l[ii + 1], )
            loop_ifgs.append(ifg)
        ifg = Edge(l[0], l[-1])
        assert ifg in available_edges, f'{ifg} not found in original list'
        loop_ifgs.append(ifg)
        loops_to_ifgs[i] = loop_ifgs

    return loops_to_ifgs


def setup_edges(ifg_files: List['str'], weighted: bool = False) -> List[Union[Edge, WeightedEdge]]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    if weighted:
        return [WeightedEdge(i.first, i.second, i.nan_fraction) for i in ifgs]
    else:
        return [Edge(i.first, i.second) for i in ifgs]


def mst_closure_wrapped(ifg_files):
    edges = setup_edges(ifg_files)
    weighted_edges = setup_edges(ifg_files, weighted=True)

    all_loops = find_closed_loops(weighted_edges)
    increasing_loops = filter_loops_to_increasing_sequence_loops(all_loops)
    loops_to_ifgs = associate_ifgs_with_loops(increasing_loops, edges)
    return loops_to_ifgs


# new_test_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*_unw.tif')
# ifg_files = [f.as_posix() for f in new_test_files if '_dem.tif' not in f.as_posix()]
# loops_to_ifgs = mst_closure_wrapped(ifg_files=ifg_files)
