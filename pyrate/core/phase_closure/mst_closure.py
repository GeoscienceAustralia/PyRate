from collections import namedtuple
from pathlib import Path
from typing import List, Dict, Union, Tuple
from datetime import datetime, date, timedelta
import numpy as np
import networkx as nx
from pyrate.core.shared import dem_or_ifg

Edge = namedtuple('Edge', ['first', 'second'])
SignedEdge = namedtuple('SignedEdge', ['edge', 'sign'])
WeightedEdge = namedtuple('WeightedEdge', ['first', 'second', 'weight'])


def discard_edges_with_same_members(simple_cycles):
    seen_sc_sets = set()
    filtered_sc = []
    for sc in simple_cycles:
        loop = sc[:]
        sc.sort()
        sc = tuple(sc)
        if sc not in seen_sc_sets:
            filtered_sc.append(loop)
        seen_sc_sets.add(sc)
    return filtered_sc


def find_closed_loops(weighted_edges: List[WeightedEdge]) -> List[List[date]]:
    g = nx.Graph()
    g.add_weighted_edges_from(weighted_edges)
    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]  # discard edges

    # also discard loops when the loop members are the same
    return discard_edges_with_same_members(simple_cycles)


def add_signs_to_loops(loops, available_edges) -> Dict[int, List[SignedEdge]]:
    signed_loops = dict()
    available_edges = set(available_edges)  # hash it once for O(1) lookup
    for i, l in enumerate(loops):
        signed_loop = []
        l.append(l[0])  # add the closure loop
        for ii, ll in enumerate(l[:-1]):
            if l[ii+1] > ll:
                edge = Edge(ll, l[ii+1])
                assert edge in available_edges
                signed_edge = SignedEdge(edge, 1)  # opposite direction of ifg
            else:
                edge = Edge(l[ii+1], ll)
                assert edge in available_edges
                signed_edge = SignedEdge(edge, -1)  # in direction of ifg
            signed_loop.append(signed_edge)

        signed_loops[i] = signed_loop

    return signed_loops


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


def mst_closure_wrapped(ifg_files: List[str]):
    available_edges = setup_edges(ifg_files)
    weighted_edges = setup_edges(ifg_files, weighted=True)

    all_loops = find_closed_loops(weighted_edges)  # find loops with weights
    signed_loops = add_signs_to_loops(all_loops, available_edges)
    return signed_loops


# new_test_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*_unw.tif')
# ifg_files = [f.as_posix() for f in new_test_files]
# signed_loops = mst_closure_wrapped(ifg_files=ifg_files)
# import IPython; IPython.embed()