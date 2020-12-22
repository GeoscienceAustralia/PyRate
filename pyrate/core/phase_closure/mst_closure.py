from collections import namedtuple
from typing import List, Union
from datetime import date
import networkx as nx
from pyrate.core.shared import dem_or_ifg

Edge = namedtuple('Edge', ['first', 'second'])
SignedEdge = namedtuple('SignedEdge', ['edge', 'sign'])
SignedWeightedEdge = namedtuple('SignedWeightedEdge', ['SignedEdge', 'weight'])


class WeightedLoop:

    def __init__(self, loop: List[SignedWeightedEdge]):
        self.loop = loop

    @property
    def weight(self):
        return sum({l.weight for l in self.loop})

    @property
    def earliest_date(self):
        return min({loop.SignedEdge.edge.first for loop in self.loop})


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


def find_closed_loops(edges: List[Edge]) -> List[List[date]]:
    g = nx.Graph()
    edges = [(we.first, we.second) for we in edges]
    g.add_edges_from(edges)
    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]  # discard edges

    # also discard loops when the loop members are the same
    return discard_edges_with_same_members(simple_cycles)


def add_signs_and_weights_to_loops(loops, available_edges) -> List[List[SignedEdge]]:
    weighted_signed_loops = []
    available_edges = set(available_edges)  # hash it once for O(1) lookup
    for i, l in enumerate(loops):
        weighted_signed_loop = []
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
            weighted_signed_edge = SignedWeightedEdge(
                signed_edge,
                (signed_edge.edge.second - signed_edge.edge.first).days
            )
            weighted_signed_loop.append(weighted_signed_edge)

        weighted_signed_loops.append(weighted_signed_loop)

    return weighted_signed_loops


def setup_edges(ifg_files: List['str']) -> List[Edge]:
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    return [Edge(i.first, i.second) for i in ifgs]


def find_signed_closed_loops(ifg_files: List[str]) -> List[List[SignedEdge]]:
    available_edges = setup_edges(ifg_files)
    all_loops = find_closed_loops(available_edges)  # find loops with weights
    signed_loops = add_signs_and_weights_to_loops(all_loops, available_edges)
    signed_loops.sort(key=lambda x: x.weight)
    return signed_loops


def sort_loops_based_on_weights_and_date(signed_loops: List[List[SignedEdge]]) -> List[WeightedLoop]:
    weighted_loops = [WeightedLoop(sl) for sl in signed_loops]
    weighted_loops.sort(key=lambda x: [x.weight, x.earliest_date])  # sort based on weights and date
    return weighted_loops
