from typing import List
from datetime import date
import pytest
import numpy as np
import networkx as nx
from pyrate.core.phase_closure.mst_closure import Edge, __setup_edges, __find_closed_loops
from pyrate.core.phase_closure.collect_loops import dedupe_loops, find_loops


def test_collect_loops():
    graph = np.array(
        [
            [0, 1, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [0, 1, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [0, 1, 0, 1, 0]
        ]
    )

    n = 4
    count, all_loops = find_loops(graph, n)
    assert count == 6
    deduped_loops = dedupe_loops(all_loops)
    np.testing.assert_array_equal(deduped_loops, [[0, 1, 2, 3], [0, 1, 4, 3], [1, 2, 3, 4]])


def test_count_loops():
    graph = np.array(
            [
                [0, 1, 1, 1],
                [1, 0, 1, 1],
                [1, 1, 0, 1],
                [1, 1, 1, 0],
            ]
    )

    n = 4
    count, all_loops = find_loops(graph, n)
    assert len(all_loops) == 6
    np.testing.assert_array_equal(all_loops, [[0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2],
                                              [0, 3, 2, 1]])
    deduped_loops = dedupe_loops(all_loops)
    np.testing.assert_array_equal(deduped_loops, [all_loops[0]])


def __find_closed_loops_nx(edges: List[Edge], max_loop_length: int) -> List[List[date]]:
    g = nx.Graph()
    edges = [(we.first, we.second) for we in edges]
    g.add_edges_from(edges)
    dg = nx.DiGraph(g)
    print(f"Evaluating all possible loops using NetworkX simple_cycles function")
    simple_cycles = list(nx.simple_cycles(dg))
    print(f"Total number of possible loops is {len(simple_cycles)}")
    print(f"Discarding loops with less than 3 edges and more than {max_loop_length} edges")
    loop_subset = [scc for scc in simple_cycles
                   if
                   (len(scc) > 2)  # three or more edges reqd for closed loop
                   and
                   (len(scc) <= max_loop_length)  # discard loops exceeding max loop length
                   ]
    print(f"Number of remaining loops is {len(loop_subset)}")

    # also discard loops when the loop members are the same
    return dedupe_loops(loop_subset)


@pytest.fixture(scope='module')
def available_edges(cropa_geotifs):
    available_edges = __setup_edges(cropa_geotifs)
    return available_edges


def max_loop_length(available_edges):
    all_possible_loops = __find_closed_loops_nx(available_edges, max_loop_length=100)
    max_length = max([len(l) for l in all_possible_loops])
    return max_length


def test_find_closed_loops_vs_networkx(available_edges):
    max_length = max_loop_length(available_edges)

    for n in range(max_length + 1):
        print(f'====checking for max_loop_length {n}=====')
        networkx_loops = __find_closed_loops_nx(available_edges, max_loop_length=n)
        networkx_loops = [sorted(l) for l in networkx_loops]
        networkx_set = {tuple(l) for l in networkx_loops}
        our_loops = __find_closed_loops(available_edges, max_loop_length=n)
        our_loops = [sorted(l) for l in our_loops]
        our_set = {tuple(l) for l in our_loops}
        assert networkx_set == our_set
