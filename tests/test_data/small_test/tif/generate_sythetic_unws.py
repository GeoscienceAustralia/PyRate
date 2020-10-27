from collections import namedtuple
import itertools
from typing import List, Dict
from datetime import datetime, date, timedelta
import numpy as np
import networkx as nx
from pyrate.core.shared import dem_or_ifg

Edge = namedtuple('Edge', ['first', 'second'])


def add_more_ifgs(ifgs):
    edges_with_weights = [(i.first, i.second) for i in ifgs]

    edge1 = (date(2006, 8, 28), date(2006, 11, 6))
    edge2 = (date(2006, 11, 6), date(2007, 7, 9))
    edge3 = (date(2007, 6, 4), date(2007, 8, 13))
    edge4 = (date(2006, 11, 6), date(2007, 6, 4))
    edge5 = (date(2006, 8, 28), date(2007, 3, 26))
    edge6 = (date(2006, 12, 11), date(2007, 6, 4))
    edge7 = (date(2007, 2, 19), date(2007, 8, 13))
    edge8 = (date(2006, 10, 2), date(2007, 6, 4))
    edge9 = (date(2006, 10, 2), date(2006, 12, 2))
    edge10 = (date(2006, 12, 2), date(2007, 2, 19))
    return edges_with_weights + [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8, edge9, edge10]


def find_closed_loops(all_edges):
    g = nx.Graph()
    g.add_edges_from(all_edges)
    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]  # discard edges
    print('len(simple_cycles): ', len(simple_cycles))
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
    print('len(filtered_loops): ', len(filtered_loops))
    return filtered_loops


def associate_ifgs_with_loops(loops, available_edges):
    loops_to_ifgs = {}

    for i, l in enumerate(loops):
        loop_ifgs = []
        for ii, ll in enumerate(l[:-1]):
            ifg = Edge(ll, l[ii + 1])
            loop_ifgs.append(ifg)
        ifg = Edge(l[0], l[-1])
        assert ifg in available_edges, f'{ifg} not found in original list'
        loop_ifgs.append(ifg)
        loops_to_ifgs[i] = loop_ifgs

    return loops_to_ifgs


def filter_to_ifgs_used_in_loops(loops_to_ifgs: Dict[(int, Edge)]):
    return set(itertools.chain(*loops_to_ifgs.values()))


def create_mock_ifgs_unless_already_available(ifgs_used, original_ifgs, all_ifgs):

    pass


def setup_data(ifg_files):
    ifg_files.sort()
    ifgs = [dem_or_ifg(i) for i in ifg_files]
    for i in ifgs:
        i.open()
        i.nodata_value = 0
    return ifgs

# g = nx.Graph()
# g.add_edges_from(all_edges)
# nx.draw_networkx(g)
import matplotlib.pyplot as plt
# plt.show()
# import IPython; IPython.embed()
from pathlib import Path
new_test_files = Path('/home/sudipta/Documents/GEOTIFF').glob('*.tif')
ifgs = setup_data([f.as_posix() for f in new_test_files if '_dem.tif' not in f.as_posix()])
edges = [Edge(i.first, i.second) for i in ifgs]
all_loops = find_closed_loops(edges)
increasing_loops = filter_loops_to_increasing_sequence_loops(all_loops)
print(len(increasing_loops))
loops_to_ifgs = associate_ifgs_with_loops(increasing_loops, edges)

