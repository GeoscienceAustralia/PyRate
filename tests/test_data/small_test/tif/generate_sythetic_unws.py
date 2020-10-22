from collections import namedtuple
import itertools
from typing import List, Dict
from datetime import datetime, date
import numpy as np
import networkx as nx
from tests import common

np.random.seed(10)
num = 200

LoopIfgMock = namedtuple('LoopIfg', ['first', 'second', 'nan_fraction'])


def generate_random_edges(ifgs):
    day = np.random.choice(range(1, 29), num)
    month = np.random.choice(range(1, 13), num)
    year = np.random.choice(range(2005, 2009), num)

    primary = [date(y, m, d) for y, m, d in zip(year, month, day)]

    day = np.random.choice(range(1, 29), num)
    month = np.random.choice(range(1, 13), num)
    year = np.random.choice(range(2005, 2009), num)

    secondary = [date(y, m, d) for y, m, d in zip(year, month, day)]

    # TODO: support ifg weights
    weights = np.ones_like(day)
    edges_with_weights = [(i.first, i.second, 1) for i in ifgs]
    original_loopifgs = [LoopIfgMock(i.first, i.second, 1) for i in ifgs]
    all_edges = edges_with_weights + [(p, s, w) for p, s, w in zip(primary, secondary, weights)]
    all_ifgs = original_loopifgs + [LoopIfgMock(p, s, w) for p, s, w in zip(primary, secondary, weights)]
    return all_edges, all_ifgs


def find_closed_loops(all_edges):
    g = nx.Graph()
    g.add_weighted_edges_from(all_edges)

    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]  # discard edges
    # TODO: also discard loops when the loop members are the same?
    return simple_cycles


def associate_ifgs_with_loops(loops):
    loops_to_ifgs = {}
    for i, l in enumerate(loops):
        loop_ifgs = []
        for ii, _ in enumerate(l[:-1]):
            ifg = LoopIfgMock(l[ii], l[ii+1], 1)
            loop_ifgs.append(ifg)

        ifg = LoopIfgMock(l[ii+1], l[0], 1)
        loop_ifgs.append(ifg)  # add the last to the first node, i.e., close the loop
        loops_to_ifgs[i] = loop_ifgs

    return loops_to_ifgs


def filter_to_ifgs_used_in_loops(loops_to_ifgs: Dict[(int, LoopIfgMock)]):
    return set(itertools.chain(*loops_to_ifgs.values()))


def create_mock_ifgs_data_unless_already_available():
    pass

_, ifgs = common.copy_and_setup_small_data()
all_edges, all_mock_ifgs = generate_random_edges(ifgs)
all_loops = find_closed_loops(all_edges)
loops_to_ifgs = associate_ifgs_with_loops(all_loops)
ifgs_used = filter_to_ifgs_used_in_loops(loops_to_ifgs)

import IPython; IPython.embed()
