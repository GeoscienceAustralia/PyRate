from datetime import datetime, date
import numpy as np
import networkx as nx
from tests import common

np.random.seed(10)
num = 200


def generate_random_edges(ifgs):
    day = np.random.choice(range(1, 29), num)
    month = np.random.choice(range(1, 13), num)
    year = np.random.choice(range(2005, 2009), num)

    primary = [date(y, m, d) for y, m, d in zip(year, month, day)]

    day = np.random.choice(range(1, 29), num)
    month = np.random.choice(range(1, 13), num)
    year = np.random.choice(range(2005, 2009), num)

    secondary = [date(y, m, d) for y, m, d in zip(year, month, day)]

    weights = np.random.rand(num) * 0.1
    edges_with_weights = [(i.first, i.second, i.nan_fraction) for i in ifgs]
    all_edges = edges_with_weights + [(p, s, w) for p, s, w in zip(primary, secondary, weights)]
    return all_edges


def find_closed_loops(all_edges):
    g = nx.Graph()
    g.add_weighted_edges_from(all_edges)

    dg = nx.DiGraph(g)
    simple_cycles = nx.simple_cycles(dg)  # will have all edges
    simple_cycles = [scc for scc in simple_cycles if len(scc) > 2]
    return simple_cycles

_, ifgs = common.copy_and_setup_small_data()
all_edges = generate_random_edges(ifgs)
all_loops = find_closed_loops(all_edges)

