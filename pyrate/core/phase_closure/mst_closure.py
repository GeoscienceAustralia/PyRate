from collections import namedtuple
import itertools
from typing import List, Dict
from datetime import datetime, date, timedelta
import numpy as np
import networkx as nx
from tests import common


def find_closed_loops(ifgs):
    edges_with_weights = [(i.first, i.second, i.nan_fraction) for i in ifgs]
    g = nx.Graph()
    g.add_weighted_edges_from(edges_with_weights)

    return