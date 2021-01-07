from datetime import date
from pathlib import Path
import random
import numpy as np
import pytest
from pyrate.constants import PYRATEPATH
from pyrate.core import config as cf
from pyrate.core.phase_closure.mst_closure import (
    find_closed_loops, Edge, SignedWeightedEdge, SignedEdge, setup_edges,
    add_signs_and_weights_to_loops, sort_loops_based_on_weights_and_date, WeightedLoop,
    find_signed_closed_loops
)
from pyrate.core.phase_closure.closure_check import (
    discard_loops_containing_max_ifg_count,
)

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs


def test_discard_loops_containing_max_ifg_count(geotiffs):
    loops1 = retain_loops(geotiffs)
    loops2 = retain_loops(geotiffs)
    m_weights = [m.weight for m in loops1]
    s_weights = [m.weight for m in loops2]
    np.testing.assert_array_equal(m_weights, s_weights)


def retain_loops(tifs):
    loops = find_signed_closed_loops(tifs)
    sorted_loops = sort_loops_based_on_weights_and_date(loops)
    params = {
        cf.MAX_LOOP_COUNT_FOR_EACH_IFGS: 2,
        cf.MAX_LOOP_LENGTH: 3
    }
    retained_loops_meeting_max_loop_criretia = [sl for sl in sorted_loops
                                                if len(sl) <= params[cf.MAX_LOOP_LENGTH]]
    msg = f"After applying MAX_LOOP_LENGTH={params[cf.MAX_LOOP_LENGTH]} criteria, " \
          f"{len(retained_loops_meeting_max_loop_criretia)} loops are retained"
    print(msg)
    retained_loops = discard_loops_containing_max_ifg_count(retained_loops_meeting_max_loop_criretia, params)
    msg = f"After applying MAX_LOOP_COUNT_FOR_EACH_IFGS={params[cf.MAX_LOOP_COUNT_FOR_EACH_IFGS]} criteria, " \
          f"{len(retained_loops)} loops are retained"
    print(msg)
    return retained_loops


# def test_drop_ifgs_if_not_part_of_any_loop(geotiffs):
#     loops1 = retain_loops(geotiffs)
#     selected_tifs1 = drop_ifgs_if_not_part_of_any_loop(geotiffs, loops1)
#
#     loops2 = retain_loops(geotiffs)
#     selected_tifs2 = drop_ifgs_if_not_part_of_any_loop(geotiffs, loops2)
#     assert all([a == b for a, b in zip(selected_tifs1, selected_tifs2)])
