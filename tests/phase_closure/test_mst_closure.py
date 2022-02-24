#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


from datetime import date
import numpy as np
import pytest
from pyrate.constants import PYRATEPATH
from pyrate.core.phase_closure.mst_closure import (
    __find_closed_loops, Edge, SignedWeightedEdge, SignedEdge, __setup_edges,
    __add_signs_and_weights_to_loops, sort_loops_based_on_weights_and_date, WeightedLoop,
    __find_signed_closed_loops
)
import pyrate.constants as C
from tests.phase_closure.common import IfgDummy

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture
def geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs


@pytest.fixture
def all_loops(geotiffs):
    edges = __setup_edges(geotiffs)
    loops = __find_closed_loops(edges, max_loop_length=100)
    assert len(loops) == 541
    return loops


@pytest.fixture
def edges(geotiffs):
    all_edges = __setup_edges(geotiffs)
    return all_edges


@pytest.fixture
def signed_loops(all_loops, edges):
    loops = __add_signs_and_weights_to_loops(all_loops, edges)
    return loops


@pytest.fixture(params=[True, False])
def weight(request):
    return request.param


def test_setup_edges(geotiffs):
    edges = __setup_edges(geotiffs)
    assert len(edges) == len(geotiffs) == 30
    assert isinstance(edges[0], Edge)


def test_associate_ifgs_with_loops(signed_loops, geotiffs):
    assert len(geotiffs) == 30
    assert len(signed_loops) == 541
    assert isinstance(signed_loops[0], WeightedLoop)
    swe = signed_loops[0].loop[0]
    assert isinstance(swe, SignedWeightedEdge)
    assert isinstance(swe.signed_edge, SignedEdge)
    assert isinstance(swe.edge, Edge)
    assert isinstance(swe.first, date)
    assert isinstance(swe.second, date)


def test_sort_loops_based_on_weights_and_date(geotiffs):
    ifg_files = [IfgDummy(ifg_path) for ifg_path in geotiffs]
    params = {
        C.INTERFEROGRAM_FILES: ifg_files,
        C.MAX_LOOP_LENGTH: 100
    }
    weighted_loops = sort_loops_based_on_weights_and_date(params)
    assert len(weighted_loops) == 541
    # order
    weights = [w.weight for w in weighted_loops]
    earliest_dates = [w.earliest_date for w in weighted_loops]
    assert np.all(np.diff(weights) >= 0)

    for i, (w, d, wl) in enumerate(zip(weights[:-1], earliest_dates[:-1], weighted_loops[:-1])):
        sub_list = [d]
        if w == weights[i+1]:
            sub_list.append(earliest_dates[i+1])
        if len(sub_list) > 1:
            tds = np.array([td.days for td in np.diff(sub_list)])
            assert np.all(tds >= 0)  # assert all dates are increasing for same weights


def test_add_signs_and_weights_to_loops(closure_params):
    geotiffs = closure_params["geotiffs"]

    """also tests find_signed_closed_loops"""
    all_edges = __setup_edges(geotiffs)
    all_loops = __find_closed_loops(all_edges, max_loop_length=closure_params[C.MAX_LOOP_LENGTH])
    loops1 = __add_signs_and_weights_to_loops(all_loops, all_edges)

    all_edges = __setup_edges(geotiffs)
    all_loops = __find_closed_loops(all_edges, closure_params[C.MAX_LOOP_LENGTH])
    loops2 = __add_signs_and_weights_to_loops(all_loops, all_edges)

    compare_loops(loops1, loops2)


def compare_loops(loops1, loops2):
    for i, (l1, l2) in enumerate(zip(loops1, loops2)):
        assert l1.weight == l2.weight
        for ll1, ll2 in zip(l1.loop, l2.loop):
            assert ll1.weight == ll2.weight
            assert ll1.first == ll2.first
            assert ll1.second == ll2.second
            assert ll1.sign == ll2.sign
    assert i == len(loops1) - 1
    assert i == 540


def test_find_signed_closed_loops(closure_params):
    loops1 = __find_signed_closed_loops(closure_params)
    loops2 = __find_signed_closed_loops(closure_params)
    compare_loops(loops1, loops2)


def test_sort_loops_based_on_weights_and_date_2(closure_params):
    sorted_loops1 = sort_loops_based_on_weights_and_date(closure_params)
    sorted_loops2 = sort_loops_based_on_weights_and_date(closure_params)
    compare_loops(sorted_loops1, sorted_loops2)
