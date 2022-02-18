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


import numpy as np

import pyrate.constants as C
from pyrate.core.phase_closure.mst_closure import (
    sort_loops_based_on_weights_and_date,
)
from pyrate.core.phase_closure.closure_check import (
    discard_loops_containing_max_ifg_count,
    __drop_ifgs_if_not_part_of_any_loop
)


def test_discard_loops_containing_max_ifg_count(closure_params):
    params = closure_params
    loops1 = retain_loops(params)
    loops2 = retain_loops(params)
    m_weights = [m.weight for m in loops1]
    s_weights = [m.weight for m in loops2]
    np.testing.assert_array_equal(m_weights, s_weights)


def retain_loops(params):
    sorted_loops = sort_loops_based_on_weights_and_date(params)
    params[C.MAX_LOOP_REDUNDANCY] = 2
    params[C.MAX_LOOP_LENGTH] = 3
    retained_loops_meeting_max_loop_criteria = [sl for sl in sorted_loops
                                                if len(sl) <= params[C.MAX_LOOP_LENGTH]]
    msg = f"After applying MAX_LOOP_LENGTH={params[C.MAX_LOOP_LENGTH]} criteria, " \
          f"{len(retained_loops_meeting_max_loop_criteria)} loops are retained"
    print(msg)
    retained_loops = discard_loops_containing_max_ifg_count(retained_loops_meeting_max_loop_criteria, params)
    msg = f"After applying MAX_LOOP_REDUNDANCY={params[C.MAX_LOOP_REDUNDANCY]} criteria, " \
          f"{len(retained_loops)} loops are retained"
    print(msg)
    return retained_loops


def test_drop_ifgs_if_not_part_of_any_loop(closure_params):
    params = closure_params
    params[C.NO_DATA_VALUE] = 0.0
    geotiffs = params['geotiffs']

    loops1 = retain_loops(params)
    selected_tifs1 = __drop_ifgs_if_not_part_of_any_loop(geotiffs, loops1, params)

    loops2 = retain_loops(params)
    selected_tifs2 = __drop_ifgs_if_not_part_of_any_loop(geotiffs, loops2, params)
    assert all([a == b for a, b in zip(selected_tifs1, selected_tifs2)])
