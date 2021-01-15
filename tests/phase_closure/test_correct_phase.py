#   This Python module is part of the PyRate software package.
#
#   Copyright 2021 Geoscience Australia
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
import pytest
from pyrate.core.phase_closure.correct_phase import find_breached_pixels, __recover_pixels


@pytest.fixture(params=list(range(1, 50, 3)))
def multiple(request):
    return request.param


@pytest.fixture
def abs_closure(multiple):
    return multiple * 2 * np.pi + np.random.rand(20, 30) - 0.5  # random numbers 2*pi +- 0.5


def test__recover_pixels(abs_closure):
    large_dev_thr = 0.4
    max_closure = np.max(abs_closure)
    recovered_pixels = __recover_pixels(abs_closure, large_dev_thr, max_closure)
    # a thresh of 0.4 should recover at around 80% of these numbers
    assert np.sum(recovered_pixels) > abs_closure.size * 0.75


def test_find_breached_pixels():
    pass
