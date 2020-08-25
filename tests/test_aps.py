#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
"""
This Python module contains tests for the APS PyRate module.
"""
import os
import shutil
import pytest
import numpy as np
from pyrate import conv2tif, prepifg, correct
from pyrate.configuration import Configuration, MultiplePaths
import pyrate.core.config as cf
from pyrate.core.aps import wrap_spatio_temporal_filter, _interpolate_nans
from pyrate.core import shared
from tests import common


@pytest.fixture(params=["linear", "nearest", "cubic"])
def slpnanfill_method(request):
    return request.param


def test_interpolate_nans(slpnanfill_method):
    arr = np.random.rand(20, 10, 5)
    arr[arr < 0.1] = np.nan  # insert some nans
    assert np.sum(np.isnan(arr)) != 0  # some nans present
    _interpolate_nans(arr, method=slpnanfill_method)
    assert np.sum(np.isnan(arr)) == 0  # should not be any nans


def test_slpfilter():
    # TODO
    pass


def test_slp_filter():
    # TODO
    pass


def test_temporal_low_pass_filter():
    # TODO
    pass


def test_tlpfilter():
    # TODO
    pass


# APS correction using spatio-temporal filter
# apsest: ON = 1, OFF = 0
# Spatial low-pass filter parameters
# slpfmethod: filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 (greater than zero) in km for both butterworth and gaussian filters
# slpforder: order n for butterworth filter (default 1)
# slpnanfill: 1 for interpolation, 0 for zero fill
# slpnanfill_method: linear, nearest, cubic; only used when slpnanfill=1
# Temporal low-pass filter parameters
# tlpfmethod: 1 = Gaussian, 2 = Triangular, 3 = Mean filter
# tlpfcutoff: cutoff t0 for gaussian filter in year;
# tlpfpthr: valid pixel threshold;
# slpfmethod:     2
# slpfcutoff:     0.001
# slpforder:      1
# slpnanfill:     1
# slpnanfill_method:  cubic
# tlpfmethod:   3
# tlpfcutoff:   0.25
# tlpfpthr:     1


@pytest.fixture(params=[1, 2])
def slpfmethod(request):
    return request.param


@pytest.fixture(params=[0.001, 0.01])
def slpfcutoff(request):
    return request.param


@pytest.fixture(params=[1, 2])
def slpforder(request):
    return request.param


class TestAPSErrorCorrectionsOnDiscReused:

    @classmethod
    def setup_method(cls):
        cls.conf = common.TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._update_params_with_tiles(cls.params)
        correct._create_ifg_dict(cls.params)
        multi_paths = cls.params[cf.INTERFEROGRAM_FILES]
        cls.ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        cls.ifgs = [shared.Ifg(i) for i in cls.ifg_paths]
        for i in cls.ifgs:
            i.open()
        shared.save_numpy_phase(cls.ifg_paths, cls.params)
        correct.mst_calc_wrapper(cls.params)

    @classmethod
    def teardown_method(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    @pytest.mark.slow
    def test_aps_error_files_on_disc(self, slpfmethod, slpfcutoff, slpforder):
        self.params[cf.SLPF_METHOD] = slpfmethod
        self.params[cf.SLPF_CUTOFF] = slpfcutoff
        self.params[cf.SLPF_ORDER] = slpforder
        wrap_spatio_temporal_filter(self.params)

        # test_orb_errors_written
        aps_error_files = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        assert all(p.exists() for p in aps_error_files)
        last_mod_times = [os.stat(o).st_mtime for o in aps_error_files]
        phase_prev = [i.phase_data for i in self.ifgs]

        # run aps error removal again
        wrap_spatio_temporal_filter(self.params)
        aps_error_files2 = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        # if files are written again - times will change
        last_mod_times_2 = [os.stat(o).st_mtime for o in aps_error_files2]
        # test_aps_error_reused_if_params_unchanged
        assert all(a == b for a, b in zip(last_mod_times, last_mod_times_2))
        phase_now = [i.phase_data for i in self.ifgs]

        # run aps error correction once mroe
        wrap_spatio_temporal_filter(self.params)
        aps_error_files3 = [MultiplePaths.aps_error_path(i, self.params) for i in self.ifg_paths]
        last_mod_times_3 = [os.stat(o).st_mtime for o in aps_error_files3]
        assert all(a == b for a, b in zip(last_mod_times, last_mod_times_3))
        phase_again = [i.phase_data for i in self.ifgs]
        np.testing.assert_array_equal(phase_prev, phase_now)
        np.testing.assert_array_equal(phase_prev, phase_again)
