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
"""
This Python module contains tests for the orbital.py PyRate module.
"""
import os
import shutil
import tempfile
import pytest
from itertools import product
from numpy import empty, dot, concatenate, float32
from numpy import nan, isnan, array
from os.path import join
from pathlib import Path
from datetime import date

import numpy as np
from numpy.linalg import pinv, inv
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_allclose
from scipy.linalg import lstsq
from osgeo import gdal, gdalconst

import pyrate.constants as C
import pyrate.core.orbital
from tests.common import small5_mock_ifgs, MockIfg
from pyrate.core.algorithm import first_second_ids, get_all_epochs
from pyrate.core.orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, \
    QUADRATIC, PART_CUBIC
from pyrate.core.orbital import OrbitalError, __orb_correction, __orb_inversion
from pyrate.core.orbital import get_design_matrix, get_network_design_matrix, orb_fit_calc_wrapper
from pyrate.core.orbital import _get_num_params, remove_orbital_error, network_orbital_correction
from pyrate.core.orbital import calc_network_orb_correction
from pyrate.core.shared import Ifg, mkdir_p
from pyrate.core.shared import nanmedian
from pyrate.core import roipac
from pyrate import correct, conv2tif, prepifg
from pyrate.configuration import Configuration, MultiplePaths
from pyrate.constants import ORB_ERROR_DIR
from tests import common
from tests.common import IFMS16, TEST_CONF_GAMMA
from tests.common import SML_TEST_LEGACY_ORBITAL_DIR
from tests.common import SML_TEST_TIF, PY37GDAL302
from tests.common import small_ifg_file_list

# TODO: Purpose of this variable? Degrees are 1, 2 and 3 respectively
DEG_LOOKUP = {
    2: PLANAR,
    5: QUADRATIC,
    6: PART_CUBIC}

NUM_COEF_LOOKUP = {
    PLANAR: 2,
    QUADRATIC: 5,
    PART_CUBIC: 6}


class TestSingleDesignMatrixTests:
    """
    Tests to verify correctness of basic planar & quadratic design matrices or
    DMs. This class serves two purposes, ensuring the independent method DMs are
    produced correctly. Secondly, these indivdual DMs are subsets of the larger
    DM 'grid' required for the networked orbital correction method.
    """

    def setup_class(cls):
        # faked cell sizes
        cls.xs = 0.75
        cls.ys = 0.8
        cls.ifg = Ifg(join(SML_TEST_TIF, 'geo_060619-061002_unw.tif'))
        cls.ifg.open()
        cls.ifg.nodata_value = 0

        cls.m = MockIfg(cls.ifg, 3, 4)
        cls.m.x_size = cls.xs
        cls.m.y_size = cls.ys

    # tests for planar model

    def test_create_planar_dm(self):
        act = get_design_matrix(self.m, PLANAR, intercept=False, scale=100)
        assert act.shape == (self.m.num_cells, 2)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset=False)
        assert_array_equal(act, exp)

    def test_create_planar_dm_offsets(self):
        act = get_design_matrix(self.m, PLANAR, intercept=True, scale=100)
        assert act.shape == (self.m.num_cells, 3)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset=True)
        assert_array_almost_equal(act, exp)

    # tests for quadratic model

    def test_create_quadratic_dm(self):
        act = get_design_matrix(self.m, QUADRATIC, intercept=False, scale=100)
        assert act.shape == (self.m.num_cells, 5)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset=False)
        assert_array_equal(act, exp)

    def test_create_quadratic_dm_offsets(self):
        act = get_design_matrix(self.m, QUADRATIC, intercept=True, scale=100)
        assert act.shape == (self.m.num_cells, 6)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset=True)
        assert_array_equal(act, exp)

    # tests for partial cubic model

    def test_create_partcubic_dm(self):
        act = get_design_matrix(self.m, PART_CUBIC, intercept=False, scale=100)
        assert act.shape == (self.m.num_cells, 6)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset=False)
        assert_array_equal(act, exp)

    def test_create_partcubic_dm_offsets(self):
        act = get_design_matrix(self.m, PART_CUBIC, intercept=True, scale=100)
        assert act.shape == (self.m.num_cells, 7)
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset=True)
        assert_array_equal(act, exp)

    # tests for unittest_dm() assuming network method

    def test_create_planar_dm_network(self):
        # networked method planar version should not have offsets col
        ncol_exp = 2
        exp = unittest_dm(self.m, NETWORK_METHOD, PLANAR, False)
        assert exp.shape == (self.m.num_cells, ncol_exp)
        exp2 = unittest_dm(self.m, NETWORK_METHOD, PLANAR, True)
        assert exp2.shape == (self.m.num_cells, ncol_exp)
        assert_array_equal(exp, exp2)

    def test_create_quadratic_dm_network(self):
        # quadratic version with networked method does not have offsets col
        ncol_exp = 5
        exp = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, False)
        assert exp.shape == (self.m.num_cells, ncol_exp)
        exp2 = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, True)
        assert exp2.shape == (self.m.num_cells, ncol_exp)
        assert_array_equal(exp, exp2)

    def test_create_partcubic_dm_network(self):
        # partial cubic version with networked method does not have offsets col
        ncol_exp = 6
        exp = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, False)
        assert exp.shape == (self.m.num_cells, ncol_exp)
        exp2 = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, True)
        assert exp2.shape == (self.m.num_cells, ncol_exp)
        assert_array_equal(exp, exp2)


class TestIndependentCorrection:
    """Test cases for the orbital correction component of PyRate."""

    @classmethod
    def setup_class(cls):
        cls.ifgs = small5_mock_ifgs()
        _add_nodata(cls.ifgs)

        for ifg in cls.ifgs:
            ifg.x_size = 90.0
            ifg.y_size = 89.5
            ifg.open()

    def alt_orbital_correction(self, ifg, deg, offset, scale):
        data = ifg.phase_data.reshape(ifg.num_cells)
        dm = get_design_matrix(ifg, deg, intercept=True, scale=scale)[~isnan(data)]
        fd = data[~isnan(data)].reshape((dm.shape[0], 1))

        dmt = dm.T
        invNbb = inv(dmt.dot(dm))
        orbparams = invNbb.dot(dmt.dot(fd))
        alt_params = lstsq(dm, fd)[0]
        # FIXME: precision
        assert_array_almost_equal(orbparams, alt_params, decimal=1)

        dm2 = get_design_matrix(ifg, deg, intercept=True, scale=scale)
        fullorb = np.reshape(np.dot(dm2, orbparams), ifg.phase_data.shape)
        if offset:
            offset_removal = nanmedian(np.ravel(ifg.phase_data - fullorb))
        else:
            offset_removal = 0

        fwd_correction = fullorb - offset_removal
        # ifg.phase_data -= (fullorb - offset_removal)
        return ifg.phase_data - fwd_correction

    def check_correction(self, degree, method, offset, decimal=2):
        orig = array([c.phase_data.copy() for c in self.ifgs])
        exp = [self.alt_orbital_correction(i, degree, offset, scale=100) for i in self.ifgs]
        params = dict()
        params[C.ORBITAL_FIT_METHOD] = method
        params[C.ORBITAL_FIT_DEGREE] = degree
        params[C.ORBFIT_OFFSET] = offset
        params[C.ORBFIT_INTERCEPT] = 1
        params[C.ORBFIT_SCALE] = 100
        params[C.PARALLEL] = False
        params[C.NO_DATA_VALUE] = 0
        params[C.NAN_CONVERSION] = False
        params[C.OUT_DIR] = tempfile.mkdtemp()
        params[C.ORBITAL_FIT_LOOKS_X] = 1
        params[C.ORBITAL_FIT_LOOKS_Y] = 1
        params[C.TEMP_MLOOKED_DIR] = tempfile.mkdtemp()
        for i in self.ifgs:
            i.mm_converted = True
        remove_orbital_error(self.ifgs, params)
        corrected = array([c.phase_data for c in self.ifgs])

        assert ~(orig == corrected).all()
        self.check_results(self.ifgs, orig)  # test shape, data is non zero

        # FIXME: is decimal=2 close enough?
        for i, (e, a) in enumerate(zip(exp, corrected)):
            assert_array_almost_equal(e, a, decimal=decimal)

    def check_results(self, ifgs, corrections):
        """Helper method for result verification"""
        for i, c in zip(ifgs, corrections):
            ys, xs = c.shape
            assert i.nrows == ys
            assert i.ncols == xs

            # ensure there is real data
            assert ~ isnan(i.phase_data).all()
            assert ~ isnan(c).all()
            assert c.ptp() != 0  # ensure range of values in grid

    def test_independent_correction_planar(self):
        self.check_correction(PLANAR, INDEPENDENT_METHOD, False)

    def test_independent_correction_planar_offsets(self):
        self.check_correction(PLANAR, INDEPENDENT_METHOD, True)

    def test_independent_correction_quadratic(self):
        self.check_correction(QUADRATIC, INDEPENDENT_METHOD, False)

    def test_independent_correction_quadratic_offsets(self):
        self.check_correction(QUADRATIC, INDEPENDENT_METHOD, True)

    def test_independent_correction_partcubic(self):
        self.check_correction(PART_CUBIC, INDEPENDENT_METHOD, False)

    def test_independent_correction_partcubic_offsets(self):
        self.check_correction(PART_CUBIC, INDEPENDENT_METHOD, True, decimal=1)


class TestError:
    """Tests for the networked correction method"""

    @classmethod
    def setup_method(cls):
        out_dir = tempfile.mkdtemp()
        cls.params = common.min_params(out_dir)
        cls.ifgs = small5_mock_ifgs()

    def test_invalid_ifgs_arg(self):
        # min requirement is 1 ifg, can still subtract one epoch from the other
        with pytest.raises(OrbitalError):
            get_network_design_matrix([], PLANAR, 100, True)

    def test_invalid_degree_arg(self):
        # test failure of a few different args for 'degree'
        for d in range(-5, 1):
            with pytest.raises(OrbitalError):
                get_network_design_matrix(self.ifgs, d, 100, True)
        for d in range(4, 7):
            with pytest.raises(OrbitalError):
                get_network_design_matrix(self.ifgs, d, 100, True)

    def test_invalid_method(self):
        # test failure of a few different args for 'method'
        for m in [None, 5, -1, -3, 45.8]:
            self.params[C.ORBITAL_FIT_METHOD] = m
            with pytest.raises(OrbitalError):
                remove_orbital_error(self.ifgs, self.params)

    def test_different_looks_raise(self):
        # different x/y looks factors should be accepted
        self.params[C.ORBITAL_FIT_LOOKS_X] = 1
        self.params[C.ORBITAL_FIT_LOOKS_Y] = 5
        try:
            remove_orbital_error(self.ifgs, self.params)
        except:
            pytest.fail

    def test_looks_as_int(self):
        self.params[C.ORBITAL_FIT_LOOKS_X] = 1.1
        self.params[C.ORBITAL_FIT_LOOKS_Y] = 5
        with pytest.raises(OrbitalError):
            remove_orbital_error(self.ifgs, self.params)
        self.params[C.ORBITAL_FIT_LOOKS_X] = 1
        self.params[C.ORBITAL_FIT_LOOKS_Y] = '5'
        with pytest.raises(OrbitalError):
            remove_orbital_error(self.ifgs, self.params)


class TestNetworkDesignMatrixTests:
    """Contains tests verifying creation of sparse network design matrix."""

    def setup_class(self):
        self.ifgs = small5_mock_ifgs()
        _add_nodata(self.ifgs)
        self.nifgs = len(self.ifgs)
        self.ncells = self.ifgs[0].num_cells
        self.date_ids = get_date_ids(self.ifgs)
        self.nepochs = len(self.date_ids)
        assert self.nepochs == 6

        for ifg in self.ifgs:
            ifg.X_SIZE = 90.0
            ifg.Y_SIZE = 89.5

    def test_planar_network_dm(self):
        ncoef = 2
        offset = False
        act = get_network_design_matrix(self.ifgs, PLANAR, 100, intercept=offset)
        assert act.shape == (self.ncells * self.nifgs, ncoef * self.nepochs)
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_planar_network_dm_offset(self):
        ncoef = 2  # NB: doesn't include offset col
        offset = True
        act = get_network_design_matrix(self.ifgs, PLANAR, 100, intercept=offset)
        assert act.shape[0] == self.ncells * self.nifgs
        assert act.shape[1] == (self.nepochs * (ncoef + offset))
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_quadratic_network_dm(self):
        ncoef = 5
        offset = False
        act = get_network_design_matrix(self.ifgs, QUADRATIC, 100, intercept=offset)
        assert act.shape == (self.ncells * self.nifgs, ncoef * self.nepochs)
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_quadratic_network_dm_offset(self):
        ncoef = 5
        offset = True
        act = get_network_design_matrix(self.ifgs, QUADRATIC, 100, intercept=offset)
        assert act.shape[0] == self.ncells * self.nifgs
        assert act.shape[1] == (self.nepochs * (ncoef + offset))
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_partcubic_network_dm(self):
        ncoef = 6
        offset = False
        act = get_network_design_matrix(self.ifgs, PART_CUBIC, 100, intercept=offset)
        assert act.shape == (self.ncells * self.nifgs, ncoef * self.nepochs)
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_partcubic_network_dm_offset(self):
        ncoef = 6
        offset = True
        act = get_network_design_matrix(self.ifgs, PART_CUBIC, 100, intercept=offset)
        assert act.shape[0] == self.ncells * self.nifgs
        assert act.shape[1] == (self.nepochs * (ncoef + offset))
        assert act.ptp() != 0
        self.check_equality(ncoef, act, self.ifgs, offset)

    def check_equality(self, ncoef, dm, ifgs, offset):
        """
        Internal test function to check subsets against network design matrix
        ncoef - base number of coefficients, without extra col for offsets
        dm - network design matrix to check the results
        ifgs - sequence of Ifg objs
        offset - boolean to include extra parameters for model offsets
        """
        deg = DEG_LOOKUP[ncoef]
        np = ncoef * self.nepochs  # index of 1st offset col

        for i, ifg in enumerate(ifgs):
            exp = unittest_dm(ifg, NETWORK_METHOD, deg, offset)
            assert exp.shape == (ifg.num_cells, ncoef)

            ib1, ib2 = [x * self.ncells for x in (i, i + 1)]  # row start/end
            jbm = (ncoef + offset) * self.date_ids[ifg.first]  # starting col index for first image
            jbs = (ncoef + offset) * self.date_ids[ifg.second]  # col start for second image
            assert_array_almost_equal(-exp, dm[ib1:ib2, jbm:jbm + ncoef])
            assert_array_almost_equal(exp, dm[ib1:ib2, jbs:jbs + ncoef])

            # ensure remaining rows/cols are zero for this ifg NOT inc offsets
            assert_array_equal(0, dm[ib1:ib2, :jbm])  # all cols leading up to first image
            assert_array_equal(0, dm[ib1:ib2, jbm + ncoef + offset:jbs])  # cols btwn mas/slv
            assert_array_equal(0, dm[ib1:ib2, jbs + ncoef + offset:np])  # to end of non offsets


# components for network correction testing
def network_correction(ifgs, deg, intercept, ml_ifgs=None, tol=1e-6):
    """
    Compares results of orbital_correction() to alternate implementation.
    deg - PLANAR, QUADRATIC or PART_CUBIC
    off - True/False to calculate correction with offsets
    """
    ncells = ifgs[0].num_cells

    if ml_ifgs:
        ml_nc = ml_ifgs[0].num_cells
        ml_data = concatenate([i.phase_data.reshape(ml_nc) for i in ml_ifgs])
        dm = get_network_design_matrix(ml_ifgs, deg, 100, intercept)[~isnan(ml_data)]
        fd = ml_data[~isnan(ml_data)].reshape((dm.shape[0], 1))
    else:
        data = concatenate([i.phase_data.reshape(ncells) for i in ifgs])
        dm = get_network_design_matrix(ifgs, deg, 100, intercept)[~isnan(data)]
        fd = data[~isnan(data)].reshape((dm.shape[0], 1))

    params = pinv(dm, tol).dot(fd)
    assert params.shape == (dm.shape[1], 1)

    # calculate forward correction
    sdm = unittest_dm(ifgs[0], NETWORK_METHOD, deg)
    ncoef = _get_num_params(deg, intercept=False)  # NB: ignore offsets for network method
    assert sdm.shape == (ncells, ncoef)
    orbs = _expand_corrections(ifgs, sdm, params, ncoef, intercept)

    # tricky: get expected result before orbital_correction() modifies ifg phase
    return [i.phase_data - orb for i, orb in zip(ifgs, orbs)]


def _expand_corrections(ifgs, dm, params, ncoef, offset):
    """
    Convenience func returns model converted to data points.
    dm: design matrix (do not filter/remove nan cells)
    params: model parameters array from pinv() * dm
    ncoef: number of model coefficients (2 planar, 5 quadratic)
    offsets: True/False to calculate correction with offsets
    """
    # NB: cannot work on singular ifgs due to date ID id/indexing requirement
    date_ids = get_date_ids(ifgs)

    corrections = []
    for ifg in ifgs:
        jbm = date_ids[ifg.first] * ncoef  # starting row index for first image
        jbs = date_ids[ifg.second] * ncoef  # row start for second image
        par = params[jbs:jbs + ncoef] - params[jbm:jbm + ncoef]

        # estimate orbital correction effects
        # corresponds to "fullorb = B*parm + offset" in orbfwd.m
        cor = dm.dot(par).reshape(ifg.phase_data.shape)

        if offset:
            off = np.ravel(ifg.phase_data - cor)
            # bring all ifgs to same base level
            cor -= nanmedian(off)

        corrections.append(cor)
    return corrections


class TestNetworkCorrectionTests:
    """Verifies orbital correction using network method and no multilooking"""

    def setup_class(cls):
        # fake some real ifg data by adding nans
        cls.ifgs = small5_mock_ifgs()
        _add_nodata(cls.ifgs)

        # use different sizes to differentiate axes results
        for ifg in cls.ifgs:
            ifg.X_SIZE = 90.0
            ifg.Y_SIZE = 89.5

        cls.nc_tol = 1e-6

    """
    this test checks that the network orbital fit will return the same
    parameters if we add a constant to every interferogram. The current
    network method actually uses the constant parameters, which are
    assigned per epoch rather than per interferogram, so this test will
    fail (and should fail).
    """

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_offset_inversion(self):
        """
        Ensure pinv(DM)*obs gives equal results given constant change to fd
        """

        def get_orbital_params():
            """Returns pseudo-inverse of the DM"""
            ncells = self.ifgs[0].num_cells
            data = concatenate([i.phase_data.reshape(ncells) for i in self.ifgs])
            dm = get_network_design_matrix(self.ifgs, PLANAR, 100, True)[~isnan(data)]
            fd = data[~isnan(data)].reshape((dm.shape[0], 1))
            return dot(pinv(dm, self.nc_tol), fd)

        tol = 1e-5
        nifgs = len(self.ifgs)
        params0 = get_orbital_params()

        # apply constant change to the observed values (fd)
        for value in [5.2, -23.5]:
            for i in self.ifgs:  # change ifgs in place
                i.phase_data += value
                assert isnan(i.phase_data).any()

            params = get_orbital_params()
            diff = params - params0
            assert (diff[:-nifgs] < tol).all()
            assert_array_almost_equal(diff[-nifgs:], value, decimal=5)

            # reset back to orig data
            for i in self.ifgs:
                i.phase_data -= value

    # These functions test full size data for orbital correction. The options
    # are separated as the ifg.phase_data arrays are modified in place, allowing
    # setUp() reset phase data between tests.

    def test_network_correction_planar(self):
        deg, intercept = PLANAR, False
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    """
    the test_network_correction_{DEGREE}_offset tests check against a method
    that fits a constant offset to each interferogram but doesn't remove it.
    This differs from the current implementation of the network correction
    so these tests will fail.
    """

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_network_correction_planar_offset(self):
        deg, intercept = PLANAR, True
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    def test_network_correction_quadratic(self):
        deg, intercept = QUADRATIC, False
        offset = intercept
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_network_correction_quadratic_offset(self):
        deg, intercept = QUADRATIC, True
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    def test_network_correction_partcubic(self):
        deg, intercept = PART_CUBIC, False
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_network_correction_partcubic_offset(self):
        deg, intercept = PART_CUBIC, True
        exp = network_correction(self.ifgs, deg, intercept)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @staticmethod
    def verify_corrections(ifgs, exp, deg, intercept):
        # checks orbital correction against unit test version
        params = dict()
        params[C.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        params[C.ORBITAL_FIT_DEGREE] = deg
        params[C.ORBITAL_FIT_LOOKS_X] = 1
        params[C.ORBITAL_FIT_LOOKS_Y] = 1
        params[C.PARALLEL] = False
        params[C.OUT_DIR] = tempfile.mkdtemp()
        params[C.ORBFIT_OFFSET] = intercept
        params[C.ORBFIT_INTERCEPT] = intercept
        params[C.ORBFIT_SCALE] = 100
        params[C.PREREAD_IFGS] = None
        mkdir_p(Path(params[C.OUT_DIR]).joinpath(C.ORB_ERROR_DIR))
        network_orbital_correction(ifgs, params)
        act = [i.phase_data for i in ifgs]
        assert_array_almost_equal(act, exp, decimal=5)


class TestNetworkCorrectionTestsMultilooking:
    'Verifies orbital correction with multilooking and network method'

    @classmethod
    def setup_class(cls):
        # fake some real ifg data by adding nans
        cls.ml_ifgs = small5_mock_ifgs()
        # 2x data of default Small mock
        cls.ifgs = small5_mock_ifgs(xs=6, ys=8)

        # use different sizes to differentiate axes results
        for ifg in cls.ifgs:
            ifg.X_SIZE = 90.0
            ifg.Y_SIZE = 89.5

        # add common nodata to all ifgs
        for i in cls.ifgs + cls.ml_ifgs:
            i.phase_data[0, :] = nan

    # These functions test multilooked data for orbital correction. The options
    # are separated as the ifg.phase_data arrays are modified in place, allowing
    # setUp() refresh phase data between tests.

    def test_mlooked_network_correction_planar(self):
        deg, intercept = PLANAR, False
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_mlooked_network_correction_planar_offset(self):
        deg, intercept = PLANAR, True
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    def test_mlooked_network_correction_quadratic(self):
        deg, intercept = QUADRATIC, False
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_mlooked_network_correction_quadratic_offset(self):
        deg, intercept = QUADRATIC, True
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    def test_mlooked_network_correction_partcubic(self):
        deg, intercept = PART_CUBIC, False
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    @pytest.mark.skip(reason="legacy test against old network method")
    def test_mlooked_network_correction_partcubic_offset(self):
        deg, intercept = PART_CUBIC, True
        exp = network_correction(self.ifgs, deg, intercept, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, intercept)

    def verify_corrections(self, ifgs, exp, deg, intercept):
        # checks orbital correction against unit test version
        params = dict()
        params[C.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        params[C.ORBITAL_FIT_DEGREE] = deg
        params[C.ORBITAL_FIT_LOOKS_X] = 1
        params[C.ORBITAL_FIT_LOOKS_Y] = 1
        params[C.PARALLEL] = False
        params[C.ORBFIT_OFFSET] = intercept
        params[C.ORBFIT_INTERCEPT] = intercept
        params[C.ORBFIT_SCALE] = 100
        params[C.PREREAD_IFGS] = None
        params[C.OUT_DIR] = tempfile.mkdtemp()
        mkdir_p(Path(params[C.OUT_DIR]).joinpath(C.ORB_ERROR_DIR))
        network_orbital_correction(ifgs, params, self.ml_ifgs)
        act = [i.phase_data for i in ifgs]
        assert_array_almost_equal(act, exp, decimal=4)


def unittest_dm(ifg, method, degree, offset=False, scale=100.0):
    '''Helper/test func to create design matrix segments. Includes handling for
    making quadratic DM segments for use in network method.
    ifg - source interferogram to model design matrix on
    method - INDEPENDENT_METHOD or NETWORK_METHOD
    degree - PLANAR, QUADRATIC or PART_CUBIC
    offset - True/False to include additional cols for offsets
    '''
    assert method in [INDEPENDENT_METHOD, NETWORK_METHOD]

    xlen = ncoef = NUM_COEF_LOOKUP[degree]
    if offset and method == INDEPENDENT_METHOD:
        ncoef += 1
    else:
        offset = False  # prevent offsets in DM sections for network method

    # NB: avoids meshgrid to prevent copying production implementation
    data = empty((ifg.num_cells, ncoef), dtype=float32)
    rows = iter(data)
    yr = range(1, ifg.nrows + 1)  # simulate meshgrid starting from 1
    xr = range(1, ifg.ncols + 1)

    xsz, ysz = [i / scale for i in [ifg.x_size, ifg.y_size]]

    if degree == PLANAR:
        for y, x in product(yr, xr):
            row = next(rows)
            row[:xlen] = [x * xsz, y * ysz]
    elif degree == QUADRATIC:
        for y, x in product(yr, xr):
            ys = y * ysz
            xs = x * xsz
            row = next(rows)
            row[:xlen] = [xs ** 2, ys ** 2, xs * ys, xs, ys]
    else:
        for y, x in product(yr, xr):
            ys = y * ysz
            xs = x * xsz
            row = next(rows)
            row[:xlen] = [xs * ys ** 2, xs ** 2, ys ** 2, xs * ys, xs, ys]

    if offset:
        data[:, -1] = 1

    return data


def get_date_ids(ifgs):
    '''
    Returns unique epoch date IDs from the given Ifgs.
    '''

    dates = []
    for ifg in ifgs:
        dates += [ifg.first, ifg.second]
    return first_second_ids(dates)


def _add_nodata(ifgs):
    """Adds some NODATA/nan cells to the small mock ifgs"""
    ifgs[0].phase_data[0, :] = nan  # 3 error cells
    ifgs[1].phase_data[2, 1:3] = nan  # 2 error cells
    ifgs[2].phase_data[3, 2:3] = nan  # 1 err
    ifgs[3].phase_data[1, 2] = nan  # 1 err
    ifgs[4].phase_data[1, 1:3] = nan  # 2 err


class TestLegacyComparisonTestsOrbfitMethod1:
    """
    This is the legacy comparison test of orbital correction functionality.
    Tests use the following config
    orbfit:        1
    orbfitmethod:  1
    orbfitdegrees: 1
    orbfitlksx:    1
    orbfitlksy:    1

    """

    @classmethod
    @pytest.fixture(autouse=True)
    def setup_class(cls, roipac_params):
        cls.params = roipac_params
        cls.BASE_DIR = cls.params[C.OUT_DIR]
        # change to orbital error correction method 1
        cls.params[C.ORBITAL_FIT_METHOD] = INDEPENDENT_METHOD
        cls.params[C.ORBITAL_FIT_LOOKS_X] = 1
        cls.params[C.ORBITAL_FIT_LOOKS_Y] = 1
        cls.params[C.PARALLEL] = False
        cls.params[C.ORBFIT_OFFSET] = True

        data_paths = [os.path.join(SML_TEST_TIF, p) for p in IFMS16]
        cls.ifg_paths = [os.path.join(cls.BASE_DIR, os.path.basename(d)) for d in data_paths]

        for d in data_paths:
            shutil.copy(d, os.path.join(cls.BASE_DIR, os.path.basename(d)))

    @classmethod
    def teardown_class(cls):
        "roipac_params fixture auto cleans"
        pass

    @pytest.mark.skipif(True, reason="Does not work anymore")
    def test_orbital_correction_legacy_equality(self):
        from pyrate import correct
        from pyrate.configuration import MultiplePaths

        multi_paths = [MultiplePaths(p, params=self.params) for p in self.ifg_paths]
        for m in multi_paths:  # cheat
            m.sampled_path = m.converted_path

        self.params[C.INTERFEROGRAM_FILES] = multi_paths
        self.params['rows'], self.params['cols'] = 2, 3
        self.params[C.ORBFIT_OFFSET] = False
        Path(self.BASE_DIR).joinpath('tmpdir').mkdir(exist_ok=True, parents=True)
        correct._copy_mlooked(self.params)
        correct._update_params_with_tiles(self.params)
        correct._create_ifg_dict(self.params)
        correct._copy_mlooked(self.params)
        pyrate.core.orbital.orb_fit_calc_wrapper(self.params)

        onlyfiles = [f for f in os.listdir(SML_TEST_LEGACY_ORBITAL_DIR)
                     if os.path.isfile(os.path.join(SML_TEST_LEGACY_ORBITAL_DIR, f))
                     and f.endswith('.csv') and f.__contains__('_method1_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(SML_TEST_LEGACY_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate([m.tmp_sampled_path for m in multi_paths]):
                ifg = Ifg(j)
                ifg.open()
                if os.path.basename(j).split('_ifg.')[0] == os.path.basename(f).split(
                        '_orb_planar_1lks_method1_geo_')[1].split('.')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data, ifg.phase_data, decimal=2)

                    # means must also be equal
                    assert np.nanmean(ifg_data) == pytest.approx(np.nanmean(ifg.phase_data), abs=1e-2)

                    # number of nans must equal
                    assert np.sum(np.isnan(ifg_data)) == np.sum(np.isnan(ifg.phase_data))
                ifg.close()

        # ensure that we have expected number of matches
        assert count == len(self.ifg_paths)

    def test_orbfit_treats_process_inputs_as_read_only(self):
        pass


class TestLegacyComparisonTestsOrbfitMethod2:
    """
    This is the legacy comparison test of orbital correction functionality.
    Tests use the following config
    orbfit:        1
    orbfitmethod:  2
    orbfitdegrees: 1
    orbfitlksx:    1
    orbfitlksy:    1

    """

    @classmethod
    def setup_class(cls):
        # change to orbital error correction method 2
        cls.params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.BASE_DIR = cls.params[C.OUT_DIR]
        cls.params[C.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        cls.params[C.ORBITAL_FIT_LOOKS_X] = 1
        cls.params[C.ORBITAL_FIT_LOOKS_Y] = 1
        cls.params[C.ORBFIT_OFFSET] = True
        cls.params[C.OUT_DIR] = cls.BASE_DIR
        data_paths = [os.path.join(SML_TEST_TIF, p) for p in small_ifg_file_list()]
        cls.new_data_paths = [os.path.join(cls.BASE_DIR, os.path.basename(d)) for d in data_paths]
        cls.params[C.INTERFEROGRAM_FILES] = [MultiplePaths(file_name=d, params=cls.params) for d in data_paths]
        for p in cls.params[C.INTERFEROGRAM_FILES]:
            p.sampled_path = p.converted_path

        # copy the files from the dir into temp dir
        for d in data_paths:
            d_copy = os.path.join(cls.BASE_DIR, os.path.basename(d))
            shutil.copy(d, d_copy)
            os.chmod(d_copy, 0o660)

        cls.headers = [roipac.roipac_header(i, cls.params) for i in cls.new_data_paths]
        cls.orb_error_dir = Path(cls.params[C.OUT_DIR]).joinpath(ORB_ERROR_DIR)
        cls.orb_error_dir.mkdir(parents=True, exist_ok=True)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.BASE_DIR, ignore_errors=True)

    def test_orbital_correction_legacy_equality_orbfit_method_2(self):
        correct._copy_mlooked(self.params)
        correct._create_ifg_dict(self.params)
        remove_orbital_error(self.new_data_paths, self.params)

        onlyfiles = [f for f in os.listdir(SML_TEST_LEGACY_ORBITAL_DIR)
                     if os.path.isfile(os.path.join(SML_TEST_LEGACY_ORBITAL_DIR, f))
                     and f.endswith('.csv') and f.__contains__('_method2_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            legacy_phase_data = np.genfromtxt(os.path.join(
                SML_TEST_LEGACY_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate(self.new_data_paths):
                if os.path.basename(j).split('_unw.')[0] == os.path.basename(f).split('_method2_')[1].split('.')[0]:
                    count += 1
                    ifg = Ifg(j)
                    ifg.open()
                    # all numbers equal
                    # Note this changed as the nodata mask in the gdal_python.gdal_average changed to nan from 0
                    # np.testing.assert_array_almost_equal(legacy_phase_data, ifg.phase_data, decimal=3)
                    # number of nans must equal
                    assert np.sum(np.isnan(legacy_phase_data)) == np.sum(np.isnan(ifg.phase_data))

        # ensure that we have expected number of matches
        assert count == len(self.new_data_paths)

    def test_orbital_error_method2_dummy(self):
        """
        does not test anything except that the method is working
        """
        # change to orbital error correction method 2
        self.params[C.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        self.params[C.ORBITAL_FIT_LOOKS_X] = 2
        self.params[C.ORBITAL_FIT_LOOKS_Y] = 2
        correct._copy_mlooked(self.params)
        correct._create_ifg_dict(self.params)
        remove_orbital_error(self.new_data_paths, self.params)

        onlyfiles = [f for f in os.listdir(SML_TEST_LEGACY_ORBITAL_DIR)
                     if os.path.isfile(os.path.join(SML_TEST_LEGACY_ORBITAL_DIR, f))
                     and f.endswith('.csv') and f.__contains__('_method2_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            legacy_phase_data = np.genfromtxt(os.path.join(SML_TEST_LEGACY_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate(self.new_data_paths):
                if os.path.basename(j).split('_unw.')[0] == os.path.basename(f).split('_method2_')[1].split('.')[0]:
                    count += 1
                    ifg = Ifg(j)
                    ifg.open()
                    # number of nans must equal
                    assert np.sum(np.isnan(legacy_phase_data)) == np.sum(np.isnan(ifg.phase_data))

        # ensure that we have expected number of matches
        assert count == len(self.new_data_paths)


# TODO: Write tests for various looks and degree combinations
# TODO: write mpi tests


class TestOrbErrorCorrectionsOnDiscReused:

    @classmethod
    def setup_class(cls):
        cls.conf = TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._create_ifg_dict(cls.params)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[C.OUT_DIR])

    def test_orb_error(self, orbfit_method, orbfit_degrees):
        self.params[C.ORBITAL_FIT_METHOD] = orbfit_method
        self.params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        multi_paths = self.params[C.INTERFEROGRAM_FILES]
        self.ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        remove_orbital_error(self.ifg_paths, self.params)

        # test_orb_errors_written
        orb_error_files = [MultiplePaths.orb_error_path(i, self.params) for i in self.ifg_paths]
        assert all(p.exists() for p in orb_error_files)

        last_mod_times = np.array([os.stat(o).st_mtime for o in orb_error_files])

        # run orbit removal again
        remove_orbital_error(self.ifg_paths, self.params)
        orb_error_files2 = [MultiplePaths.orb_error_path(i, self.params) for i in self.ifg_paths]
        # if files are written again - times will change
        last_mod_times_2 = np.array([os.stat(o).st_mtime for o in orb_error_files2])

        # test_orb_error_reused_if_params_unchanged
        assert all(a == b for a, b in zip(last_mod_times, last_mod_times_2))

        # change one of the params
        _degrees = set(C.ORB_DEGREE_NAMES.keys())
        _degrees.discard(orbfit_degrees)

        # test_orb_errors_recalculated_if_params_change
        self.params[C.ORBITAL_FIT_DEGREE] = _degrees.pop()
        import time
        time.sleep(0.1)
        remove_orbital_error(self.ifg_paths, self.params)
        orb_error_files3 = [MultiplePaths.orb_error_path(i, self.params) for i in self.ifg_paths]
        last_mod_times_3 = np.array([os.stat(o).st_mtime for o in orb_error_files3])
        assert all(a != b for a, b in zip(last_mod_times, last_mod_times_3))


class TestOrbErrorCorrectionsReappliedDoesNotChangePhaseData:

    @classmethod
    def setup_method(cls):
        cls.conf = TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._create_ifg_dict(cls.params)
        multi_paths = cls.params[C.INTERFEROGRAM_FILES]
        cls.ifg_paths = [p.tmp_sampled_path for p in multi_paths]

    @classmethod
    def teardown_method(cls):
        shutil.rmtree(cls.params[C.OUT_DIR])

    def test_orb_error_multiple_run_does_not_change_phase_data(self, orbfit_method, orbfit_degrees):
        self.params[C.ORBITAL_FIT_METHOD] = orbfit_method
        self.params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        remove_orbital_error(self.ifg_paths, self.params)
        ifgs = [Ifg(i) for i in self.ifg_paths]
        for i in ifgs:
            i.open()

        phase_prev = [i.phase_data for i in ifgs]

        # orb correct once more
        correct._copy_mlooked(self.params)
        remove_orbital_error(self.ifg_paths, self.params)

        # and again
        correct._copy_mlooked(self.params)
        remove_orbital_error(self.ifg_paths, self.params)
        ifgs = [Ifg(i) for i in self.ifg_paths]
        for i in ifgs:
            i.open()
        phase_now = [i.phase_data for i in ifgs]
        np.testing.assert_array_equal(phase_now, phase_prev)


@pytest.fixture(params=[2, 3, 4])
def orbfit_looks(request):
    x_lk = request.param
    y_lk = np.random.choice([2, 3, 4])
    return x_lk, y_lk


class TestOrbfitIndependentMethodWithMultilooking:

    @classmethod
    def setup_class(cls):
        cls.conf = TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._create_ifg_dict(cls.params)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[C.OUT_DIR])

    def test_independent_method_works_with_multilooking(self, orbfit_looks, orbfit_degrees, orbfit_method=1):
        """
        tests when multilooking is used in orbfit method 1 correction
        also tests that multilooking factors in x and y can be different
        """
        xlks, ylks = orbfit_looks
        self.params[C.ORBITAL_FIT_METHOD] = orbfit_method
        self.params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        self.params[C.ORBITAL_FIT_LOOKS_Y] = int(ylks)
        self.params[C.ORBITAL_FIT_LOOKS_X] = int(xlks)
        multi_paths = self.params[C.INTERFEROGRAM_FILES]
        self.ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        remove_orbital_error(self.ifg_paths, self.params)
        ifgs = [Ifg(p) for p in self.ifg_paths]
        for i in ifgs:
            i.open()
            assert i.shape == (72, 47)  # shape should not change


from pyrate.core.shared import cell_size


class SyntheticIfg:
    """
    This class will generate a mock interferogram whose signal consists entirely
    of a synthetic orbital error. The orbital error is generated as a 2D
    polynomial signal with zero noise.
    """

    def __init__(self, orbfit_degrees):
        self.x_step = 0.001388888900000  # pixel size - same as cropA
        self.y_step = 0.001388888900000
        self.nrows = 100
        self.ncols = 100
        self.num_cells = self.nrows * self.ncols
        self.is_open = False
        self.orbfit_degrees = orbfit_degrees
        self.first = None
        self.second = None
        self._phase_data = None
        self._phase_data_first = None
        self._phase_data_second = None
        self.y_first = 0
        self.x_first = 0
        self.add_geographic_data()

    def add_geographic_data(self):
        """
        Determine and add geographic data to object
        """
        # add some geographic data
        self.x_centre = int(self.ncols / 2)
        self.y_centre = int(self.nrows / 2)
        self.lat_centre = self.y_first + (self.y_step * self.y_centre)
        self.long_centre = self.x_first + (self.x_step * self.x_centre)
        # use cell size from centre of scene
        self.x_size, self.y_size = cell_size(self.lat_centre, self.long_centre, self.x_step, self.y_step)

    @property
    def phase_data(self):
        """
        Returns phase band as an array.
        """
        if self._phase_data is None:
            self.open()
        return self._phase_data

    def open(self):
        x, y = np.meshgrid(np.arange(self.nrows) * self.x_step, np.arange(self.ncols) * self.y_step)
        x += self.x_step
        y += self.y_step

        # define some random coefficients, different for each date
        x_slope, y_slope, x2_slope, y2_slope, x_y_slope, x_y2_slope, const = np.ravel(np.random.rand(1, 7))
        x_slope_, y_slope_, x2_slope_, y2_slope_, x_y_slope_, x_y2_slope_, const_ = np.ravel(np.random.rand(1, 7))

        # compute the 2D polynomial separately for first and second dates
        self._phase_data_first = x_slope * x + y_slope * y + const  # planar
        self._phase_data_second = x_slope_ * x + y_slope_ * y + const_  # planar
        if self.orbfit_degrees == QUADRATIC:
            self._phase_data_first += x2_slope * x ** 2 + y2_slope * y ** 2 + x_y_slope * x * y
            self._phase_data_second += x2_slope_ * x ** 2 + y2_slope_ * y ** 2 + x_y_slope_ * x * y
        elif self.orbfit_degrees == PART_CUBIC:
            self._phase_data_first += x2_slope * x ** 2 + y2_slope * y ** 2 + x_y_slope * x * y + \
                                      x_y2_slope * x * (y ** 2)
            self._phase_data_second += x2_slope_ * x ** 2 + y2_slope_ * y ** 2 + x_y_slope_ * x * y + \
                                       x_y2_slope_ * x * (y ** 2)

        # combine orbit error for first and second dates to give synthetic phase data for this ifg
        self._phase_data = self._phase_data_first - self._phase_data_second
        self.is_open = True


from pyrate.core.gdal_python import _gdalwarp_width_and_height
from pyrate.core.orbital import __orb_inversion


# helper function to multilook a synthetic ifg with gdal
def mlk_ifg(ifg, nlooks):
    src = gdal.GetDriverByName('MEM').Create('', ifg.ncols, ifg.nrows, 1, gdalconst.GDT_Float32)
    gt = (0, ifg.x_step, 0, 0, 0, ifg.y_step)
    src.SetGeoTransform(gt)
    src.GetRasterBand(1).WriteArray(ifg.phase_data)
    resampled_gt = (0, ifg.x_step * nlooks, 0, 0, 0, ifg.y_step * nlooks)
    min_x, min_y = 0, 0
    max_x, max_y = ifg.x_step * ifg.ncols, ifg.y_step * ifg.nrows

    px_height, px_width = _gdalwarp_width_and_height(max_x, max_y, min_x, min_y, resampled_gt)

    dst = gdal.GetDriverByName('MEM').Create('', px_height, px_width, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform(resampled_gt)

    gdal.ReprojectImage(src, dst, '', '', gdal.GRA_Average)

    mlooked = Ifg(dst)
    mlooked.first = ifg.first
    mlooked.second = ifg.second
    return mlooked


@pytest.fixture(params=[1, 2, 3, 4])
def orb_lks(request):
    return request.param


def test_single_synthetic_ifg_independent_method(orbfit_degrees, orb_lks, ifg=None):
    """
    These tests are checking that perfect orbital errors, those matching the assumed orbital error model, can be
    completely removed by the independent orbital correction with and without multilooking.

    These tests also prove that orbital error estimates using orbfit multilooking is a valid approach and matches the
    modelled error with acceptable numerical accuracy, with the accuracy depending on the multilooking factor used.

    These tests also prove that the orbital error parameters can be approximated using multilooking in the
    independent method.
    """
    if ifg is None:
        ifg = SyntheticIfg(orbfit_degrees)
    fullres_dm = get_design_matrix(ifg, orbfit_degrees, intercept=True, scale=1)

    m_looked_ifg = mlk_ifg(ifg, orb_lks)

    mlooked_phase = np.reshape(m_looked_ifg.phase_data, m_looked_ifg.num_cells)
    mlooked_dm = get_design_matrix(m_looked_ifg, orbfit_degrees, intercept=True, scale=1)

    orb_corr = __orb_correction(fullres_dm, mlooked_dm, ifg.phase_data, mlooked_phase, offset=True)
    if orb_lks == 1:
        assert_array_almost_equal(fullres_dm, mlooked_dm)
        decimal = 4
    else:
        decimal = 2
    assert_array_almost_equal(ifg.phase_data, orb_corr, decimal=decimal)


@pytest.mark.slow
@pytest.mark.skipif((not PY37GDAL302), reason="Only run in one CI env")
def test_set_synthetic_ifgs_independent_method(mexico_cropa_params, orbfit_degrees, orb_lks):
    """
    Test that the independent method can generate a set of orbital corrections
    that matches a set of synthetic ifg for a range of multi-look factors and
    polynomial degrees.
    """
    # Use the CropA ifg network configuration
    ifgs = [Ifg(i.converted_path) for i in mexico_cropa_params[C.INTERFEROGRAM_FILES]]
    for i in ifgs:
        i.open()
        test_ifg = SyntheticIfg(orbfit_degrees)
        test_single_synthetic_ifg_independent_method(orbfit_degrees, orb_lks, test_ifg)


# an in-memory "open" interferogram that we can pass to top-level orb correction methods
class FakeIfg:
    def __init__(self, orbfit_deg, model_params, date_first, date_second):
        self.x_step = 0.001388888900000  # pixel size - same as cropA
        self.y_step = 0.001388888900000
        self.nrows = 100
        self.ncols = 100
        self.num_cells = self.nrows * self.ncols
        self.is_open = False
        self.orbfit_degrees = orbfit_deg
        self.model_params = model_params
        self.first = date_first
        self.second = date_second
        self._phase_data = None
        self.y_first = 0
        self.x_first = 0
        self.nan_fraction = 0
        self.add_geographic_data()

    def add_geographic_data(self):
        """
        Determine and add geographic data to object
        """
        # add some geographic data
        self.x_centre = int(self.ncols / 2)
        self.y_centre = int(self.nrows / 2)
        self.lat_centre = self.y_first + (self.y_step * self.y_centre)
        self.long_centre = self.x_first + (self.x_step * self.x_centre)
        # use cell size from centre of scene
        self.x_size, self.y_size = cell_size(self.lat_centre, self.long_centre, self.x_step, self.y_step)

    @property
    def phase_data(self):
        """
        Returns phase band as an array.
        """
        if self._phase_data is None:
            self.open()
        return self._phase_data

    def open(self):
        x, y = np.meshgrid(np.arange(self.nrows) * self.x_step, np.arange(self.ncols) * self.y_step)
        x += self.x_step
        y += self.y_step

        # use provided coefficients
        if self.orbfit_degrees == PLANAR:
            mx, my = self.model_params
            self._phase_data = mx * x + my * y
        elif self.orbfit_degrees == QUADRATIC:
            mx, my, mx2, my2, mxy = self.model_params
            self._phase_data = mx * x + my * y + mx2 * x ** 2 + my2 * y ** 2 + mxy * x * y
        else:
            mx, my, mx2, my2, mxy, mxy2 = self.model_params
            self._phase_data = mx * x + my * y + mx2 * x ** 2 + my2 * y ** 2 + mxy * x * y + mxy2 * x * y ** 2

        self.is_open = True


# tests for network method to recover synthetic orbital error
class SyntheticNetwork:
    """
    This class will generate a network of synthetic ifgs, based on
    orbital errors for each epoch. The signal will be purely from the synthetic
    orbital error with no noise.
    """

    def __init__(self, orbfit_deg, epochs, network, model_params):
        """
        orbfit_deg: synthesise ifgs with planar, quadratic, or part cubic
        orbit error models.
        epochs: list of epoch dates in the network
        network: list of lists, spec of the ifgs to generate for each epoch as
        primary
        model_params: list of iterable - model parameters of correct degree for
        each epoch
        """
        ifgs = []
        for i, e1 in enumerate(epochs):
            for j in network[i]:
                ifg_err_model = [model_params[j][k] - model_params[i][k] for k in range(len(model_params[0]))]
                ifgs.append(FakeIfg(orbfit_deg, ifg_err_model, e1, epochs[j]))
        self.ifgs = ifgs
        self.epochs = epochs

@pytest.mark.skip(reason="test is non-deterministic due to float calculation errors in array comparison")
def test_synthetic_network_correction(orbfit_degrees, orb_lks):
    epochs = [
        date(2000, 1, 1),
        date(2000, 1, 13),
        date(2000, 1, 25),
        date(2000, 2, 6),
        date(2000, 2, 18),
        date(2000, 3, 1)
        ]
    # start with the network as a connected tree so mst does nothing
    network = [[2], [2], [3], [4, 5], [], []]
    # six sets of model parameters - one for each epoch
    model_params = [[-1, 1, -1, 1, -1, 1],
                    [0, 1, 2, 3, 4, 5],
                    [5, 4, 3, 2, 1, 0],
                    [3, 6, 9, 6, 3, 0],
                    [9, 4, 1, 0, 1, 4],
                    [1, 1, 1, 1, 1, 1]]
    if orbfit_degrees == PLANAR:
        nparam = 2
    elif orbfit_degrees == QUADRATIC:
        nparam = 5
    else:
        nparam = 6
    model_params = [mi[:nparam] for mi in model_params]

    # network method uses a hard coded scale of 100
    scale = 100

    syn_data = SyntheticNetwork(orbfit_degrees, epochs, network, model_params)
#    id_dict = {date: i for date, i in enumerate(epochs)}
    nepochs = len(epochs)

    mlk_ifgs = [mlk_ifg(ifg, orb_lks) for ifg in syn_data.ifgs]

    coeffs = calc_network_orb_correction(mlk_ifgs, orbfit_degrees, scale, nepochs, intercept=True)

    # reconstruct correction
    reconstructed = []
    # ifgs are built with lat/long metadata,
    # orbfit modelling is done with metres coordinates
    csx = syn_data.ifgs[0].x_size
    csy = syn_data.ifgs[0].y_size
    x, y = (coord + 1 for coord in np.meshgrid(np.arange(100, dtype=float), np.arange(100, dtype=float)))
    x *= csx
    y *= csy
    x /= scale
    y /= scale

    for i, js in enumerate(network):
        for j in js:
            cpair = [cj - ci for ci, cj in zip(coeffs[i], coeffs[j])]
            if orbfit_degrees == PLANAR:
                reconstructed.append(cpair[0] * x + cpair[1] * y)
            elif orbfit_degrees == QUADRATIC:
                reconstructed.append(cpair[0] * x ** 2 + cpair[1] * y ** 2 + cpair[2] * x * y \
                                     + cpair[3] * x + cpair[4] * y)
            else:
                reconstructed.append(cpair[0] * x * y ** 2 + cpair[1] * x ** 2 + cpair[2] * y ** 2 + \
                                     cpair[3] * x * y + cpair[4] * x + cpair[5] * y)

    for orig, recon in zip(syn_data.ifgs, reconstructed):
        assert_array_almost_equal(orig.phase_data, recon, decimal=2)


def test_orbital_inversion():
    """Small unit to test the application of numpy pseudoinverse"""
    A = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
    d = np.array([2, 4, 3])
    exp = np.array([1.5, 0.5, 2.5])
    res = __orb_inversion(A, d)
    assert_array_almost_equal(res, exp, decimal=9)
