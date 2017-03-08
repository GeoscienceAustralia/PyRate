#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
import unittest
from itertools import product
from numpy import empty, dot, concatenate, float32
from numpy import nan, isnan, array
from os.path import join

import numpy as np
from numpy.linalg import pinv, inv
from numpy.testing import assert_array_equal, assert_array_almost_equal
from scipy.linalg import lstsq

import pyrate.orbital
import tests.common
from .common import sydney5_mock_ifgs, MockIfg
from pyrate import algorithm
from pyrate import config as cf
from pyrate.orbital import INDEPENDENT_METHOD, NETWORK_METHOD, PLANAR, \
    QUADRATIC, PART_CUBIC
from pyrate.orbital import OrbitalError, orbital_correction
from pyrate.orbital import get_design_matrix, get_network_design_matrix, \
    _get_num_params
from pyrate.shared import Ifg
from pyrate.shared import nanmedian
from tests import common
from tests.common import SYD_TEST_DIR
from tests.common import SYD_TEST_MATLAB_ORBITAL_DIR
from tests.common import SYD_TEST_TIF, sydney_data_setup
from tests.common import sydney_ifg_file_list

DEG_LOOKUP = {
    2: PLANAR,
    5: QUADRATIC,
    6: PART_CUBIC}

NUM_COEF_LOOKUP = {
    PLANAR: 2,
    QUADRATIC: 5,
    PART_CUBIC: 6}


class SingleDesignMatrixTests(unittest.TestCase):
    """
    Tests to verify correctness of basic planar & quadratic design matrices or
    DMs. This class serves two purposes, ensuring the independent method DMs are
    produced correctly. Secondly, these indivdual DMs are subsets of the larger
    DM 'grid' required for the networked orbital correction method.
    """

    def setUp(self):
        # faked cell sizes
        self.xs = 0.75
        self.ys = 0.8
        self.ifg = Ifg(join(SYD_TEST_TIF, 'geo_060619-061002_unw.tif'))
        self.ifg.open()
        self.ifg.nodata_value = 0

        self.m = MockIfg(self.ifg, 3, 4)
        self.m.x_size = self.xs
        self.m.y_size = self.ys

    # tests for planar model

    def test_create_planar_dm(self):
        offset = False
        act = get_design_matrix(self.m, PLANAR, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 2))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset)
        assert_array_equal(act, exp)

    def test_create_planar_dm_offsets(self):
        offset = True
        act = get_design_matrix(self.m, PLANAR, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 3))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PLANAR, offset)
        assert_array_almost_equal(act, exp)


    # tests for quadratic model

    def test_create_quadratic_dm(self):
        offset = False
        act = get_design_matrix(self.m, QUADRATIC, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 5))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset)
        assert_array_equal(act, exp)

    def test_create_quadratic_dm_offsets(self):
        offset = True
        act = get_design_matrix(self.m, QUADRATIC, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 6))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, QUADRATIC, offset)
        assert_array_equal(act, exp)

    # tests for partial cubic model

    def test_create_partcubic_dm(self):
        offset = False
        act = get_design_matrix(self.m, PART_CUBIC, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 6))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset)
        assert_array_equal(act, exp)


    def test_create_partcubic_dm_offsets(self):
        offset = True
        act = get_design_matrix(self.m, PART_CUBIC, offset)
        self.assertEqual(act.shape, (self.m.num_cells, 7))
        exp = unittest_dm(self.m, INDEPENDENT_METHOD, PART_CUBIC, offset)
        assert_array_equal(act, exp)


    # tests for unittest_dm() assuming network method

    def test_create_planar_dm_network(self):
        # networked method planar version should not have offsets col
        ncol_exp = 2
        exp = unittest_dm(self.m, NETWORK_METHOD, PLANAR, False)
        self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
        exp2 = unittest_dm(self.m, NETWORK_METHOD, PLANAR, True)
        self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
        assert_array_equal(exp, exp2)

    def test_create_quadratic_dm_network(self):
        # quadratic version with networked method does not have offsets col
        ncol_exp = 5
        exp = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, False)
        self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
        exp2 = unittest_dm(self.m, NETWORK_METHOD, QUADRATIC, True)
        self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
        assert_array_equal(exp, exp2)

    def test_create_partcubic_dm_network(self):
        # partial cubic version with networked method does not have offsets col
        ncol_exp = 6
        exp = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, False)
        self.assertEqual(exp.shape, (self.m.num_cells, ncol_exp))
        exp2 = unittest_dm(self.m, NETWORK_METHOD, PART_CUBIC, True)
        self.assertEqual(exp2.shape, (self.m.num_cells, ncol_exp))
        assert_array_equal(exp, exp2)


class IndependentCorrectionTests(unittest.TestCase):
    """Test cases for the orbital correction component of PyRate."""

    def setUp(self):
        self.ifgs = sydney5_mock_ifgs()
        _add_nodata(self.ifgs)

        for ifg in self.ifgs:
            ifg.x_size = 90.0
            ifg.y_size = 89.5
            ifg.open()

    def alt_orbital_correction(self, ifg, deg, offset):
        # almost an exact translation of MATLAB version
        data = ifg.phase_data.reshape(ifg.num_cells)
        dm = get_design_matrix(ifg, deg, offset)[~isnan(data)]
        fd = data[~isnan(data)].reshape((dm.shape[0], 1))

        dmt = dm.T
        invNbb = inv(dmt.dot(dm))
        orbparams = invNbb.dot(dmt.dot(fd))
        alt_params = lstsq(dm, fd)[0]
        # FIXME: precision
        assert_array_almost_equal(orbparams, alt_params, decimal=2)

        dm2 = get_design_matrix(ifg, deg, offset)

        if offset:
            fullorb = np.reshape(np.dot(dm2[:, :-1], orbparams[:-1]),
                             ifg.phase_data.shape)
        else:
            fullorb = np.reshape(np.dot(dm2, orbparams), ifg.phase_data.shape)

        offset_removal = nanmedian(
            np.reshape(ifg.phase_data - fullorb, (1, -1)))
        fwd_correction = fullorb - offset_removal
        # ifg.phase_data -= (fullorb - offset_removal)
        return ifg.phase_data - fwd_correction

    def check_correction(self, degree, method, offset, decimal=2):
        orig = array([c.phase_data.copy() for c in self.ifgs])
        exp = [self.alt_orbital_correction(i, degree, offset) for i in self.ifgs]
        params = dict()
        params[cf.ORBITAL_FIT_METHOD] = method
        params[cf.ORBITAL_FIT_DEGREE] = degree
        params[cf.PARALLEL] = False
        params[cf.NO_DATA_VALUE] = 0
        params[cf.NAN_CONVERSION] = False
        for i in self.ifgs:
            i.mm_converted = True
        orbital_correction(self.ifgs, params, None, offset)
        corrected = array([c.phase_data for c in self.ifgs])

        self.assertFalse((orig == corrected).all())
        self.check_results(self.ifgs, orig)  # test shape, data is non zero

        # FIXME: is decimal=2 close enough?
        for i, (e, a) in enumerate(zip(exp, corrected)):
            assert_array_almost_equal(e, a, decimal=decimal)

    def check_results(self, ifgs, corrections):
        """Helper method for result verification"""
        for i, c in zip(ifgs, corrections):
            ys, xs = c.shape
            self.assertEqual(i.nrows, ys)
            self.assertEqual(i.ncols, xs)

            # ensure there is real data
            self.assertFalse(isnan(i.phase_data).all())
            self.assertFalse(isnan(c).all())
            self.assertTrue(c.ptp() != 0)  # ensure range of values in grid

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


class ErrorTests(unittest.TestCase):
    """Tests for the networked correction method"""

    def test_invalid_ifgs_arg(self):
        # min requirement is 1 ifg, can still subtract one epoch from the other
        self.assertRaises(OrbitalError, get_network_design_matrix, [], PLANAR, True)

    def test_invalid_degree_arg(self):
        # test failure of a few different args for 'degree'
        ifgs = sydney5_mock_ifgs()
        for d in range(-5, 1):
            self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)
        for d in range(4, 7):
            self.assertRaises(OrbitalError, get_network_design_matrix, ifgs, d, True)

    def test_invalid_method(self):
        # test failure of a few different args for 'method'
        ifgs = sydney5_mock_ifgs()
        params = dict()
        params[cf.ORBITAL_FIT_DEGREE] = PLANAR
        params[cf.PARALLEL] = False
        for m in [None, 5, -1, -3, 45.8]:
            params[cf.ORBITAL_FIT_METHOD] = m
            self.assertRaises(OrbitalError, orbital_correction, ifgs, params, None)

    def test_multilooked_ifgs_arg(self):
        # check some bad args for network method with multilooked ifgs
        ifgs = sydney5_mock_ifgs()
        args = [[None, None, None, None, None], ["X"] * 5]
        params = dict()
        params[cf.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        params[cf.PARALLEL] = False
        params[cf.ORBITAL_FIT_DEGREE] = PLANAR
        for a in args:
            args = (ifgs, params, a)
            self.assertRaises(OrbitalError, orbital_correction, *args)

        # ensure failure if # ifgs doesn't match # mlooked ifgs
        args = (ifgs, params, ifgs[:4])
        self.assertRaises(OrbitalError, orbital_correction, *args)


class NetworkDesignMatrixTests(unittest.TestCase):
    """Contains tests verifying creation of sparse network design matrix."""

    def setUp(self):
        self.ifgs = sydney5_mock_ifgs()
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
        act = get_network_design_matrix(self.ifgs, PLANAR, offset)
        self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
        self.assertNotEqual(act.ptp(), 0)
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_planar_network_dm_offset(self):
        ncoef = 2 # NB: doesn't include offset col
        offset = True
        act = get_network_design_matrix(self.ifgs, PLANAR, offset)
        self.assertEqual(act.shape[0], self.ncells * self.nifgs)
        self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
        self.assertNotEqual(act.ptp(), 0)
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_quadratic_network_dm(self):
        ncoef = 5
        offset = False
        act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
        self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
        self.assertNotEqual(act.ptp(), 0)
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_quadratic_network_dm_offset(self):
        ncoef = 5
        offset = True
        act = get_network_design_matrix(self.ifgs, QUADRATIC, offset)
        self.assertEqual(act.shape[0], self.ncells * self.nifgs)
        self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
        self.assertNotEqual(act.ptp(), 0)
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_partcubic_network_dm(self):
        ncoef = 6
        offset = False
        act = get_network_design_matrix(self.ifgs, PART_CUBIC, offset)
        self.assertEqual(act.shape, (self.ncells * self.nifgs, ncoef * self.nepochs))
        self.assertNotEqual(act.ptp(), 0)
        self.check_equality(ncoef, act, self.ifgs, offset)

    def test_partcubic_network_dm_offset(self):
        ncoef = 6
        offset = True
        act = get_network_design_matrix(self.ifgs, PART_CUBIC, offset)
        self.assertEqual(act.shape[0], self.ncells * self.nifgs)
        self.assertEqual(act.shape[1], (self.nepochs * ncoef) + self.nifgs)
        self.assertNotEqual(act.ptp(), 0)
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
        np = ncoef * self.nepochs # index of 1st offset col

        for i, ifg in enumerate(ifgs):
            exp = unittest_dm(ifg, NETWORK_METHOD, deg, offset)
            self.assertEqual(exp.shape, (ifg.num_cells, ncoef))

            # NB: this is Hua Wang's MATLAB code slightly modified for Py
            ib1, ib2 = [x * self.ncells for x in (i, i+1)] # row start/end
            jbm = ncoef * self.date_ids[ifg.master] # starting col index for master
            jbs = ncoef * self.date_ids[ifg.slave] # col start for slave
            assert_array_almost_equal(-exp, dm[ib1:ib2, jbm:jbm+ncoef])
            assert_array_almost_equal( exp, dm[ib1:ib2, jbs:jbs+ncoef])

            # ensure remaining rows/cols are zero for this ifg NOT inc offsets
            assert_array_equal(0, dm[ib1:ib2, :jbm]) # all cols leading up to master
            assert_array_equal(0, dm[ib1:ib2, jbm + ncoef:jbs]) # cols btwn mas/slv
            assert_array_equal(0, dm[ib1:ib2, jbs + ncoef:np]) # to end of non offsets

            # check offset cols for 1s and 0s
            if offset is True:
                ip1 = i + np # offset column index
                assert_array_equal(1, dm[ib1:ib2, ip1])
                assert_array_equal(0, dm[ib1:ib2, np:ip1]) # cols before offset col
                assert_array_equal(0, dm[ib1:ib2, ip1 + 1:]) # cols after offset col


# components for network correction testing
def network_correction(ifgs, deg, off, ml_ifgs=None, tol=1e-6):
    """
    Compares results of orbital_correction() to alternate implementation.
    deg - PLANAR, QUADRATIC or PART_CUBIC
    off - True/False to calculate correction with offsets
    """
    ncells = ifgs[0].num_cells

    if ml_ifgs:
        ml_nc = ml_ifgs[0].num_cells
        ml_data = concatenate([i.phase_data.reshape(ml_nc) for i in ml_ifgs])
        dm = get_network_design_matrix(ml_ifgs, deg, off)[~isnan(ml_data)]
        fd = ml_data[~isnan(ml_data)].reshape((dm.shape[0], 1))
    else:
        data = concatenate([i.phase_data.reshape(ncells) for i in ifgs])
        dm = get_network_design_matrix(ifgs, deg, off)[~isnan(data)]
        fd = data[~isnan(data)].reshape((dm.shape[0], 1))

    params = pinv(dm, tol).dot(fd)
    assert params.shape == (dm.shape[1], 1)

    # calculate forward correction
    sdm = unittest_dm(ifgs[0], NETWORK_METHOD, deg)
    ncoef = _get_num_params(deg, offset=False)  # NB: ignore offsets for network method
    assert sdm.shape == (ncells, ncoef)
    orbs = _expand_corrections(ifgs, sdm, params, ncoef, off)

    # tricky: get expected result before orbital_correction() modifies ifg phase
    return [i.phase_data - orb for i, orb in zip(ifgs, orbs)]


def _expand_corrections(ifgs, dm, params, ncoef, offsets):
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
        jbm = date_ids[ifg.master] * ncoef # starting row index for master
        jbs = date_ids[ifg.slave] * ncoef # row start for slave
        par = params[jbs:jbs + ncoef] - params[jbm:jbm + ncoef]

        # estimate orbital correction effects
        # corresponds to "fullorb = B*parm + offset" in orbfwd.m
        cor = dm.dot(par).reshape(ifg.phase_data.shape)

        if offsets:
            off = np.ravel(ifg.phase_data - cor)
            # bring all ifgs to same base level
            cor -= nanmedian(off)

        corrections.append(cor)
    return corrections


class NetworkCorrectionTests(unittest.TestCase):
    """Verifies orbital correction using network method and no multilooking"""

    def setUp(self):
        # fake some real ifg data by adding nans
        self.ifgs = sydney5_mock_ifgs()
        _add_nodata(self.ifgs)

        # use different sizes to differentiate axes results
        for ifg in self.ifgs:
            ifg.X_SIZE = 90.0
            ifg.Y_SIZE = 89.5

        self.nc_tol = 1e-6

    def test_offset_inversion(self):
        """
        Ensure pinv(DM)*obs gives equal results given constant change to fd
        """
        def get_orbital_params():
            """Returns pseudo-inverse of the DM"""
            ncells = self.ifgs[0].num_cells
            data = concatenate([i.phase_data.reshape(ncells) for i in self.ifgs])
            dm = get_network_design_matrix(self.ifgs, PLANAR, True)[~isnan(data)]
            fd = data[~isnan(data)].reshape((dm.shape[0], 1))
            return dot(pinv(dm, self.nc_tol), fd)

        tol = 1e-5
        nifgs = len(self.ifgs)
        params0 = get_orbital_params()

        # apply constant change to the observed values (fd)
        for value in [5.2, -23.5]:
            for i in self.ifgs: # change ifgs in place
                i.phase_data += value
                self.assertTrue(isnan(i.phase_data).any())

            params = get_orbital_params()
            diff = params - params0
            self.assertTrue((diff[:-nifgs] < tol).all())
            assert_array_almost_equal(diff[-nifgs:], value, decimal=5)

            # reset back to orig data
            for i in self.ifgs:
                i.phase_data -= value

    # These functions test full size data for orbital correction. The options
    # are separated as the ifg.phase_data arrays are modified in place, allowing
    # setUp() reset phase data between tests.

    def test_network_correction_planar(self):
        deg, offset = PLANAR, False
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_network_correction_planar_offset(self):
        deg, offset = PLANAR, True
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_network_correction_quadratic(self):
        deg, offset = QUADRATIC, False
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_network_correction_quadratic_offset(self):
        deg, offset = QUADRATIC, True
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_network_correction_partcubic(self):
        deg, offset = PART_CUBIC, False
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_network_correction_partcubic_offset(self):
        deg, offset = PART_CUBIC, True
        exp = network_correction(self.ifgs, deg, offset)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    @staticmethod
    def verify_corrections(ifgs, exp, deg, offset):
        # checks orbital correction against unit test version
        params = dict()
        params[cf.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        params[cf.ORBITAL_FIT_DEGREE] = deg
        params[cf.PARALLEL] = False
        orbital_correction(ifgs, params, None, offset)
        act = [i.phase_data for i in ifgs]
        assert_array_almost_equal(act, exp, decimal=5)


class NetworkCorrectionTestsMultilooking(unittest.TestCase):
    'Verifies orbital correction with multilooking and network method'

    def setUp(self):
        # fake some real ifg data by adding nans
        self.ml_ifgs = sydney5_mock_ifgs()
        # 2x data of default Syd mock
        self.ifgs = sydney5_mock_ifgs(xs=6, ys=8)

        # use different sizes to differentiate axes results
        for ifg in self.ifgs:
            ifg.X_SIZE = 90.0
            ifg.Y_SIZE = 89.5

        # add common nodata to all ifgs
        for i in self.ifgs + self.ml_ifgs:
            i.phase_data[0,:] = nan

    # These functions test multilooked data for orbital correction. The options
    # are separated as the ifg.phase_data arrays are modified in place, allowing
    # setUp() refresh phase data between tests.

    def test_mlooked_network_correction_planar(self):
        deg, offset = PLANAR, False
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_mlooked_network_correction_planar_offset(self):
        deg, offset = PLANAR, True
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_mlooked_network_correction_quadratic(self):
        deg, offset = QUADRATIC, False
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_mlooked_network_correction_quadratic_offset(self):
        deg, offset = QUADRATIC, True
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_mlooked_network_correction_partcubic(self):
        deg, offset = PART_CUBIC, False
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def test_mlooked_network_correction_partcubic_offset(self):
        deg, offset = PART_CUBIC, True
        exp = network_correction(self.ifgs, deg, offset, self.ml_ifgs)
        self.verify_corrections(self.ifgs, exp, deg, offset)

    def verify_corrections(self, ifgs, exp, deg, offset):
        # checks orbital correction against unit test version
        params = dict()
        params[cf.ORBITAL_FIT_METHOD] = NETWORK_METHOD
        params[cf.ORBITAL_FIT_DEGREE] = deg
        params[cf.PARALLEL] = False
        orbital_correction(ifgs, params, self.ml_ifgs, offset)
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
        offset = False # prevent offsets in DM sections for network method

    # NB: avoids meshgrid to prevent copying production implementation
    data = empty((ifg.num_cells, ncoef), dtype=float32)
    rows = iter(data)
    yr = range(1, ifg.nrows+1)  # simulate meshgrid starting from 1
    xr = range(1, ifg.ncols+1)

    xsz, ysz = [i/scale for i in [ifg.x_size, ifg.y_size]]

    if degree == PLANAR:
        for y, x in product(yr, xr):
            row = next(rows)
            row[:xlen] = [x * xsz, y * ysz]
    elif degree == QUADRATIC:
        for y, x in product(yr, xr):
            ys = y * ysz
            xs = x * xsz
            row = next(rows)
            row[:xlen] = [xs**2, ys**2, xs*ys, xs, ys]
    else:
        for y, x in product(yr, xr):
            ys = y * ysz
            xs = x * xsz
            row = next(rows)
            row[:xlen] = [xs*ys**2, xs**2, ys**2, xs*ys, xs, ys]

    if offset:
        data[:, -1] = 1

    return data


def get_date_ids(ifgs):
    '''
    Returns unique master/slave date IDs from the given Ifgs.
    '''

    dates = []
    for ifg in ifgs:
        dates += [ifg.master, ifg.slave]
    return algorithm.master_slave_ids(dates)


def _add_nodata(ifgs):
    """Adds some NODATA/nan cells to the sydney mock ifgs"""
    ifgs[0].phase_data[0, :] = nan # 3 error cells
    ifgs[1].phase_data[2, 1:3] = nan # 2 error cells
    ifgs[2].phase_data[3, 2:3] = nan # 1 err
    ifgs[3].phase_data[1, 2] = nan # 1 err
    ifgs[4].phase_data[1, 1:3] = nan # 2 err


class MatlabComparisonTestsOrbfitMethod1(unittest.TestCase):
    """
    This is the matlab comparison test of orbital correction functionality.
    Tests use the following config
    orbfit:        1
    orbfitmethod:  2
    orbfitdegrees: 1
    orbfitlksx:    2
    orbfitlksy:    2

    """

    def setUp(self):
        self.BASE_DIR = tempfile.mkdtemp()

        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))

        # change to orbital error correction method 1
        self.params[cf.ORBITAL_FIT_METHOD] = 1
        self.params[cf.ORBITAL_FIT_LOOKS_X] = 2
        self.params[cf.ORBITAL_FIT_LOOKS_Y] = 2
        self.params[cf.PARALLEL] = False

        data_paths = [os.path.join(SYD_TEST_TIF, p) for p in common.IFMS16]
        self.ifg_paths = [os.path.join(self.BASE_DIR, os.path.basename(d))
                          for d in data_paths]

        for d in data_paths:
            shutil.copy(d, os.path.join(self.BASE_DIR, os.path.basename(d)))

    def tearDown(self):
        shutil.rmtree(self.BASE_DIR)

    def test_orbital_correction_matlab_equality(self):
        from pyrate.scripts import run_pyrate

        run_pyrate.orb_fit_calc(self.ifg_paths, self.params)

        onlyfiles = [f for f in os.listdir(SYD_TEST_MATLAB_ORBITAL_DIR)
            if os.path.isfile(os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, f))
            and f.endswith('.csv') and f.__contains__('_method1_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                SYD_TEST_MATLAB_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifg_paths):
                ifg = Ifg(j)
                ifg.open()
                if os.path.basename(j).split('_unw.')[0] == \
                        os.path.basename(f).split(
                            '_orb_planar_1lks_method1_')[1].split('.')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        ifg.phase_data, decimal=2)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data),
                        np.nanmean(ifg.phase_data), places=2)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                                np.sum(np.isnan(ifg.phase_data)))
                ifg.close()

        # ensure that we have expected number of matches
        self.assertEqual(count, len(self.ifg_paths))


class MatlabComparisonTestsOrbfitMethod2(unittest.TestCase):
    """
    This is the matlab comparison test of orbital correction functionality.
    Tests use the following config
    orbfit:        1
    orbfitmethod:  2
    orbfitdegrees: 1
    orbfitlksx:    1
    orbfitlksy:    1

    """
    def setUp(self):
        self.BASE_DIR = tempfile.mkdtemp()

        self.params = cf.get_config_params(
            os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))

        # change to orbital error correction method 2
        self.params[cf.ORBITAL_FIT_METHOD] = 2
        self.params[cf.ORBITAL_FIT_LOOKS_X] = 1
        self.params[cf.ORBITAL_FIT_LOOKS_Y] = 1

        data_paths = [os.path.join(SYD_TEST_TIF, p) for p in
                      sydney_ifg_file_list()]
        self.new_data_paths = [os.path.join(self.BASE_DIR, os.path.basename(d))
                          for d in data_paths]
        for d in data_paths:
            d_copy = os.path.join(self.BASE_DIR, os.path.basename(d))
            shutil.copy(d, d_copy)
            os.chmod(d_copy, 0o660)

        self.ifgs = sydney_data_setup(datafiles=self.new_data_paths)

        for i in self.ifgs:
            if not i.is_open:
                i.open()
            if not i.nan_converted:
                i.convert_to_nans()

            if not i.mm_converted:
                i.convert_to_mm()
                i.write_modified_phase()

    def tearDown(self):
        for i in self.ifgs:
            i.close()
        shutil.rmtree(self.BASE_DIR)

    def test_orbital_correction_matlab_equality_orbfit_method_2(self):
        pyrate.orbital.remove_orbital_error(self.ifgs, self.params)

        onlyfiles = [f for f in os.listdir(SYD_TEST_MATLAB_ORBITAL_DIR)
            if os.path.isfile(os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, f))
            and f.endswith('.csv') and f.__contains__('_method2_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            matlab_phase_data = np.genfromtxt(os.path.join(
                SYD_TEST_MATLAB_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if os.path.basename(j.data_path).split('_unw.')[0] == \
                        os.path.basename(f).split(
                            '_method2_')[1].split('.')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(matlab_phase_data,
                        j.phase_data, decimal=3)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(matlab_phase_data)),
                                np.sum(np.isnan(j.phase_data)))

        # ensure that we have expected number of matches
        self.assertEqual(count, len(self.ifgs))

    def test_orbital_error_method2_dummy(self):
        """
        does not test anything except that the method is working
        """
        # change to orbital error correction method 2
        self.params[cf.ORBITAL_FIT_METHOD] = 2
        self.params[cf.ORBITAL_FIT_LOOKS_X] = 2
        self.params[cf.ORBITAL_FIT_LOOKS_Y] = 2

        pyrate.orbital.remove_orbital_error(self.ifgs, self.params)

        onlyfiles = [f for f in os.listdir(SYD_TEST_MATLAB_ORBITAL_DIR)
            if os.path.isfile(os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, f))
            and f.endswith('.csv') and f.__contains__('_method2_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            matlab_phase_data = np.genfromtxt(os.path.join(
                SYD_TEST_MATLAB_ORBITAL_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if os.path.basename(j.data_path).split('_unw.')[0] == \
                        os.path.basename(f).split(
                            '_method2_')[1].split('.')[0]:
                    count += 1
                    # # all numbers equal
                    # np.testing.assert_array_almost_equal(matlab_phase_data,
                    #     j.phase_data, decimal=3)
                    #
                    # # number of nans must equal
                    # self.assertEqual(np.sum(np.isnan(matlab_phase_data)),
                    #             np.sum(np.isnan(j.phase_data)))

        # ensure that we have expected number of matches
        self.assertEqual(count, len(self.ifgs))

# TODO: Write tests for various looks and degree combinations
# TODO: write mpi tests

if __name__ == "__main__":
    unittest.main()
