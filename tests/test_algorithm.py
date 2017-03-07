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
This Python module contains tests for the algorithm.py PyRate module.
"""
import unittest
from datetime import date
from math import pi, cos, sin, radians
from numpy import array, reshape, squeeze
from os.path import join
from unittest import TestCase

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_allclose

from pyrate.algorithm import (least_squares_covariance,
                              is_square,
                              unit_vector,
                              ifg_date_lookup,
                              get_all_epochs,
                              get_epochs,
                              master_slave_ids,
                              )

from pyrate.config import parse_namelist
from pyrate.shared import Ifg, convert_radians_to_mm
from tests.common import sydney5_mock_ifgs, SYD_TEST_TIF


class LeastSquaresTests(TestCase):
    """
    Unit tests for the PyRate least_squares_covariance() implementation.
    """

    @staticmethod
    def test_least_squares_covariance():
        b = array([[13, 7.2, 5.7]]).T
        A = array([[1, 0.4, 0.3], [1, 1, 1]]).T
        v = array([[1, 1, 1]]).T
        r = least_squares_covariance(A, b, v)
        exp = [10.1628, 2.8744]
        assert_array_almost_equal(r.T.squeeze(), exp, decimal=4)

    def test_least_squares_covariance_overdetermined(self):
        # must be overdetermined, ie. more observations than params
        b = array([[10]]).T
        A = array([[1]]).T
        v = array([[1]]).T
        self.assertRaises(ValueError,
                          least_squares_covariance, A, b, v)

        # try non transposed style
        b = array([[10]])
        A = array([[1]])
        v = array([[1]])
        self.assertRaises(ValueError,
                          least_squares_covariance, A, b, v)


class AlgorithmTests(TestCase):
    """
    Misc unittests for functions in the algorithm module.
    """

    def test_is_square(self):
        self.assertTrue(is_square(np.empty((2, 2))))

    def test_is_not_square(self):
        for shape in [(3, 2), (2, 3)]:
            self.assertFalse(is_square(np.empty(shape)))

    @staticmethod
    def test_phase_conversion():
        # ROIPAC interferograms in units of radians, verify conversion to mm
        xs, ys = 5, 7
        data = (np.arange(xs * ys) - 1.7) * 0.1 # fake a range of values
        data = np.where(data == 0, np.nan, data)
        wavelen = 0.0562356424
        exp = (data * wavelen * 1000) / (4 * pi)
        act = convert_radians_to_mm(data, wavelen)
        assert_allclose(exp, act)

    def test_unit_vector(self):
        # last values here simulate a descending pass
        incidence = [radians(x) for x in (34.3, 39.3, 29.3, 34.3)]
        azimuth = [radians(x) for x in (77.8, 77.9, 80.0, 282.2)]

        vert, ns, ew = [], [], []
        for i, a in zip(incidence, azimuth):
            vert.append(cos(i))
            ns.append(sin(i) * cos(a))
            ew.append(sin(i) * sin(a))

        sh = 4
        unitv = [array(ew), array(ns), array(vert)]
        unitv = [a.reshape(sh) for a in unitv]

        # NB: assumes radian inputs
        act = unit_vector(reshape(incidence, sh),
                                    reshape(azimuth, sh))
        for a, e in zip(act, unitv):
            assert_array_almost_equal(squeeze(a), e)

        # check unit vec components have correct signs
        E, N, V = act
        # test E/W component of ascending is +ve
        self.assertTrue((E[:-2]).all() > 0)
        self.assertTrue(E[-1] < 0)  # test E/W component of descending is -ve
        self.assertTrue((N > 0).all())  # ensure all north values are positive

        # check unit vec components have correct magnitudes
        self.assertTrue((abs(V) > abs(E)).all())
        self.assertTrue((abs(V) > abs(N)).all())
        self.assertTrue((abs(E) > abs(N)).all())


class DateLookupTests(TestCase):
    """
    Tests for the algorithm.ifg_date_lookup() function.
    """

    def setUp(self):
        self.ifgs = sydney5_mock_ifgs()

    def test_ifg_date_lookup(self):
        # check reverse lookup of ifg given a master and slave date tuple
        date_pair = (date(2006, 8, 28), date(2006, 12, 11))
        i = ifg_date_lookup(self.ifgs, date_pair)
        self.assertEqual(self.ifgs[0], i)

        # test with reversed date tuple, should reorder it according to age
        date_pair = (date(2006, 12, 11), date(2006, 11, 6))
        i = ifg_date_lookup(self.ifgs, date_pair)
        self.assertEqual(self.ifgs[1], i)

    def test_ifg_date_lookup_failure(self):
        # error when lookup cannot find an ifg given a date pair
        dates = (date(2006, 12, 11), date(2007, 3, 26))
        self.assertRaises(ValueError,
                          ifg_date_lookup, self.ifgs, dates)

    def test_date_lookup_bad_inputs(self):
        # test some bad inputs to date lookup
        inputs = [(None, None), (1, 10), (34.56, 345.93),
                  (date(2007, 3, 26), ""), (date(2007, 3, 26), None)]

        for d in inputs:
            self.assertRaises(ValueError,
                              ifg_date_lookup, self.ifgs, d)


# TODO: InitialModelTests
#class InitialModelTests(unittest.TestCase):

#    def test_initial_model(self):
        # 1. fake an RSC file with coords
        # 2. fake a ones(shape)  # could also make a ramp etc
        # data is single band of DISPLACEMENT
        #raise NotImplementedError


class EpochsTests(TestCase):
    """
    Unittests for the EpochList class.
    """

    def test_get_epochs(self):
        def str2date(s):
            segs = s[:4], s[4:6], s[6:] # year, month, day
            return date(*[int(sg) for sg in segs])

        raw_date = ['20060619', '20060828', '20061002', '20061106', '20061211',
                    '20070115', '20070219', '20070326', '20070430', '20070604',
                    '20070709', '20070813', '20070917']

        exp_dates = [str2date(d) for d in raw_date]
        exp_repeat = [1, 1, 3, 3, 4, 3, 3, 3, 3, 3, 3, 2, 2]
        exp_spans = [0, 0.1916, 0.2875, 0.3833, 0.4791, 0.5749, 0.6708, 0.7666,
                            0.8624, 0.9582, 1.0541, 1.1499, 1.2457]

        # test against Hua's results
        ifms = join(SYD_TEST_TIF, "ifms_17")
        ifgs = [Ifg(join(SYD_TEST_TIF, p)) for p in parse_namelist(ifms)]
        for i in ifgs:
            i.open()

        epochs = get_epochs(ifgs)[0]

        self.assertTrue((exp_dates == epochs.dates).all())
        self.assertTrue((exp_repeat == epochs.repeat).all())
        assert_array_almost_equal(exp_spans, epochs.spans, decimal=4)

    def test_get_all_epochs(self):
        # test function to extract all dates from sequence of ifgs
        ifgs = sydney5_mock_ifgs()
        for i in ifgs:
            i.nodata_value = 0
        dates = [date(2006, 8, 28), date(2006, 11, 6), date(2006, 12, 11),
                 date(2007, 1, 15), date(2007, 3, 26), date(2007, 9, 17)]

        self.assertEqual(dates, sorted(set(get_all_epochs(ifgs))))

    def test_get_epoch_count(self):
        self.assertEqual(6, len(set(get_all_epochs(sydney5_mock_ifgs()))))

    def test_master_slave_ids(self):
        d0 = date(2006, 6, 19)
        d1 = date(2006, 8, 28)
        d2 = date(2006, 10, 2)
        d3 = date(2006, 11, 6)
        exp = {d0: 0, d1: 1, d2: 2, d3: 3}

        # test unordered and with duplicates
        self.assertEqual(exp, master_slave_ids([d3, d0, d2, d1]))
        self.assertEqual(exp,
                         master_slave_ids([d3, d0, d2, d1, d3, d0]))

if __name__ == "__main__":
    unittest.main()
