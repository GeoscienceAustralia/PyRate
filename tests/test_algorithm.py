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

from . import common

from core.algorithm import get_all_epochs, get_epochs, master_slave_ids
from core.config import parse_namelist
from core.shared import Ifg, convert_radians_to_mm


class AlgorithmTests(TestCase):
    """Misc unittests for functions in the algorithm module."""

    @staticmethod
    def test_phase_conversion():
        """ """
        # ROIPAC interferograms in units of radians, verify conversion to mm
        xs, ys = 5, 7
        data = (np.arange(xs * ys) - 1.7) * 0.1  # fake a range of values
        data = np.where(data == 0, np.nan, data)
        wavelen = 0.0562356424
        exp = (data * wavelen * 1000) / (4 * pi)
        act = convert_radians_to_mm(data, wavelen)
        assert_allclose(exp, act)


class EpochsTests(TestCase):
    """Unittests for the EpochList class."""

    def test_get_epochs(self):
        """ """
        def str2date(s):
            """

            Args:
              s: 

            Returns:

            """
            segs = s[:4], s[4:6], s[6:]  # year, month, day
            return date(*[int(sg) for sg in segs])

        raw_date = [
            "20060619",
            "20060828",
            "20061002",
            "20061106",
            "20061211",
            "20070115",
            "20070219",
            "20070326",
            "20070430",
            "20070604",
            "20070709",
            "20070813",
            "20070917",
        ]

        exp_dates = [str2date(d) for d in raw_date]
        exp_repeat = [1, 1, 3, 3, 4, 3, 3, 3, 3, 3, 3, 2, 2]
        exp_spans = [0, 0.1916, 0.2875, 0.3833, 0.4791, 0.5749, 0.6708, 0.7666, 0.8624, 0.9582, 1.0541, 1.1499, 1.2457]

        ifms = join(common.SML_TEST_TIF, "ifms_17")
        ifgs = [Ifg(join(common.SML_TEST_TIF, p)) for p in parse_namelist(ifms)]
        for i in ifgs:
            i.open()

        epochs = get_epochs(ifgs)[0]

        self.assertTrue((exp_dates == epochs.dates).all())
        self.assertTrue((exp_repeat == epochs.repeat).all())
        assert_array_almost_equal(exp_spans, epochs.spans, decimal=4)

    def test_get_all_epochs(self):
        """ """
        # test function to extract all dates from sequence of ifgs
        ifgs = common.small5_mock_ifgs()
        for i in ifgs:
            i.nodata_value = 0
        dates = [date(2006, 8, 28), date(2006, 11, 6), date(2006, 12, 11), date(2007, 1, 15), date(2007, 3, 26), date(2007, 9, 17)]

        self.assertEqual(dates, sorted(set(get_all_epochs(ifgs))))

    def test_get_epoch_count(self):
        """ """
        self.assertEqual(6, len(set(get_all_epochs(common.small5_mock_ifgs()))))

    def test_master_slave_ids(self):
        """ """
        d0 = date(2006, 6, 19)
        d1 = date(2006, 8, 28)
        d2 = date(2006, 10, 2)
        d3 = date(2006, 11, 6)
        exp = {d0: 0, d1: 1, d2: 2, d3: 3}

        # test unordered and with duplicates
        self.assertEqual(exp, master_slave_ids([d3, d0, d2, d1]))
        self.assertEqual(exp, master_slave_ids([d3, d0, d2, d1, d3, d0]))


if __name__ == "__main__":
    unittest.main()
