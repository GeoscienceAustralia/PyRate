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
'''
This module contains tests for the mst.py PyRate module.
'''

import glob
import os
import shutil
import subprocess
import tempfile
import unittest
from itertools import product
from numpy import empty, array, nan, isnan, sum as nsum

import numpy as np
from tests.common import MockIfg, sydney5_mock_ifgs, sydney_data_setup

from pyrate import algorithm
from pyrate import config as cf
from pyrate import mst
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate.shared import IfgPart, Tile, create_tiles
from tests import common


class MSTTests(unittest.TestCase):
    '''Basic verification of minimum spanning tree (MST) functionality.'''

    def setUp(self):
        self.ifgs = sydney_data_setup()

    def test_mst_matrix_as_array(self):
        # Verifies MST matrix func returns array with dict/trees in each cell
        for i in self.ifgs[3:]:
            i.phase_data[0, 1] = 0  # partial stack of NODATA to one cell

        for i in self.ifgs:
            i.convert_to_nans() # zeros to NaN/NODATA

        epochs = algorithm.get_epochs(self.ifgs)[0]
        res = mst.mst_matrix_as_array(self.ifgs)
        ys, xs = res.shape

        for y, x in product(range(ys), range(xs)):
            r = res[y, x]
            num_nodes = len(r)
            self.assertTrue(num_nodes < len(epochs.dates))

            stack = array([i.phase_data[y, x] for i in self.ifgs])  # 17 ifg stack
            self.assertTrue(0 == nsum(stack == 0))  # all 0s should be converted
            nc = nsum(isnan(stack))
            exp_count = len(epochs.dates) - 1

            if nc == 0:
                self.assertEqual(num_nodes, exp_count)
            elif nc > 5:
                # rough test: too many nans must reduce the total tree size
                self.assertTrue(num_nodes <= (17-nc))

    def test_mst_matrix_as_ifgs(self):
        # ensure only ifgs are returned, not individual MST graphs
        ifgs = sydney5_mock_ifgs()
        nifgs = len(ifgs)
        ys, xs = ifgs[0].shape
        result = mst.mst_matrix_ifgs_only(ifgs)

        for coord in product(range(ys), range(xs)):
            stack = (i.phase_data[coord] for i in self.ifgs)
            nc = nsum([isnan(n) for n in stack])
            self.assertTrue(len(result[coord]) <= (nifgs - nc))

            # HACK: type testing here is a bit grubby
            self.assertTrue(all([isinstance(i, MockIfg) for i in ifgs]))

    def test_partial_nan_pixel_stack(self):
        # Ensure a limited # of coherent cells results in a smaller MST tree
        num_coherent = 3

        def assert_equal():
            res = mst.mst_matrix_as_array(mock_ifgs)
            self.assertEqual(len(res[0,0]), num_coherent)

        mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
        for m in mock_ifgs[num_coherent:]:
            m.phase_data[:] = nan
        assert_equal()

        # fill in more nans leaving only one ifg
        for m in mock_ifgs[1:num_coherent]:
            m.phase_data[:] = nan
        num_coherent = 1
        assert_equal()

    def test_all_nan_pixel_stack(self):
        # ensure full stack of NaNs in an MST pixel classifies to NaN
        mock_ifgs = [MockIfg(i, 1, 1) for i in self.ifgs]
        for m in mock_ifgs:
            m.phase_data[:] = nan

        res = mst.mst_matrix_as_array(mock_ifgs)
        exp = empty((1,1), dtype=object)
        exp[:] = nan

        shape = (mock_ifgs[0].nrows, mock_ifgs[0].ncols)
        self.assertTrue(res.shape == shape)
        self.assertEqual(exp, res)


class DefaultMSTTests(unittest.TestCase):

    def test_default_mst(self):
        # default MST from full set of Ifgs shouldn't drop any nodes
        ifgs = sydney5_mock_ifgs()
        dates = [(i.master, i.slave) for i in ifgs]

        res = mst.mst_from_ifgs(ifgs)[0]
        num_edges = len(res)
        self.assertEqual(num_edges, len(ifgs))

        # test edges, note node order can be reversed
        for edge in res:
            self.assertTrue(edge in dates or (edge[1], edge[0]) in dates)

        # check all nodes exist in this default tree
        mst_dates = set(res)
        mst_dates = list(sum(mst_dates, ()))
        for i in ifgs:
            for node in (i.master, i.slave):
                self.assertIn(node, mst_dates)


class NetworkxMSTTreeCheck(unittest.TestCase):
    def setUp(self):
        self.ifgs = sydney_data_setup()

    def test_assert_is_not_tree(self):
        non_overlapping = [1, 2, 5, 6, 12, 13, 14, 15, 16, 17]
        ifgs_non_overlapping = [ifg for i, ifg in enumerate(self.ifgs)
                                if i+1 in non_overlapping]
        edges, is_tree, ntrees, _ = mst.mst_from_ifgs(ifgs_non_overlapping)
        self.assertFalse(is_tree)
        self.assertEqual(4, ntrees)

    def test_sydney_data_tree(self):
        self.assertTrue(mst.mst_from_ifgs(self.ifgs)[1])

    def test_assert_is_tree(self):
        overlapping = [1, 2, 3, 4, 6, 7, 10, 11, 16, 17]

        ifgs_overlapping = [ifg for i, ifg in enumerate(self.ifgs)
                                if (i+1 in overlapping)]
        edges, is_tree, ntrees, _ = mst.mst_from_ifgs(ifgs_overlapping)
        self.assertFalse(is_tree)
        self.assertEqual(4, ntrees)

    def test_assert_two_trees_overlapping(self):
        overlapping = [3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17]

        ifgs_overlapping = [ifg for i, ifg in enumerate(self.ifgs)
                                if (i+1 in overlapping)]
        edges, is_tree, ntrees, _ = mst.mst_from_ifgs(ifgs_overlapping)
        self.assertFalse(is_tree)
        self.assertEqual(2, ntrees)

    def test_assert_two_trees_non_overlapping(self):
        non_overlapping = [2, 5, 6, 12, 13, 15]
        ifgs_non_overlapping = [ifg for i, ifg in enumerate(self.ifgs)
                                if i+1 in non_overlapping]
        edges, is_tree, ntrees, _ = mst.mst_from_ifgs(ifgs_non_overlapping)
        self.assertFalse(is_tree)
        self.assertEqual(2, ntrees)


class IfgPartTest(unittest.TestCase):

    def setUp(self):
        self.ifgs = sydney_data_setup()
        self.params = cf.get_config_params(
            os.path.join(common.SYD_TEST_DIR, 'pyrate_system_test.conf'))

    def test_ifg_part_shape_and_slice(self):
        r_start = 0
        r_end = 10
        for i in self.ifgs:
            tile = Tile(0, top_left=(r_start, 0), bottom_right=(r_end, i.ncols))
            ifg_part = IfgPart(i.data_path, tile)
            self.assertEqual(ifg_part.phase_data.shape,
                             (r_end-r_start, i.phase_data.shape[1]))
            np.testing.assert_array_equal(ifg_part.phase_data,
                                          i.phase_data[r_start:r_end, :])

    def test_mst_multiprocessing_serial(self):
        self.params[cf.PARALLEL] = False
        original_mst = mst.mst_boolean_array(self.ifgs)
        parallel_mst = mst.mst_parallel(self.ifgs, params=self.params)
        np.testing.assert_array_equal(original_mst, parallel_mst)

    def test_mst_multiprocessing(self):
        self.params[cf.PARALLEL] = True
        original_mst = mst.mst_boolean_array(self.ifgs)
        parallel_mst = mst.mst_parallel(self.ifgs, params=self.params)
        np.testing.assert_array_equal(original_mst, parallel_mst)


if __name__ == "__main__":
    unittest.main()
