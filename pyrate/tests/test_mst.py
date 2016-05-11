'''
Tests for Minimum Spanning Tree (MST) functionality in PyRate.

.. codeauthor:: Ben Davies
'''

import unittest
import os
from itertools import product
from numpy import empty, array, nan, isnan, sum as nsum
import numpy as np
import tempfile
import shutil
import sys
import glob
import subprocess

from pyrate import mst
from pyrate import algorithm
from pyrate.tests.common import MockIfg, sydney5_mock_ifgs, sydney_data_setup
from pyrate.shared import Ifg, IfgPart
from pyrate.mst import mst_parallel
from pyrate import config as cf
from pyrate.tests import common
from pyrate.scripts import run_pyrate, run_prepifg


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

        epochs = algorithm.get_epochs(self.ifgs)
        res = mst.mst_matrix_as_array(self.ifgs)
        ys, xs = res.shape

        for y, x in product(xrange(ys), xrange(xs)):
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

        for coord in product(xrange(ys), xrange(xs)):
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
            ifg_part = IfgPart(i.data_path,
                               r_start=r_start, r_end=r_end,
                               c_start=0, c_end=i.ncols)
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


class MPITest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        params = cf.get_config_params(cls.test_conf)
        params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        params[cf.PROCESSOR] = 1  # gamma
        params[cf.IFG_FILE_LIST] = os.path.join(
            common.SYD_TEST_GAMMA, 'ifms_17')
        params[cf.OUT_DIR] = cls.tif_dir
        params[cf.PARALLEL] = 1

        xlks, ylks, crop = run_pyrate.transform_params(params)

        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        # dest_paths are tifs that have been geotif converted and multilooked
        dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, params, xlks)

        run_prepifg.gamma_prepifg(base_unw_paths, params)
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

        cls.log_file = os.path.join(cls.tif_dir, 'mst_mpi.log')
        # Calc mst using MPI
        cls.conf_file = tempfile.mktemp(suffix='.conf', dir=cls.tif_dir)
        cf.write_config_file(params, cls.conf_file)
        str = 'mpirun -np 2 python pyrate/nci/run_pyrate_pypar.py ' + cls.conf_file
        cmd = str.split()

        subprocess.check_call(cmd)

        mst_mat_binary_file = os.path.join(params[cf.OUT_DIR], 'mst_mat.npy')
        cls.mst = np.load(mst_mat_binary_file)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir)

    def test_mpi_mst_sigmle_processor(self):
        mlooked_ifgs = glob.glob(os.path.join(self.tif_dir, '*cr.tif'))
        self.assertEqual(len(mlooked_ifgs), 17)
        original_mst = mst.mst_boolean_array(self.ifgs)
        np.testing.assert_array_equal(original_mst, self.mst)

    def test_mst_log_written(self):
        log_file = glob.glob(os.path.join(self.tif_dir, '*.log'))[0]
        self.assertTrue(os.path.exists(log_file))

if __name__ == "__main__":
    unittest.main()
