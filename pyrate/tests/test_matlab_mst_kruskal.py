__author__ = 'sudipta'

import unittest
import numpy as np
import os

from pyrate.tests.common import sydney_data_setup
from pyrate.matlab_mst_kruskal import get_nml
from pyrate.matlab_mst_kruskal import sort_list, get_sub_structure
from pyrate.matlab_mst_kruskal import matlab_mst_kruskal
from pyrate.matlab_mst_kruskal import matlab_mst, matlab_mst_boolean_array
from pyrate.tests.common import sydney_data_setup_ifg_file_list
from pyrate.matlab_mst_kruskal import IfgListPyRate as IfgList
from pyrate.matlab_mst_kruskal import calculate_connect_and_ntrees


class IfgListTest(unittest.TestCase):

    def setUp(self):

        self.matlab_n = [1,  2,  3,  3,  4,  4,  4,  5, 5,  6,  6,  7,  7,  8,
                         9, 10, 11,  3,  5,  7,  9,  5,  6,  8, 11, 12, 8, 13,
                         9, 10, 13, 10, 11, 12]
        self.matlab_masnum = [1, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6,
                              7, 7, 8, 9, 10, 11]
        self.matlab_slvnum = [3, 5, 7, 9, 5, 6, 8, 11, 12,  8, 13, 9,
                              10, 13, 10, 11, 12]

        ifg_instance = IfgList(sydney_data_setup_ifg_file_list())
        self.ifg_list, self.epoch_list = get_nml(ifg_instance)

    def test_matlab_n(self):
        # add 1 to ifg_list.n as matlab indices start from 1.
        np.testing.assert_array_equal(self.matlab_n, self.ifg_list.n + 1)

    def test_matlab_masnum(self):
        # add 1 to ifg_list.masnum as matlab indices start from 1.
        np.testing.assert_array_equal(self.matlab_masnum,
                                      self.ifg_list.master_num + 1)

    def test_matlab_slvnum(self):
        # add 1 to ifg_list.slvnum as matlab indices start from 1.
        np.testing.assert_array_equal(self.matlab_slvnum,
                                      self.ifg_list.slave_num + 1)


class MatlabMstKruskalTest(unittest.TestCase):

    def setUp(self):
        ifg_instance = IfgList(sydney_data_setup_ifg_file_list())
        self.ifg_list, _ = get_nml(ifg_instance)
        self.sorted_list = sort_list(self.ifg_list.id, self.ifg_list.master_num,
                                self.ifg_list.slave_num, self.ifg_list.nan_frac)
        self.matlab_sorted_list_zero_nan_frac = [   (1, 1, 3, 0.0),
                                                    (2, 2, 5, 0.0),
                                                    (3, 3, 7, 0.0),
                                                    (4, 3, 9, 0.0),
                                                    (5, 4, 5, 0.0),
                                                    (6, 4, 6, 0.0),
                                                    (7, 4, 8, 0.0),
                                                    (8, 5, 11, 0.0),
                                                    (9, 5, 12, 0.0),
                                                    (10, 6, 8, 0.0),
                                                    (11, 6, 13, 0.0),
                                                    (12, 7, 9, 0.0),
                                                    (13, 7, 10, 0.0),
                                                    (14, 8, 13, 0.0),
                                                    (15, 9, 10, 0.0),
                                                    (16, 10, 11, 0.0),
                                                    (17, 11, 12, 0.0)]
        self.ifg_list_mst_matlab = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 16]

    def test_sorted_list_equal_matlab(self):
        for s, m in zip(self.sorted_list,
                        self.matlab_sorted_list_zero_nan_frac):
            self.assertEquals((s[0]+1, s[1]+1, s[2]+1, s[3]), m)

    def test_mst_kruskal_matlab(self):
        edges = get_sub_structure(self.ifg_list,
                                  np.zeros(len(self.ifg_list.id), dtype=bool))
        ifg_list_mst = matlab_mst_kruskal(edges)
        ifg_list_mst = [i + 1 for i in ifg_list_mst]  # add 1 to each index
        self.assertSequenceEqual(ifg_list_mst, self.ifg_list_mst_matlab)


class MSTKruskalCalcConnectAndNTrees(unittest.TestCase):
    def test_calculate_connect_and_ntrees(self):
        connect_start = np.ones(shape=(10, 10), dtype=bool)
        connect_start[1, :-1] = False
        _, connect, ntrees = calculate_connect_and_ntrees(connect_start, [])
        self.assertEqual(ntrees, 9)
        np.testing.assert_array_equal(connect,
                                      np.ones(shape=(9, 10), dtype=bool))
        connect_start[2:5, :-1] = False

        _, connect, ntrees = calculate_connect_and_ntrees(connect_start, [])
        self.assertEqual(ntrees, 6)
        np.testing.assert_array_equal(connect,
                                      np.ones(shape=(6, 10), dtype=bool))

    def test_calculate_connect_and_ntrees_raise_1(self):
        connect_start = np.eye(10, dtype=bool)
        connect_start[1, :] = False
        with self.assertRaises(ValueError):
            calculate_connect_and_ntrees(connect_start, [])

    def test_calculate_connect_and_ntrees_raise_2(self):
        connect_start = np.eye(10, dtype=bool)
        with self.assertRaises(ValueError):
            calculate_connect_and_ntrees(connect_start, [])


class MSTKruskalConnectAndTresSydneyData(unittest.TestCase):

    def test_calculate_connect_and_ntrees_sydney_data(self):
        ifg_instance = IfgList(sydney_data_setup_ifg_file_list())
        ifg_list, _ = get_nml(ifg_instance)
        edges = get_sub_structure(ifg_list, np.zeros(len(ifg_list.id), dtype=bool))
        mst, connected, ntrees = matlab_mst_kruskal(edges,
            ntrees=True
            )

        self.assertTrue(connected[0].all())
        self.assertEqual(ntrees, 1)
        self.assertEqual(len(connected[0]), len(mst) + 1)

    def test_assert_is_not_tree(self):
        non_overlapping = [1, 2, 5, 6, 12, 13, 14, 15, 16, 17]
        sydney_files = sydney_data_setup_ifg_file_list()
        datafiles = [f for i, f in enumerate(sydney_files)
                     if i+1 in non_overlapping]
        non_overlapping_ifg_isntance = IfgList(datafiles)
        ifg_list, _ = get_nml(non_overlapping_ifg_isntance)
        edges = get_sub_structure(ifg_list,
                                  np.zeros(len(ifg_list.id), dtype=bool))
        mst, connected, ntrees = matlab_mst_kruskal(edges, ntrees=True)
        self.assertEqual(connected.shape[0], 4)
        self.assertEqual(ntrees, 4)

    def test_assert_is_tree(self):
        overlapping = [1, 2, 3, 4, 6, 7, 10, 11, 16, 17]
        sydney_files = sydney_data_setup_ifg_file_list()
        datafiles = [f for i, f in enumerate(sydney_files)
                     if i+1 in overlapping]

        overlapping_ifg_isntance = IfgList(datafiles)

        ifg_list, _ = get_nml(overlapping_ifg_isntance)
        edges = get_sub_structure(ifg_list,
                                  np.zeros(len(ifg_list.id), dtype=bool))
        _, connected, ntrees = matlab_mst_kruskal(edges, ntrees=True)
        self.assertEqual(ntrees, 4)
        self.assertEqual(connected.shape[0], 4)


class MatlabMSTTests(unittest.TestCase):
    """
    Tests to ensure matlab and python mst outputs are the same.
    """
    def setUp(self):
        self.ifgs = sydney_data_setup()
        self.ifg_file_list = sydney_data_setup_ifg_file_list()

        self.matlab_mst_list = ['geo_060619-061002.tif',
                                'geo_060828-061211.tif',
                                'geo_061002-070219.tif',
                                'geo_061002-070430.tif',
                                'geo_061106-061211.tif',
                                'geo_061106-070115.tif',
                                'geo_061106-070326.tif',
                                'geo_061211-070709.tif',
                                'geo_061211-070813.tif',
                                'geo_070115-070917.tif',
                                'geo_070219-070604.tif',
                                'geo_070604-070709.tif']

        # reorder ifgs as per the matlab list
        self.ifgs = sydney_data_setup()

    def test_matlab_mst_kruskal(self):
        """
        test that the matlab and python mst algos outputs are the same
        """
        ifg_instance = IfgList(datafiles=self.ifg_file_list)
        ifg_list, _ = get_nml(ifg_instance)
        edges = get_sub_structure(ifg_list,
                                  np.zeros(len(ifg_list.id), dtype=bool))
        ifg_list_mst_id = matlab_mst_kruskal(edges)

        self.assertEquals(len(self.matlab_mst_list),
                          len(ifg_list_mst_id))

        for i in ifg_list_mst_id:
            self.assertIn(os.path.split(ifg_list.nml[i])[-1],
                          self.matlab_mst_list)

    def test_matlab_make_mstmat(self):
        """
        tests equality of boolean mst arrays of both python and matlab.
        """
        ifg_instance = IfgList(datafiles=self.ifg_file_list)
        ifg_list, _ = get_nml(ifg_instance)
        mst_mat = matlab_mst(ifg_list, p_threshold=1)

        # path to csv folders from matlab output
        from pyrate.tests.common import SYD_TEST_MATLAB_MST_DIR

        onlyfiles = [f for f in os.listdir(SYD_TEST_MATLAB_MST_DIR)
                if os.path.isfile(os.path.join(SYD_TEST_MATLAB_MST_DIR, f))]

        for i, f in enumerate(onlyfiles):
            mst_f = np.genfromtxt(os.path.join(SYD_TEST_MATLAB_MST_DIR, f),
                                  delimiter=',')
            for k, j in enumerate(self.ifg_file_list):
                if f.split('matlab_')[-1].split('.')[0] == \
                        os.path.split(j)[-1].split('.')[0]:
                    np.testing.assert_array_equal(mst_f, mst_mat[k, :, :])

    def test_matlab_make_mstmat_boolean_array(self):
        """
        tests equality of boolean mst arrays of both python and matlab.
        """
        ifg_instance = IfgList(datafiles=self.ifg_file_list)
        ifg_list, _ = get_nml(ifg_instance)
        mst_mat = matlab_mst_boolean_array(ifg_list, p_threshold=1)

        # path to csv folders from matlab output
        from pyrate.tests.common import SYD_TEST_MATLAB_MST_DIR

        onlyfiles = [f for f in os.listdir(SYD_TEST_MATLAB_MST_DIR)
                if os.path.isfile(os.path.join(SYD_TEST_MATLAB_MST_DIR, f))]

        for i, f in enumerate(onlyfiles):
            mst_f = np.genfromtxt(os.path.join(SYD_TEST_MATLAB_MST_DIR, f),
                                  delimiter=',')
            for k, j in enumerate(self.ifg_file_list):
                if f.split('matlab_')[-1].split('.')[0] == \
                        os.path.split(j)[-1].split('.')[0]:
                    np.testing.assert_array_equal(mst_f, mst_mat[k, :, :])

    def test_mas_mat_vs_mst_mat_generator(self):
        ifg_instance = IfgList(datafiles=self.ifg_file_list)
        ifg_list, _ = get_nml(ifg_instance, nan_conversion=True)
        mst_mat1 = matlab_mst(ifg_list)
        mst_mat2 = matlab_mst_boolean_array(ifg_list)

        np.testing.assert_array_equal(mst_mat2, mst_mat1)


if __name__ == '__main__':
    unittest.main()
