__author__ = 'sudipta'

import unittest
import numpy as np
from pyrate.algorithm import get_epochs
from pyrate.tests.common import sydney_data_setup
from pyrate.matlab_mst_kruskal import get_nml
from pyrate.matlab_mst_kruskal import sort_list
from pyrate.matlab_mst_kruskal import matlab_mst_kruskal


class IfgListTest(unittest.TestCase):

    def setUp(self):
        self.ifgs = self.ifg = sydney_data_setup()
        self.matlab_n = [1,  2,  3,  3,  4,  4,  4,  5, 5,  6,  6,  7,  7,  8,
                         9, 10, 11,  3,  5,  7,  9,  5,  6,  8, 11, 12, 8, 13,
                         9, 10, 13, 10, 11, 12]
        self.matlab_masnum = [1, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6,
                              7, 7, 8, 9, 10, 11]
        self.matlab_slvnum = [3, 5, 7, 9, 5, 6, 8, 11, 12,  8, 13, 9,
                              10, 13, 10, 11, 12]

        ifgs = sydney_data_setup()
        self.ifg_list, self.epoch_list = get_nml(ifgs)

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
        self.ifg_list, _ = get_nml()
        self.sorted_list = sort_list(self.ifg_list)
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
        ifg_list_mst = matlab_mst_kruskal(self.ifg_list)
        ifg_list_mst = [i + 1 for i in ifg_list_mst]  # add 1 to each index
        self.assertSequenceEqual(ifg_list_mst, self.ifg_list_mst_matlab)

if __name__ == '__main__':
    unittest.main()
