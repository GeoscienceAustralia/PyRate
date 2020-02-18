# coding: utf-8
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
This Python module contains tests for mpi operations in PyRate.
Tun this module as 'mpirun -n 4 pytest tests/test_mpi.py'
"""
import os
import unittest

from __main__ import conv2tif_handler, prepifg_handler, process_handler, merge_handler


class SystemTest(unittest.TestCase):
    def setUp(self):
        self.root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        self.rows = 1
        self.cols = 1

    def test_roipac_workflow(self):
        input_config_path = os.path.join(self.root_path, "tests", "test_data", "system", "roipac",
                                         "input_parameters.conf")

        conv2tif_handler(input_config_path)
        prepifg_handler(input_config_path)
        process_handler(input_config_path, self.rows, self.cols)
        merge_handler(input_config_path, self.rows, self.cols)

        self.assertTrue(True)

    def test_gamma_workflow(self):
        input_config_path = os.path.join(self.root_path, "tests", "test_data", "system", "gamma",
                                         "input_parameters.conf")
        conv2tif_handler(input_config_path)
        prepifg_handler(input_config_path)
        process_handler(input_config_path, self.rows, self.cols)
        merge_handler(input_config_path, self.rows, self.cols)
        self.assertTrue(True)

    def test_geotiff_workflow(self):
        input_config_path = os.path.join(self.root_path, "tests", "test_data", "system", "geotiff",
                                         "input_parameters.conf")
        prepifg_handler(input_config_path)
        process_handler(input_config_path, self.rows, self.cols)
        merge_handler(input_config_path, self.rows, self.cols)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
