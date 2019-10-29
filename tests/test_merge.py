# coding: utf-8
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
This Python module contains tests for mpi operations in PyRate.
Tun this module as 'mpirun -n 4 pytest tests/test_mpi.py'
"""
import glob
import shutil
import numpy as np
import pytest
import os
import tempfile
import random
import string

import pyrate.core.orbital
import pyrate.core.shared
import tests.common
from pyrate import (
    process, prepifg, merge, conv2tif)
from tests.common import (small_data_setup, reconstruct_mst,
    reconstruct_linrate, SML_TEST_DEM_HDR_GAMMA, pre_prepare_ifgs)
from tests import common
from tests.test_covariance import legacy_maxvar
from pyrate.core import algorithm, ref_phs_est as rpe, mpiops, config as cf, covariance, refpixel
from pyrate.merge import create_png_from_tif
import unittest


class MergingTest(unittest.TestCase):

    def test_png_creation(self):

        output_folder_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"tests", "test_data", "merge")
        create_png_from_tif(output_folder_path)

        # check if color map is created
        output_color_map_path = os.path.join(output_folder_path, "colormap.txt")
        if not os.path.isfile(output_color_map_path):
            self.assertTrue(False, "Output color map file not found at: " + output_color_map_path)

        # check if png is created
        output_image_path = os.path.join(output_folder_path, "linrate.png")
        if not os.path.isfile(output_image_path):
            self.assertTrue(False, "Output png file not found at: " + output_image_path)

        # check if kml is created
        output_kml_path = os.path.join(output_folder_path, "linrate.kml")
        if not os.path.isfile(output_kml_path):
            self.assertTrue(False, "Output kml file not found at: " + output_kml_path)

        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
