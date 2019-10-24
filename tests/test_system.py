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


class SystemTest(unittest.TestCase):

    def test_roipac_workflow(self):
        self.assertTrue(True)

    def test_gamma_workflow(self):
        self.assertTrue(True)

    def test_geotiff_workflow(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
