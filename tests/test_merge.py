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
This Python module contains tests for the Merge step of PyRate.
"""
import os
import unittest
from pathlib import Path
from pyrate.merge import create_png_from_tif
from tests.common import TESTDIR



class MergingTest(unittest.TestCase):

    def test_png_creation(self):

        output_folder_path = Path(TESTDIR).joinpath("test_data", "merge")
        create_png_from_tif(output_folder_path)

        # check if color map is created
        output_color_map_path = os.path.join(output_folder_path, "colourmap.txt")
        if not os.path.isfile(output_color_map_path):
            self.assertTrue(False, "Output color map file not found at: " + output_color_map_path)

        # check if png is created
        output_image_path = os.path.join(output_folder_path, "stack_rate.png")
        if not os.path.isfile(output_image_path):
            self.assertTrue(False, "Output png file not found at: " + output_image_path)

        # check if kml is created
        output_kml_path = os.path.join(output_folder_path, "stack_rate.kml")
        if not os.path.isfile(output_kml_path):
            self.assertTrue(False, "Output kml file not found at: " + output_kml_path)

        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
