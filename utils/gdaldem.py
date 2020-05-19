# This Python module is part of the PyRate software package
#
# Copyright 2020 Geoscience Australia
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
Based on: http://gis.stackexchange.com/questions/130199/changing-color-of-raster-images-based-on-their-data-values-gdal
Utility script to convert grayscale geotiff into colour.
================================================================================
Usage:
-------------
Custom color:
python gdaldem.py input.tif color.txt output_color.tif
Auto color:
python gdaldem.py input.tif auto output_color.tif

Examples:
-------------
Custom color:
python utils/gdaldem.py tests/test_data/small_test/dem/roipac_test_trimmed.tif utils/tiff_colour_map.txt out.tif

Auto color:
python utils/gdaldem.py tests/test_data/small_test/dem/roipac_test_trimmed.tif auto out.tif

"""
import subprocess
import sys
import os
import tempfile
import numpy as np
from pyrate.core.shared import DEM


def main(input_file, color_file, output_file):
    cmd = "gdaldem color-relief " + input_file \
          + ' ' + color_file + ' ' + output_file
    subprocess.check_call(cmd, shell=True)


def gen_color_file(input_file):
    fp, temp_file = tempfile.mkstemp(suffix='.txt')

    dem = DEM(input_file)
    dem.open()
    phase_data = dem.height_band.ReadAsArray()

    max_ph = np.nanmax(phase_data)
    min_ph = np.nanmin(phase_data)
    range_ph = max_ph-min_ph
    colors = ['black', 'blue', 'yellow', 'orange', 'red', 'white']
    with open(temp_file, 'w') as f:
        for i, c in enumerate(colors[:-1]):
            f.write(str(int(min_ph + (i + 1)*range_ph/len(colors))) +
                    ' ' + c + '\n')
        f.write(str(int(max_ph - range_ph/len(colors))) +
                ' ' + colors[-1] + '\n')
    os.close(fp)
    return temp_file


if __name__ == '__main__':
    input_file = sys.argv[1]
    color_file = sys.argv[2]
    output_file = sys.argv[3]
    if color_file == 'auto':
        print('\nauto generating color file')
        color_file = gen_color_file(input_file)
    with open(color_file, 'r') as f:
        print('\ncolor file contents')
        for l in f.readlines():
            print(l)
    main(input_file, color_file, output_file)
