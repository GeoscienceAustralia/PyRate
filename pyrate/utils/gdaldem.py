"""
Based on: http://gis.stackexchange.com/questions/130199/changing-color-of-raster-images-based-on-their-data-values-gdal
Utility script to convert gray tif into color.
================================================================================
Usage:
-------------
Custom color:
python gdaldem.py input_tif.tif color.txt output_color.tif
Auto color:
python gdaldem.py input_tif.tif auto output_color.tif

Examples:
-------------
Custom color:
python pyrate/utils/gdaldem.py tests/sydney_test/dem/sydney_trimmed.tif pyrate/utils/sydney_dem_color.txt out_tif.tif

Auto color:
python pyrate/utils/gdaldem.py tests/sydney_test/dem/sydney_trimmed.tif auto out_tif.tif

"""
import subprocess
import sys
import os
import tempfile
import numpy as np
from pyrate.shared import Ifg, DEM


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
        print '\nauto generating color file'
        color_file = gen_color_file(input_file)
    with open(color_file, 'r') as f:
        print '\ncolor file contents'
        print '='*50
        for l in f.readlines():
            print l
    print '='*50
    main(input_file, color_file, output_file)
