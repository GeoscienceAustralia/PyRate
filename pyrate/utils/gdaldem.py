"""
Based on: http://gis.stackexchange.com/questions/130199/changing-color-of-raster-images-based-on-their-data-values-gdal
Utility script to convert gray tif into color.
Usage: python gdaldem.py input_tif.tif color.txt output_color.tif

python pyrate/utils/gdaldem.py tests/sydney_test/dem/sydney_trimmed.tif pyrate/utils/sydney_dem_color.txt out_tif.tif

"""
import subprocess
import sys


def main(input_file, color_file, output_file):
    # TODO: implement auto color scaling option
    cmd = "gdaldem color-relief " + input_file \
          + ' ' + color_file + ' ' + output_file
    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    input_file = sys.argv[1]
    color_file = sys.argv[2]
    output_file = sys.argv[3]
    main(input_file, color_file, output_file)
