__author__ = 'sudipta'

import unittest
import subprocess
import tempfile
import numpy as np
from osgeo import gdal
from pyrate import gdal_python as gdalwarp
from pyrate.tests import common


class TestCrop(unittest.TestCase):
    def test_sydney_data_cropping(self):
        sydney_test_ifgs = common.sydney_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        extents_str = [str(e) for e in extents]
        cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None', '-q', '-te']

        for s in sydney_test_ifgs:
            temp_tif = tempfile.mktemp(suffix='.tif')
            t_cmd = cmd + extents_str + [s.data_path,  temp_tif]
            subprocess.check_call(t_cmd)
            clipped_ref = gdal.Open(temp_tif).ReadAsArray()
            rast = gdal.Open(s.data_path)
            clipped = gdalwarp.crop_raster(rast, extents)[0]
            np.testing.assert_array_almost_equal(clipped_ref, clipped)


if __name__ == '__main__':
    unittest.main()
