__author__ = 'sudipta'

import unittest
import subprocess
import tempfile
import os
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
        cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None', '-q', '-te'] \
              + extents_str

        for s in sydney_test_ifgs:
            temp_tif = tempfile.mktemp(suffix='.tif')
            t_cmd = cmd + [s.data_path,  temp_tif]
            subprocess.check_call(t_cmd)
            clipped_ref = gdal.Open(temp_tif).ReadAsArray()
            rast = gdal.Open(s.data_path)
            clipped = gdalwarp.crop(rast, extents)[0]
            np.testing.assert_array_almost_equal(clipped_ref, clipped)
            rast = None  # manual close
            os.remove(temp_tif)


class TestResample(unittest.TestCase):

    def test_sydney_data_resampling(self):
        sydney_test_ifgs = common.sydney_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        extents_str = [str(e) for e in extents]
        resolutions = [0.001666666, .001, 0.002, 0.0025, .01]

        for res in resolutions:
            res = [res, -res]
            cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None', '-q', '-te']\
                  + extents_str + ['-tr']
            new_res_str = [str(r) for r in res]
            cmd += new_res_str
            for s in sydney_test_ifgs:
                temp_tif = tempfile.mktemp(suffix='.tif')
                t_cmd = cmd + [s.data_path,  temp_tif]
                subprocess.check_call(t_cmd)
                resampled_ds = gdal.Open(temp_tif)
                resampled_ref = resampled_ds.ReadAsArray()

                rast = gdal.Open(s.data_path)
                reampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                    prefix='resampled_')
                resampled = gdalwarp.resample(s.data_path, extents, res,
                                                    reampled_temp_tif)
                np.testing.assert_array_almost_equal(resampled_ref, resampled)
                rast = None  # manual close
                os.remove(temp_tif)
                os.remove(reampled_temp_tif)


if __name__ == '__main__':
    unittest.main()
