__author__ = 'sudipta'

import unittest
import subprocess
import tempfile
import os
import numpy as np
from numpy import where, nan, isclose
from osgeo import gdal, gdalconst
from pyrate import gdal_python as gdalwarp
from pyrate.tests import common
from pyrate import prepifg
from pyrate.shared import Ifg


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
            clipped = gdalwarp.crop(s.data_path, extents)[0]
            np.testing.assert_array_almost_equal(clipped_ref, clipped)
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
            self.check_same_resampled_output(extents, extents_str, res,
                                             sydney_test_ifgs)

    def check_same_resampled_output(self, extents, extents_str, res,
                                    sydney_test_ifgs):
        cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None',
               '-q', '-r', 'near', '-te'] \
              + extents_str

        if res[0]:
            new_res_str = [str(r) for r in res]
            cmd += ['-tr'] + new_res_str
        for s in sydney_test_ifgs:
            temp_tif = tempfile.mktemp(suffix='.tif')
            t_cmd = cmd + [s.data_path, temp_tif]
            subprocess.check_call(t_cmd)
            resampled_ds = gdal.Open(temp_tif)
            resampled_ref = resampled_ds.ReadAsArray()

            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                 prefix='resampled_')
            resampled = gdalwarp.resample_nearest_neighbour(s.data_path,
                                                            extents, res,
                                                            resampled_temp_tif)
            np.testing.assert_array_almost_equal(resampled_ref,
                                                 resampled[0, :, :])
            os.remove(temp_tif)
            os.remove(resampled_temp_tif)  # also proves file was written

    def test_none_resolution_output(self):
        sydney_test_ifgs = common.sydney_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        extents_str = [str(e) for e in extents]

        self.check_same_resampled_output(extents, extents_str, [None, None],
                                         sydney_test_ifgs)

    def test_output_file_written(self):
        sydney_test_ifgs = common.sydney_data_setup()
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        resolutions = [0.001666666, .001, 0.002, 0.0025, .01]
        for res in resolutions:
            for s in sydney_test_ifgs:
                resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                    prefix='resampled_')
                gdalwarp.resample_nearest_neighbour(s.data_path, extents,
                                                    [res, -res],
                                                    resampled_temp_tif)
                self.assertTrue(os.path.exists(resampled_temp_tif))
                os.remove(resampled_temp_tif)

    def test_sydney_data_crop_vs_resample(self):
        sydney_test_ifgs = common.sydney_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        for s in sydney_test_ifgs:
            clipped = gdalwarp.crop(s.data_path, extents)[0]
            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                prefix='resampled_')
            resampled = gdalwarp.resample_nearest_neighbour(
                s.data_path, extents, [None, None], resampled_temp_tif)
            self.assertTrue(os.path.exists(resampled_temp_tif))
            np.testing.assert_array_almost_equal(resampled[0, :, :], clipped)
            os.remove(resampled_temp_tif)


class BasicReampleTests(unittest.TestCase):

    def test_reproject_with_no_data(self):

        data = np.array([[2, 7],
                         [2, 7]])
        src_ds = gdal.GetDriverByName('MEM').Create('', 2, 2)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 1, 1)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_NearestNeighbour)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[7]])
        np.testing.assert_array_equal(got_data, expected_data)

    def test_reproject_with_no_data_2(self):

        data = np.array([[2, 7, 7, 7],
                         [2, 7, 7, 2]])
        height, width = data.shape
        src_ds = gdal.GetDriverByName('MEM').Create('', width, height)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 2, 1)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_NearestNeighbour)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[7, 3]])
        np.testing.assert_array_equal(got_data, expected_data)

    def test_reproject_with_no_data_3(self):

        data = np.array([[2, 7, 7, 7],
                         [2, 7, 7, 7],
                         [2, 7, 7, 7],
                         [2, 7, 7, 2],
                         [2, 7, 7, 2]])
        src_ds = gdal.GetDriverByName('MEM').Create('', 4, 5)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 2, 2)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_NearestNeighbour)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[7, 7],
                                  [7, 3]])
        np.testing.assert_array_equal(got_data, expected_data)

    def test_reproject_with_no_data_4(self):

        data = np.array([[2, 7, 7, 7, 2],
                         [2, 7, 7, 7, 2],
                         [2, 7, 7, 7, 2],
                         [2, 7, 7, 2, 2],
                         [2, 7, 7, 2, 2]])
        src_ds = gdal.GetDriverByName('MEM').Create('', 5, 5)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 2, 2)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_NearestNeighbour)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[7, 7],
                                  [7, 3]])
        np.testing.assert_array_equal(got_data, expected_data)


    def test_reproject_with_no_data_5(self):

        data = np.array([[2, 7, 7, 7, 2],
                         [2, 7, 7, 7, 2],
                         [2, 7, 7, 7, 2],
                         [2, 7, 7, 2, 2],
                         [2, 7, 7, 2, 2],
                         [2, 7, 7, 2, 2]])
        src_ds = gdal.GetDriverByName('MEM').Create('', 5, 6)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 2, 3)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_NearestNeighbour)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[7, 7],
                                  [7, 3],
                                  [7, 3]])
        np.testing.assert_array_equal(got_data, expected_data)


    def test_reproject_with_no_data_6(self):

        data = np.array([[4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 2, 2, 7.],
                         [4, 7, 7, 2, 2, 7.],
                         [4, 7, 7, 10, 2, 7.]], dtype=np.float32)
        src_ds = gdal.GetDriverByName('MEM').Create('', 6, 6, 1,
                                                    gdalconst.GDT_Float32)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 3, 3, 1,
                                                    gdalconst.GDT_Float32)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_Average)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[5.5, 7, 7],
                                  [5.5, 7, 7],
                                  [5.5, 8, 7]])
        np.testing.assert_array_equal(got_data, expected_data)

    def test_reproject_with_no_data_2bands(self):

        data = np.array([[[4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 7, 2, 7.],
                         [4, 7, 7, 2, 2, 7.],
                         [4, 7, 7, 2, 2, 7.],
                         [4, 7, 7, 10, 2, 7.]],
                        [[2, 0, 0, 0, 0, 0.],
                         [2, 0, 0, 0, 0, 2.],
                         [0, 1., 0, 0, 0, 1.],
                         [0, 0, 0, 0, 0, 2],
                         [0, 0, 0, 0, 0, 0.],
                         [0, 0, 0, 0, 0, 0.]]], dtype=np.float32)
        src_ds = gdal.GetDriverByName('MEM').Create('', 6, 6, 2,
                                                    gdalconst.GDT_Float32)

        src_ds.GetRasterBand(1).WriteArray(data[0, :, :])
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.GetRasterBand(2).WriteArray(data[1, :, :])
        # src_ds.GetRasterBand(1).SetNoDataValue()
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

        dst_ds = gdal.GetDriverByName('MEM').Create('', 3, 3, 2,
                                                    gdalconst.GDT_Float32)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_Average)
        got_data = dst_ds.GetRasterBand(1).ReadAsArray()
        expected_data = np.array([[5.5, 7, 7],
                                  [5.5, 7, 7],
                                  [5.5, 8, 7]])
        np.testing.assert_array_equal(got_data, expected_data)
        band2 = dst_ds.GetRasterBand(2).ReadAsArray()
        np.testing.assert_array_equal(band2, np.array([[1., 0., 0.5],
                                                       [0.25, 0., 0.75],
                                                       [0., 0., 0.]]))


if __name__ == '__main__':
    unittest.main()
