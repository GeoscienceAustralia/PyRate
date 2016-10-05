__author__ = 'sudipta'

import unittest
import subprocess
import tempfile
import os
import shutil
import glob
import numpy as np
from numpy import where, nan, isclose
from osgeo import gdal, gdalconst, osr
from pyrate import gdal_python as gdalwarp
from pyrate.tests import common
from pyrate import prepifg
from pyrate.shared import Ifg, DEM


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

    def test_resampled_tif_has_metadata(self):
        sydney_test_ifgs = common.sydney_data_setup()

        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        for s in sydney_test_ifgs:

            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                prefix='resampled_')
            gdalwarp.resample_nearest_neighbour(
                s.data_path, extents, [None, None], resampled_temp_tif)
            dst_ds = gdal.Open(resampled_temp_tif)
            md = dst_ds.GetMetadata()
            self.assertDictEqual(md, s.meta_data)
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


    def test_reproject_average_resampling(self):

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

    def test_reproject_average_resampling_with_2bands(self):

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


class TestOldPrepifgVsGdalPython(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mktemp()
        shutil.copytree(common.SYD_TEST_TIF, self.test_dir)
        self.ifgs = common.sydney_data_setup(
            datafiles=glob.glob(os.path.join(self.test_dir, "*.tif")))
        self.ref_gtif = gdal.Open(self.ifgs[0].data_path, gdalconst.GA_ReadOnly)
        self.ref_proj = self.ref_gtif.GetProjection()
        self.ref_gt = self.ref_gtif.GetGeoTransform()
        self.data = self.ref_gtif.ReadAsArray()
        self.md = self.ref_gtif.GetMetadata()
        self.temp_tif = tempfile.mktemp(suffix='.tif')
        self.out_tif = tempfile.mktemp(suffix='.tif')

    def tearDown(self):
        if os.path.exists(self.temp_tif):
            os.remove(self.temp_tif)
        shutil.rmtree(self.test_dir)

    def test_gdal_python_vs_old_prepifg_prep(self):

        for i in range(10):
            data = np.array(np.random.randint(0, 3, size=(10, 10)),
                            dtype=np.float32)

            src_ds = gdal.GetDriverByName('GTiff').Create(self.temp_tif,
                                                          10, 10, 2,
                                                          gdalconst.GDT_Float32)

            src_ds.GetRasterBand(1).WriteArray(data)
            src_ds.GetRasterBand(1).SetNoDataValue(0)
            nan_matrix = where(data == 0, nan, data)
            src_ds.GetRasterBand(2).WriteArray(np.isnan(nan_matrix))
            src_ds.GetRasterBand(2).SetNoDataValue(-100)
            src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])

            dst_ds = gdal.GetDriverByName('MEM').Create('', 5, 5, 2,
                                                        gdalconst.GDT_Float32)
            dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

            for k, v in self.md.iteritems():
                src_ds.SetMetadataItem(k, v)
                dst_ds.SetMetadataItem(k, v)

            src_ds.FlushCache()
            gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_Average)
            nan_frac = dst_ds.GetRasterBand(2).ReadAsArray()
            avg = dst_ds.GetRasterBand(1).ReadAsArray()
            thresh = 0.5
            avg[nan_frac >= thresh] = np.nan

            ifg = Ifg(self.temp_tif)
            x_looks = y_looks = 2
            res = 2
            resolution = [res, -res]
            # minX, minY, maxX, maxY = extents
            # adfGeoTransform[0] /* top left x */
            # adfGeoTransform[1] /* w-e pixel resolution */
            # adfGeoTransform[2] /* 0 */
            # adfGeoTransform[3] /* top left y */
            # adfGeoTransform[4] /* 0 */
            # adfGeoTransform[5] /* n-s pixel resolution (negative value) */
            extents = [str(e) for e in [10, 0, 20, 10]]

            # only band 1 is resapled in warp_old
            data, self.old_prepifg_path = prepifg.warp_old(
                ifg, x_looks, y_looks, extents, resolution,
                thresh=thresh, crop_out=4, verbose=False)

            np.testing.assert_array_equal(data, avg)

    def test_gdal_python_vs_old_prepifg_prep2(self):

        for i in range(10):
            thresh = 0.5
            x_looks = y_looks = 2
            res = 2
            extents = [10, 0, 20, 10]
            extents_str = [str(e) for e in extents]
            data = np.array(np.random.randint(0, 3, size=(10, 10)),
                            dtype=np.float32)

            # create the  self.temp_tiff
            self.manipulation(data, self.temp_tif, self.md)

            # only band 1 is resapled in warp_old
            averaged_and_resampled, _ = gdalwarp.crop_resample_average(
                self.temp_tif, extents, [res, -res], self.out_tif, thresh)
            ifg = Ifg(self.temp_tif)
            # only band 1 is resampled in warp_old
            data, self.old_prepifg_path = prepifg.warp_old(
                ifg, x_looks, y_looks, extents_str, [res, -res],
                thresh=thresh, crop_out=4, verbose=False)

            np.testing.assert_array_almost_equal(
                data,
                averaged_and_resampled, decimal=4)

    @staticmethod
    def manipulation(data, tiff, md):
        src_ds = gdal.GetDriverByName('GTiff').Create(tiff,
                                                      10, 10, 2,
                                                      gdalconst.GDT_Float32)
        src_ds.GetRasterBand(1).WriteArray(data)
        src_ds.GetRasterBand(1).SetNoDataValue(0)
        nan_matrix = where(data == 0, nan, data)
        src_ds.GetRasterBand(2).WriteArray(np.isnan(nan_matrix))
        src_ds.GetRasterBand(2).SetNoDataValue(-100)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])
        dst_ds = gdal.GetDriverByName('MEM').Create('', 5, 5, 2,
                                                    gdalconst.GDT_Float32)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])
        for k, v in md.iteritems():
            src_ds.SetMetadataItem(k, v)
            dst_ds.SetMetadataItem(k, v)
        src_ds.FlushCache()
        return dst_ds, src_ds

    def test_gdal_python_vs_old_prepifg_no_match_pirate(self):

        for ifg in self.ifgs:
            extents = [150.91, -34.229999976, 150.949166651, -34.17]
            extents_str = [str(e) for e in extents]
            orig_res = 0.000833333
            thresh = 0.5

            for looks in range(10):
                x_looks = y_looks = looks
                res = orig_res*x_looks
                averaged_and_resapled, out_ds = gdalwarp.crop_resample_average(
                    ifg.data_path, extents, new_res=[res, -res],
                    output_file=self.temp_tif, thresh=thresh, match_pirate=False)

                # only band 1 is resampled in warp_old
                data, self.old_prepifg_path = prepifg.warp_old(
                    ifg, x_looks, y_looks, extents_str, [res, -res],
                    thresh=thresh, crop_out=4, verbose=False)
                yres, xres = data.shape

                # old_prepifg warp resample method loses one row at the bottom if
                # nrows % 2 == 1
                averaged_and_resapled = averaged_and_resapled[:yres, :xres]
                nrows, ncols = averaged_and_resapled.shape
                np.testing.assert_array_almost_equal(
                    data, averaged_and_resapled, decimal=4)

                # make sure they are the same after they are opened again
                # Last [yres:nrows, xres:ncols] won't match due to pirate
                # dropping last few rows/colums depending on resolution/looks
                data_from_file = gdal.Open(self.old_prepifg_path).ReadAsArray()
                new_from_file = out_ds.ReadAsArray()
                np.testing.assert_array_almost_equal(data_from_file[:yres, :xres],
                                                     new_from_file[:yres, :xres])
                out_ds = None  # manual close

    def test_gdal_python_vs_old_prepifg(self):

        for ifg in self.ifgs:
            extents = [150.91, -34.229999976, 150.949166651, -34.17]
            extents_str = [str(e) for e in extents]
            orig_res = 0.000833333
            thresh = 0.5
            for looks in range(10):
                x_looks = y_looks = looks
                res = orig_res*x_looks
                averaged_and_resampled = gdalwarp.crop_resample_average(
                    ifg.data_path, extents, new_res=[res, -res],
                    output_file=self.temp_tif, thresh=thresh)[0]

                # only band 1 is resampled in warp_old
                data, self.old_prepifg_path = prepifg.warp_old(
                    ifg, x_looks, y_looks, extents_str, [res, -res],
                    thresh=thresh, crop_out=4, verbose=False)
                yres, xres = data.shape

                # old_prepifg warp resample method loses one row at the bottom if
                # nrows % 2 == 1
                np.testing.assert_array_almost_equal(data, averaged_and_resampled[:yres, :xres],
                                                     decimal=4)
                nrows, ncols = averaged_and_resampled.shape
                # make sure they are the same after they are opened again
                data_from_file = gdal.Open(self.old_prepifg_path).ReadAsArray()
                self.assertTrue(os.path.exists(self.temp_tif))
                new_from_file = gdal.Open(self.temp_tif).ReadAsArray()
                np.testing.assert_array_almost_equal(data_from_file,
                                                     new_from_file)


    def test_no_out_file_when_driver_type_is_mem(self):
        for ifg in self.ifgs:
            extents = [150.91, -34.229999976, 150.949166651, -34.17]
            extents_str = [str(e) for e in extents]
            orig_res = 0.000833333
            thresh = 0.5
            x_looks = y_looks = 6
            res = orig_res*x_looks
            averaged_and_resampled, out_ds = gdalwarp.crop_resample_average(
                ifg.data_path, extents, new_res=[res, -res],
                output_file=self.temp_tif, thresh=thresh,
                out_driver_type='MEM')

            # only band 1 is resampled in warp_old
            data, self.old_prepifg_path = prepifg.warp_old(
                ifg, x_looks, y_looks, extents_str, [res, -res],
                thresh=thresh, crop_out=4, verbose=False)
            yres, xres = data.shape

            # old_prepifg warp resample method loses one row at the bottom if
            # nrows % 2 == 1
            averaged_and_resampled = averaged_and_resampled[:yres, :xres]
            np.testing.assert_array_almost_equal(data, averaged_and_resampled,
                                                 decimal=4)

            # make sure they are the same after they are opened again
            data_from_file = gdal.Open(self.old_prepifg_path).ReadAsArray()
            new_from_file = out_ds.ReadAsArray()
            np.testing.assert_array_almost_equal(data_from_file,
                                                 new_from_file)
            self.assertFalse(os.path.exists(self.temp_tif))
            out_ds = None  # manual close


class TestMEMVsGTiff(unittest.TestCase):

    @staticmethod
    def check(driver_type):

        temp_tif = tempfile.mktemp(suffix='.tif')

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
        src_ds = gdal.GetDriverByName(driver_type).Create(temp_tif, 6, 6, 2,
                                                    gdalconst.GDT_Float32)

        src_ds.GetRasterBand(1).WriteArray(data[0, :, :])
        src_ds.GetRasterBand(1).SetNoDataValue(2)
        src_ds.GetRasterBand(2).WriteArray(data[1, :, :])
        src_ds.GetRasterBand(2).SetNoDataValue(3)
        src_ds.SetGeoTransform([10, 1, 0, 10, 0, -1])
        src_ds.FlushCache()

        dst_ds = gdal.GetDriverByName('MEM').Create('', 3, 3, 2,
                                                    gdalconst.GDT_Float32)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.GetRasterBand(1).Fill(3)
        dst_ds.SetGeoTransform([10, 2, 0, 10, 0, -2])

        gdal.ReprojectImage(src_ds, dst_ds, '', '', gdal.GRA_Average)
        band1 = dst_ds.GetRasterBand(1).ReadAsArray()
        np.testing.assert_array_equal(band1, np.array([[5.5, 7, 7],
                                                       [5.5, 7, 7],
                                                       [5.5, 8, 7]]))
        band2 = dst_ds.GetRasterBand(2).ReadAsArray()
        np.testing.assert_array_equal(band2, np.array([[1., 0., 0.5],
                                                       [0.25, 0., 0.75],
                                                       [0., 0., 0.]]))
        if os.path.exists(temp_tif):
            os.remove(temp_tif)

    def test_mem(self):
        self.check('MEM')

    def test_gtiff(self):
        self.check('GTiff')


class TestGDalAverageResampleing(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()
