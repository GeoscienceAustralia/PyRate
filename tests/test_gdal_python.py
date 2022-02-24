#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
This Python module contains tests for the gdal_python.py PyRate module.
"""
import os
import shutil
import subprocess
import tempfile
from copy import copy
from pathlib import Path

import numpy as np
from osgeo import gdal, gdalconst, osr

import pyrate.core
from pyrate import constants as c, conv2tif

from pyrate.configuration import MultiplePaths
from pyrate.core import gdal_python, ifgconstants as ifc
from pyrate.core.shared import Ifg
from tests import common


class TestResample(common.UnitTestAdaptation):

    def test_small_data_resampling(self):
        small_test_ifgs = common.small_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        extents_str = [str(e) for e in extents]
        resolutions = [0.001666666, .001, 0.002, 0.0025, .01]

        for res in resolutions:
            res = [res, -res]
            self.check_same_resampled_output(extents, extents_str, res,
                                             small_test_ifgs)

    def check_same_resampled_output(self, extents, extents_str, res,
                                    small_test_ifgs):
        cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None',
               '-q', '-r', 'near', '-te'] \
              + extents_str

        if res[0]:
            new_res_str = [str(r) for r in res]
            cmd += ['-tr'] + new_res_str
        for s in small_test_ifgs:
            temp_tif = tempfile.mktemp(suffix='.tif')
            t_cmd = cmd + [s.data_path, temp_tif]
            subprocess.check_call(t_cmd)
            resampled_ds = gdal.Open(temp_tif)
            resampled_ref = resampled_ds.ReadAsArray()

            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                 prefix='resampled_')
            resampled = gdal_python.resample_nearest_neighbour(s.data_path,
                                                               extents, res,
                                                               resampled_temp_tif)
            np.testing.assert_array_almost_equal(resampled_ref, resampled[0, :, :])
            try:
                os.remove(temp_tif)
            except PermissionError:
                print("File opened by another process.")

            try:
                os.remove(resampled_temp_tif)  # also proves file was written
            except PermissionError:
                print("File opened by another process.")

    def test_none_resolution_output(self):
        small_test_ifgs = common.small_data_setup()
        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        extents_str = [str(e) for e in extents]

        self.check_same_resampled_output(extents, extents_str, [None, None],
                                         small_test_ifgs)

    def test_output_file_written(self):
        small_test_ifgs = common.small_data_setup()
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        resolutions = [0.001666666, .001, 0.002, 0.0025, .01]
        for res in resolutions:
            for s in small_test_ifgs:
                resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                    prefix='resampled_')
                gdal_python.resample_nearest_neighbour(s.data_path, extents,
                                                       [res, -res],
                                                       resampled_temp_tif)
                self.assertTrue(os.path.exists(resampled_temp_tif))
                os.remove(resampled_temp_tif)

    def test_resampled_tif_has_metadata(self):
        small_test_ifgs = common.small_data_setup()

        # minX, minY, maxX, maxY = extents
        extents = [150.91, -34.229999976, 150.949166651, -34.17]
        for s in small_test_ifgs:

            resampled_temp_tif = tempfile.mktemp(suffix='.tif',
                                                prefix='resampled_')
            gdal_python.resample_nearest_neighbour(
                s.data_path, extents, [None, None], resampled_temp_tif)
            dst_ds = gdal.Open(resampled_temp_tif)
            md = dst_ds.GetMetadata()
            self.assertDictEqual(md, s.meta_data)
            try:
                os.remove(resampled_temp_tif)
            except PermissionError:
                print("File opened by another process.")


class TestBasicReampleTests(common.UnitTestAdaptation):

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


class TestMEMVsGTiff(common.UnitTestAdaptation):

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
            try:
                os.remove(temp_tif)
            except PermissionError:
                print("File opened by another process.")

    def test_mem(self):
        self.check('MEM')

    def test_gtiff(self):
        self.check('GTiff')


def test_coherence_files_not_converted():
    # define constants
    NO_DATA_VALUE = 0
    driver = gdal.GetDriverByName('GTiff')

    # create a sample gdal dataset

    # sample gdal dataset
    sample_gdal_filename = "dataset_01122000.tif"
    options = ['PROFILE=GeoTIFF']
    sample_gdal_dataset = driver.Create(sample_gdal_filename, 5, 5, 1, gdal.GDT_Float32, options=options)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    sample_gdal_dataset.SetProjection(wkt_projection)

    sample_gdal_band = sample_gdal_dataset.GetRasterBand(1)
    sample_gdal_band.SetNoDataValue(NO_DATA_VALUE)
    sample_gdal_band.WriteArray(np.arange(25).reshape(5, 5))
    sample_gdal_dataset.SetMetadataItem(ifc.FIRST_DATE, '2019-10-20')
    sample_gdal_dataset.SetMetadataItem(ifc.SECOND_DATE, '2019-11-01')
    sample_gdal_dataset.SetMetadataItem(ifc.PYRATE_WAVELENGTH_METRES, '10.05656')
    sample_gdal_dataset.FlushCache()
    sample_gdal_dataset = None
    ifg = Ifg(sample_gdal_filename)
    ifg.open()

    # create a coherence mask dataset
    tmpdir = tempfile.mkdtemp()
    out_dir = Path(tmpdir)  # we won't be creating any output coherence mask files as there are already GeoTIFFs
    params = common.min_params(out_dir)
    coherence_mask_filename = MultiplePaths(Path("mask_dataset_01122000-02122000.tif").as_posix(), params)
    coherence_mask_dataset = driver.Create(coherence_mask_filename.converted_path, 5, 5, 1, gdal.GDT_Float32)
    srs = osr.SpatialReference()
    wkt_projection = srs.ExportToWkt()
    coherence_mask_dataset.SetProjection(wkt_projection)
    coherence_mask_band = coherence_mask_dataset.GetRasterBand(1)
    coherence_mask_band.SetNoDataValue(NO_DATA_VALUE)
    arr = np.arange(0, 75, 3).reshape(5, 5) / 100.0
    arr[3, 4] = 0.25  # insert some random lower than threshold number
    arr[4, 2] = 0.20  # insert some random lower than threshold number

    coherence_mask_band.WriteArray(arr)
    # del the tmp handler datasets created
    del coherence_mask_dataset
    # create an artificial masked dataset
    expected_result_array = np.nan_to_num(
        np.array(
            [
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [np.nan, np.nan, np.nan, np.nan, np.nan],
                [10.0, 11.0, 12.0, 13.0, 14.0],
                [15.0, 16.0, 17.0, 18.0, np.nan],
                [20.0, 21.0, np.nan, 23.0, 24.0],
            ]
        )
    )

    # use the gdal_python.coherence_masking to find the actual mask dataset
    coherence_thresh = 0.3

    gdal_python.coherence_masking(ifg.dataset, coherence_mask_filename.converted_path, coherence_thresh)

    sample_gdal_array = np.nan_to_num(ifg.phase_data)

    # compare the artificial masked and actual masked datasets
    np.testing.assert_array_equal(sample_gdal_array, expected_result_array)

    # del the tmp datasets created
    os.remove(coherence_mask_filename.converted_path)

    ifg.close()
    os.remove(sample_gdal_filename)


def test_small_data_coherence(gamma_or_mexicoa_conf):
    work_dir = Path(tempfile.mkdtemp())
    params = common.manipulate_test_conf(conf_file=gamma_or_mexicoa_conf, work_dir=work_dir)

    params[c.COH_MASK] = 1

    ifg_multilist = copy(params[c.INTERFEROGRAM_FILES])
    conv2tif.main(params)

    for i in ifg_multilist:
        p = Path(i.converted_path)
        p.chmod(0o664)  # assign write permission as conv2tif output is readonly
        ifg = pyrate.core.shared.dem_or_ifg(data_path=p.as_posix())
        if not isinstance(ifg, Ifg):
            continue
        ifg.open()
        # now do coherence masking and compare
        ifg = pyrate.core.shared.dem_or_ifg(data_path=p.as_posix())
        ifg.open()
        converted_coh_file_path = pyrate.core.prepifg_helper.coherence_paths_for(p, params, tif=True)
        gdal_python.coherence_masking(ifg.dataset,
                                      coh_file_path=converted_coh_file_path,
                                      coh_thr=params[c.COH_THRESH]
                                      )
        nans = np.isnan(ifg.phase_data)
        coherence_path = pyrate.core.prepifg_helper.coherence_paths_for(p, params, tif=True)
        cifg = Ifg(coherence_path)
        cifg.open()
        cifg_below_thrhold = cifg.phase_data < params[c.COH_THRESH]
        np.testing.assert_array_equal(nans, cifg_below_thrhold)
    shutil.rmtree(work_dir)
