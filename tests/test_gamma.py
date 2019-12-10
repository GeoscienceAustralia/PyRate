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
This Python module contains tests for the gamma.py PyRate module.
"""
import glob
import os
import shutil
import tempfile
import unittest
from datetime import date, time
from os.path import join

import numpy as np
from numpy.testing import assert_array_almost_equal
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from osgeo import gdalconst
from osgeo import gdal_array

import core.ifgconstants as ifc
from core import shared, config as cf, gamma
from core.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR)
import prepifg, conv2tif
from core.shared import write_fullres_geotiff, GeotiffException
import common
from common import GAMMA_TEST_DIR, SML_TEST_GAMMA
from common import TEST_CONF_GAMMA, TEMPDIR
from common import small_data_setup
gdal.UseExceptions()

LIGHTSPEED = 3e8  # approx
PYRATEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class GammaCommandLineTests(unittest.TestCase):

    def setUp(self):
        self.base = join(PYRATEPATH,'tests', 'test_data', 'gamma')
        self.hdr = join(self.base, 'dem16x20raw.dem.par')
        temp_text = tempfile.mktemp()
        self.confFile = os.path.join(TEMPDIR,'{}/gamma_test.cfg'.format(temp_text))
        self.ifgListFile = os.path.join(TEMPDIR, '{}/gamma_ifg.list'.format(temp_text))
        self.base_dir = os.path.dirname(self.confFile)
        shared.mkdir_p(self.base_dir)

    def tearDown(self):
        try:
            os.remove(self.exp_path)
        except:
            pass
        shutil.rmtree(self.base_dir)

    def makeInputFiles(self, data):
        with open(self.confFile, 'w') as conf:
            conf.write('{}: {}\n'.format(DEM_HEADER_FILE, self.hdr))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '1'))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(SLC_DIR, ''))
        with open(self.ifgListFile, 'w') as ifgl:
            ifgl.write(data)

class GammaToGeoTiffTests(unittest.TestCase):
    """Tests conversion of GAMMA rasters to custom PyRate GeoTIFF"""

    @classmethod
    def setUpClass(cls):
        # create common combined header obj so the headers are only read once
        # tricker: needs both ifg headers, and DEM one for the extents
        filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
        hdr_paths = [join(GAMMA_TEST_DIR, f) for f in filenames]
        hdrs = [gamma.parse_epoch_header(p) for p in hdr_paths]
        dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')

        cls.DEM_HDR = gamma.parse_dem_header(dem_hdr_path)
        cls.COMBINED = gamma.combine_headers(*hdrs, dem_hdr=cls.DEM_HDR)

    def tearDown(self):
        if os.path.exists(self.dest):
            os.remove(self.dest)

    def test_to_geotiff_dem(self):
        hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
        hdr = gamma.parse_dem_header(hdr_path)
        data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
        self.dest = os.path.join(TEMPDIR, "tmp_gamma_dem.tif")

        write_fullres_geotiff(hdr, data_path, self.dest, nodata=0)
        exp_path = join(GAMMA_TEST_DIR, 'dem16x20_subset_from_gamma.tif')
        exp_ds = gdal.Open(exp_path)
        ds = gdal.Open(self.dest)

        # compare data and geographic headers
        # HACK: round  expected to nearest integer
        assert_array_almost_equal(np.rint(exp_ds.ReadAsArray()),
                                  ds.ReadAsArray())
        self.compare_rasters(exp_ds, ds)
        md = ds.GetMetadata()
        self.assertTrue(md['AREA_OR_POINT'] == 'Area')

    def test_to_geotiff_ifg(self):
        self.dest = os.path.join(TEMPDIR, 'tmp_gamma_ifg.tif')
        data_path = join(GAMMA_TEST_DIR,
                         '16x20_20090713-20090817_VV_4rlks_utm.unw')
        write_fullres_geotiff(self.COMBINED, data_path, self.dest, nodata=0)

        ds = gdal.Open(self.dest)
        exp_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.tif')
        exp_ds = gdal.Open(exp_path)

        # compare data and geographic headers
        assert_array_almost_equal(exp_ds.ReadAsArray(), ds.ReadAsArray())
        self.compare_rasters(ds, exp_ds)

        md = ds.GetMetadata()
        self.assertEqual(len(md), 11) # 11 metadata items
        self.assertTrue(md[ifc.MASTER_DATE] == str(date(2009, 7, 13)))
        self.assertTrue(md[ifc.SLAVE_DATE] == str(date(2009, 8, 17)))
        self.assertTrue(md[ifc.PYRATE_TIME_SPAN] == str(35 / ifc.DAYS_PER_YEAR))
        # TODO test failing
        # self.assertTrue(md[ifc.MASTER_TIME] == str(12))
        # TODO test failing
        # self.assertTrue(md[ifc.SLAVE_TIME] == str(time(12)))

        wavelen = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        self.assertAlmostEqual(wavelen, 0.05627457792190739)

    def test_to_geotiff_wrong_input_data(self):
        # use TIF, not UNW for data
        self.dest = os.path.join(TEMPDIR, 'tmp_gamma_ifg.tif')
        data_path = join(GAMMA_TEST_DIR,
                         '16x20_20090713-20090817_VV_4rlks_utm.tif')
        self.assertRaises(GeotiffException, write_fullres_geotiff,
                            self.COMBINED, data_path, self.dest, nodata=0)

    def test_mismatching_cell_resolution(self):
        hdrs = self.DEM_HDR.copy()
        hdrs[ifc.PYRATE_X_STEP] = 0.1  # fake a mismatch
        data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.unw')
        self.dest = os.path.join(TEMPDIR, 'fake')

        self.assertRaises(GeotiffException, write_fullres_geotiff, hdrs,
                            data_path, self.dest, 0)

    def compare_rasters(self, ds, exp_ds):
        band = ds.GetRasterBand(1)
        exp_band = exp_ds.GetRasterBand(1)

        nodata = band.GetNoDataValue()
        self.assertFalse(nodata is None)
        self.assertEqual(exp_band.GetNoDataValue(), nodata)

        pj = ds.GetProjection()
        self.assertTrue('WGS 84' in pj)
        self.assertEqual(exp_ds.GetProjection(), pj)
        for exp, act in zip(exp_ds.GetGeoTransform(), ds.GetGeoTransform()):
            self.assertAlmostEqual(exp, act, places=4)

    def test_bad_projection(self):
        hdr = self.DEM_HDR.copy()
        hdr[ifc.PYRATE_DATUM] = 'nonexistent projection'
        data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
        self.dest = os.path.join(TEMPDIR, 'tmp_gamma_dem2.tif')
        self.assertRaises(GeotiffException, write_fullres_geotiff, hdr,
                            data_path, self.dest, nodata=0)


class GammaHeaderParsingTests(unittest.TestCase):
    'Tests conversion of GAMMA headers to Py dicts'

    def test_parse_gamma_epoch_header(self):
        # minimal required headers are:
        # date:      2009  7 13
        # radar_frequency:        5.3310040e+09   Hz
        path = join(GAMMA_TEST_DIR, 'r20090713_VV.slc.par')
        hdrs = gamma.parse_epoch_header(path)

        exp_date = date(2009, 7, 13)
        self.assertEqual(hdrs[ifc.MASTER_DATE], exp_date)

        exp_wavelen = LIGHTSPEED / 5.3310040e+09
        self.assertEqual(hdrs[ifc.PYRATE_WAVELENGTH_METRES], exp_wavelen)

        incidence_angle = 22.9671
        self.assertEqual(hdrs[ifc.PYRATE_INCIDENCE_DEGREES], incidence_angle)

    def test_parse_gamma_dem_header(self):
        path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
        hdrs = gamma.parse_dem_header(path)

        self.assertEqual(hdrs[ifc.PYRATE_NCOLS], 16)
        self.assertEqual(hdrs[ifc.PYRATE_NROWS], 20)
        self.assertEqual(hdrs[ifc.PYRATE_LAT], -33.3831945)
        self.assertEqual(hdrs[ifc.PYRATE_LONG], 150.3870833)
        self.assertEqual(hdrs[ifc.PYRATE_X_STEP], 6.9444445e-05)
        self.assertEqual(hdrs[ifc.PYRATE_Y_STEP], -6.9444445e-05)


# Test data for the epoch header combination
H0 = {ifc.MASTER_DATE : date(2009, 7, 13),
      ifc.MASTER_TIME : time(12),
      ifc.PYRATE_WAVELENGTH_METRES: 1.8,
      ifc.PYRATE_INCIDENCE_DEGREES: 35.565,
      }

H1 = {ifc.MASTER_DATE : date(2009, 8, 17),
      ifc.MASTER_TIME : time(12, 10, 10),
      ifc.PYRATE_WAVELENGTH_METRES: 1.8,
      ifc.PYRATE_INCIDENCE_DEGREES: 35.56,
      }

H1_ERR1 = {ifc.MASTER_DATE : date(2009, 8, 17),
          ifc.MASTER_TIME : time(12),
          ifc.PYRATE_WAVELENGTH_METRES: 2.4,
          ifc.PYRATE_INCIDENCE_DEGREES: 35.56,
          }

H1_ERR2 = {ifc.MASTER_DATE : date(2009, 8, 17),
          ifc.MASTER_TIME : time(12),
          ifc.PYRATE_WAVELENGTH_METRES: 1.8,
          ifc.PYRATE_INCIDENCE_DEGREES: 35.76,
          }


class HeaderCombinationTests(unittest.TestCase):
    'Tests GAMMA epoch and DEM headers can be combined into a single Py dict'

    def setUp(self):
        self.err = gamma.GammaException
        dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
        self.dh = gamma.parse_dem_header(dem_hdr_path)

    def test_combine_headers(self):
        filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
        paths = [join(GAMMA_TEST_DIR, p) for p in filenames]
        hdr0, hdr1 = [gamma.parse_epoch_header(p) for p in paths]

        chdr = gamma.combine_headers(hdr0, hdr1, self.dh)

        exp_timespan = (18 + 17) / ifc.DAYS_PER_YEAR
        self.assertEqual(chdr[ifc.PYRATE_TIME_SPAN], exp_timespan)

        exp_date = date(2009, 7, 13)
        self.assertEqual(chdr[ifc.MASTER_DATE], exp_date)
        exp_date2 = date(2009, 8, 17)
        self.assertEqual(chdr[ifc.SLAVE_DATE], exp_date2)

        exp_wavelen = LIGHTSPEED / 5.3310040e+09
        self.assertEqual(chdr[ifc.PYRATE_WAVELENGTH_METRES], exp_wavelen)

    def test_fail_non_dict_header(self):
        self.assertRaises(self.err, gamma.combine_headers, H0, '', self.dh)
        self.assertRaises(self.err, gamma.combine_headers, '', H0, self.dh)
        self.assertRaises(self.err, gamma.combine_headers, H0, H1, None)
        self.assertRaises(self.err, gamma.combine_headers, H0, H1, '')

    def test_fail_mismatching_wavelength(self):
        self.assertRaises(self.err, gamma.combine_headers, H0, H1_ERR1, self.dh)
        
    def test_fail_mismatching_incidence(self):
        self.assertRaises(self.err, gamma.combine_headers, H0, H1_ERR2, self.dh)

    def test_fail_same_date(self):
        self.assertRaises(self.err, gamma.combine_headers, H0, H0, self.dh)

    def test_fail_bad_date_order(self):
        self.assertRaises(self.err, gamma.combine_headers, H1, H0, self.dh)


class TestGammaParallelVsSerial(unittest.TestCase):
    """
    Test Gamma prepifg produces correct results when run in serial and 
    parallel and that metadata is correctly set by both methods. These tests
    exclude the comparison of DEM files.
    """
    @classmethod
    def setUpClass(cls):
        # read in the params
        _, _, params = cf.get_ifg_paths(TEST_CONF_GAMMA)
        glob_prefix = "*utm_unw_1rlks_1cr.tif"

        # SERIAL
        cls.serial_dir = tempfile.mkdtemp()
        params[cf.OUT_DIR] = cls.serial_dir
        params[cf.PARALLEL] = False
        shared.mkdir_p(cls.serial_dir)

        gtif_paths = conv2tif.main(params)
        prepifg.main(params)

        serial_df = glob.glob(os.path.join(cls.serial_dir, glob_prefix))
        cls.serial_ifgs = small_data_setup(datafiles=serial_df)

        # Clean up serial converted tifs so we can test parallel conversion
        common.remove_tifs(SML_TEST_GAMMA)

        # PARALLEL 
        cls.parallel_dir = tempfile.mkdtemp()
        params[cf.OUT_DIR] = cls.parallel_dir
        params[cf.PARALLEL] = True
        shared.mkdir_p(cls.parallel_dir)

        gtif_paths = conv2tif.main(params)
        prepifg.main(params)
        
        para_df = glob.glob(os.path.join(cls.parallel_dir, glob_prefix))
        cls.para_ifgs = small_data_setup(datafiles=para_df)

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.parallel_dir)
        except PermissionError:
            print("File opened by another process.")

        try:
            shutil.rmtree(cls.serial_dir)
        except PermissionError:
            print("File opened by another process.")

        try:
            common.remove_tifs(SML_TEST_GAMMA)
        except PermissionError:
            print("File opened by another process.")

    def test_equality(self):
        for s, p in zip(self.serial_ifgs, self.para_ifgs):
            np.testing.assert_array_almost_equal(s.phase_data, p.phase_data)

    def test_meta_data_exist(self):
        for s, p in zip(self.serial_ifgs, self.para_ifgs):
            # all metadata equal
            self.assertDictEqual(s.meta_data, p.meta_data)
            # test that DATA_TYPE exists in metadata
            self.assertIn(ifc.DATA_TYPE, s.meta_data.keys())
            # test that DATA_TYPE is MULTILOOKED
            self.assertEqual(s.meta_data[ifc.DATA_TYPE], ifc.MULTILOOKED)


if __name__ == "__main__":
    unittest.main()
