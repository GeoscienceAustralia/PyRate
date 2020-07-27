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
This Python module contains tests for the roipac.py PyRate module.
"""

import os
import shutil
import sys
import tempfile
from datetime import date
from os.path import exists, join

from numpy.testing import assert_array_almost_equal
from osgeo import gdal

import pyrate.core.ifgconstants as ifc
from pyrate.core import shared, roipac
from pyrate.core.config import (
    INPUT_IFG_PROJECTION,
    NO_DATA_VALUE,
    OBS_DIR,
    OUT_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    DEM_HEADER_FILE
)
# from pyrate.scripts.conv2tif import main as roipacMain
from pyrate.core.shared import GeotiffException
from pyrate.core.shared import write_fullres_geotiff
from tests.common import HEADERS_TEST_DIR, PREP_TEST_OBS, PREP_TEST_TIF
from tests.common import SML_TEST_DEM_DIR, SML_TEST_OBS, TEMPDIR, UnitTestAdaptation
from tests.common import SML_TEST_DEM_ROIPAC, SML_TEST_DEM_HDR

gdal.UseExceptions()


#  sanity checking
if not exists(HEADERS_TEST_DIR):
    sys.exit("ERROR: Missing the 'headers' data for unittests\n")

# constants
SHORT_HEADER_PATH = join(SML_TEST_OBS, 'geo_060619-061002.unw.rsc')
FULL_HEADER_PATH = join(HEADERS_TEST_DIR, "geo_060619-060828.unw.rsc")


class TestRoipacToGeoTiff(UnitTestAdaptation):
    """Tests conversion of GAMMA rasters to custom PyRate GeoTIFF"""

    @classmethod
    def setup_class(cls):
        hdr_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw.rsc')
        cls.HDRS = roipac.parse_header(hdr_path)

    @classmethod
    def teardown_class(cls):
        pass

    def test_to_geotiff_dem(self):
        hdr = roipac.parse_header(SML_TEST_DEM_HDR)        
        self.dest = os.path.join(TEMPDIR, "tmp_roipac_dem.tif")

        write_fullres_geotiff(hdr, SML_TEST_DEM_ROIPAC, self.dest, nodata=0)
        exp_path = join(SML_TEST_DEM_DIR, 'roipac_test_trimmed.tif')
        exp_ds = gdal.Open(exp_path)
        ds = gdal.Open(self.dest)

        # compare data and geographic headers
        assert_array_almost_equal(exp_ds.ReadAsArray(), ds.ReadAsArray())
        self.compare_rasters(ds, exp_ds)
        self.assertIsNotNone(ds.GetMetadata())

    def test_to_geotiff_ifg(self):
        # tricker: needs ifg header, and DEM one for extents
        hdrs = self.HDRS.copy()
        hdrs[ifc.PYRATE_DATUM] = 'WGS84'
        hdrs[ifc.DATA_TYPE] = ifc.ORIG

        self.dest = os.path.join('tmp_roipac_ifg.tif')
        data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
        write_fullres_geotiff(hdrs, data_path, self.dest, nodata=0)

        ds = gdal.Open(self.dest)
        band = ds.GetRasterBand(1)
        self.assertEqual(ds.RasterCount, 1)

        exp_path = join(PREP_TEST_TIF, 'geo_060619-061002_unw.tif')
        exp_ds = gdal.Open(exp_path)
        exp_band = exp_ds.GetRasterBand(1)

        # compare data and geographic headers
        assert_array_almost_equal(exp_band.ReadAsArray(), band.ReadAsArray())
        self.compare_rasters(ds, exp_ds)

        md = ds.GetMetadata()
        date1 = date(2006, 6, 19)
        date2 = date(2006, 10, 2)
        diff = (date2 - date1).days
        self.assertTrue(md[ifc.FIRST_DATE] == str(date1))
        self.assertTrue(md[ifc.SECOND_DATE] == str(date2))
        self.assertTrue(md[ifc.PYRATE_TIME_SPAN] == str(diff / ifc.DAYS_PER_YEAR))

        wavelen = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        self.assertAlmostEqual(wavelen, 0.0562356424)

    def test_to_geotiff_wrong_input_data(self):
        # ensure failure if TIF/other file used instead of binary UNW data
        self.dest = os.path.join(TEMPDIR, 'tmp_roipac_ifg.tif')
        data_path = join(PREP_TEST_TIF, 'geo_060619-061002_unw.tif')
        self.assertRaises(
            GeotiffException,
            write_fullres_geotiff,
            self.HDRS,
            data_path,
            self.dest,
            nodata=0)

    def test_bad_projection(self):
        hdrs = self.HDRS.copy()
        hdrs[ifc.PYRATE_DATUM] = 'bad datum string'
        hdrs[ifc.DATA_TYPE] = ifc.ORIG
        self.dest = os.path.join(TEMPDIR, 'tmp_roipac_ifg2.tif')
        data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
        self.assertRaises(GeotiffException, write_fullres_geotiff, hdrs, data_path, self.dest, 0)

    def test_mismatching_cell_resolution(self):
        hdrs = self.HDRS.copy()
        hdrs[ifc.PYRATE_X_STEP] = 0.1 # fake a mismatch
        hdrs[ifc.PYRATE_DATUM] = 'WGS84'
        data_path = join(PREP_TEST_OBS, 'geo_060619-061002.unw')
        self.dest = os.path.join(TEMPDIR, 'fake')

        self.assertRaises(GeotiffException, write_fullres_geotiff, hdrs, data_path, self.dest, 0)

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


class TestDateParsing(UnitTestAdaptation):

    def test_parse_short_date_pre2000(self):
        dstr = "980416"
        self.assertEqual(date(1998, 4, 16), roipac.parse_date(dstr))

    def test_parse_short_date_post2000(self):
        dstr = "081006"
        self.assertEqual(date(2008, 10, 6), roipac.parse_date(dstr))

    def test_parse_date_range(self):
        dstr = "980416-081006"
        exp = (date(1998, 4, 16), date(2008, 10, 6))
        self.assertEqual(exp, roipac.parse_date(dstr))


class TestHeaderParsingTests(UnitTestAdaptation):
    """Verifies ROIPAC headers are parsed correctly."""

    # short format header tests
    def test_parse_short_roipac_header(self):
        hdrs = roipac.parse_header(SHORT_HEADER_PATH)
        self.assertEqual(hdrs[ifc.PYRATE_NCOLS], 47)
        self.assertEqual(hdrs[ifc.PYRATE_NROWS], 72)
        self.assertAlmostEqual(hdrs[ifc.PYRATE_LONG], 150.910)
        self.assertEqual(hdrs[ifc.PYRATE_X_STEP], 0.000833333)
        self.assertEqual(hdrs[ifc.PYRATE_LAT], -34.170000000)
        self.assertEqual(hdrs[ifc.PYRATE_Y_STEP], -0.000833333)
        self.assertEqual(hdrs[ifc.PYRATE_WAVELENGTH_METRES], 0.0562356424)

    def test_parse_short_header_has_timespan(self):
        # Ensures TIME_SPAN_YEAR field is added during parsing
        hdrs = roipac.parse_header(SHORT_HEADER_PATH)
        self.assertIn(roipac.TIME_SPAN_YEAR, hdrs.keys())

        # check time span calc
        first = date(2006, 6, 19)
        second = date(2006, 10, 2)
        diff = (second - first).days / ifc.DAYS_PER_YEAR
        self.assertEqual(diff, hdrs[roipac.TIME_SPAN_YEAR])

    # long format header tests
    def test_parse_full_roipac_header(self):
        # Ensures "long style" original header can be parsed correctly
        hdrs = roipac.parse_header(FULL_HEADER_PATH)

        # check DATE/ DATE12 fields are parsed correctly
        date0 = date(2006, 6, 19) # from "DATE 060619" header
        date2 = date(2006, 8, 28) # from DATE12 060619-060828
        self.assertEqual(hdrs[ifc.FIRST_DATE], date0)
        self.assertEqual(hdrs[ifc.SECOND_DATE], date2)

    def test_read_full_roipac_header2(self):
        # Tests header from cropped original dataset is parsed correctly
        hdrs = roipac.parse_header(FULL_HEADER_PATH)
        self.assertTrue(len(hdrs) is not None)

    def test_xylast(self):
        # Test the X_LAST and Y_LAST header elements are calculated
        hdrs = roipac.parse_header(FULL_HEADER_PATH)
        self.assertAlmostEqual(hdrs[roipac.X_LAST], 151.8519444445)
        self.assertAlmostEqual(hdrs[roipac.Y_LAST], -34.625)
