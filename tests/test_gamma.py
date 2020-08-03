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
This Python module contains tests for the gamma.py PyRate module.
"""
import os
import shutil
import tempfile
from datetime import date, time
from os.path import join
import pytest
from pathlib import Path
import numpy as np
from numpy.testing import assert_array_almost_equal
from osgeo import gdal

import pyrate.core.ifgconstants as ifc
from pyrate.core import shared, config as cf, gamma
from pyrate.core.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR)
from pyrate import prepifg, conv2tif
from pyrate.core.shared import write_fullres_geotiff, GeotiffException
from pyrate.constants import PYRATEPATH

from tests.common import manipulate_test_conf
from pyrate.configuration import Configuration
from tests.common import GAMMA_TEST_DIR
from tests.common import TEMPDIR
from tests.common import small_data_setup

gdal.UseExceptions()

LIGHTSPEED = 3e8  # approx


class TestGammaCommandLineTests:

    def setup_method(self):
        self.base = join(PYRATEPATH, 'tests', 'test_data', 'gamma')
        self.hdr = join(self.base, 'dem16x20raw.dem.par')
        temp_text = tempfile.mktemp()
        self.confFile = os.path.join(TEMPDIR,'{}/gamma_test.cfg'.format(temp_text))
        self.ifgListFile = os.path.join(TEMPDIR, '{}/gamma_ifg.list'.format(temp_text))
        self.base_dir = os.path.dirname(self.confFile)
        shared.mkdir_p(self.base_dir)

    def teardown_method(self):
        try:
            os.remove(self.exp_path)
        except:
            pass
        shutil.rmtree(self.base_dir, ignore_errors=True)

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


class TestGammaToGeoTiff:
    """Tests conversion of GAMMA rasters to custom PyRate GeoTIFF"""

    @classmethod
    def setup_method(cls):
        # create common combined header obj so the headers are only read once
        # tricker: needs both ifg headers, and DEM one for the extents
        filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
        hdr_paths = [join(GAMMA_TEST_DIR, f) for f in filenames]
        hdrs = [gamma.parse_epoch_header(p) for p in hdr_paths]
        dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')

        cls.DEM_HDR = gamma.parse_dem_header(dem_hdr_path)
        cls.COMBINED = gamma.combine_headers(*hdrs, dem_hdr=cls.DEM_HDR)

    def teardown_method(self):
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
        assert md['AREA_OR_POINT'] == 'Area'

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
        assert len(md) == 11 # 11 metadata items
        assert md[ifc.FIRST_DATE] == str(date(2009, 7, 13))
        assert md[ifc.SECOND_DATE] == str(date(2009, 8, 17))
        assert md[ifc.PYRATE_TIME_SPAN] == str(35 / ifc.DAYS_PER_YEAR)
        # TODO test failing
        # self.assertTrue(md[ifc.FIRST_DATE] == str(12))
        # TODO test failing
        # self.assertTrue(md[ifc.SECOND_DATE] == str(time(12)))

        wavelen = float(md[ifc.PYRATE_WAVELENGTH_METRES])
        assert wavelen == pytest.approx(0.05627457792190739)

    def test_to_geotiff_wrong_input_data(self):
        # use TIF, not UNW for data
        self.dest = os.path.join(TEMPDIR, 'tmp_gamma_ifg.tif')
        data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.tif')
        with pytest.raises(GeotiffException):
            write_fullres_geotiff(self.COMBINED, data_path, self.dest, nodata=0)

    def test_mismatching_cell_resolution(self):
        hdrs = self.DEM_HDR.copy()
        hdrs[ifc.PYRATE_X_STEP] = 0.1  # fake a mismatch
        data_path = join(GAMMA_TEST_DIR, '16x20_20090713-20090817_VV_4rlks_utm.unw')
        self.dest = os.path.join(TEMPDIR, 'fake')

        with pytest.raises(GeotiffException):
            write_fullres_geotiff(hdrs, data_path, self.dest, 0)

    def compare_rasters(self, ds, exp_ds):
        band = ds.GetRasterBand(1)
        exp_band = exp_ds.GetRasterBand(1)

        nodata = band.GetNoDataValue()
        assert ~(nodata is None)
        assert exp_band.GetNoDataValue() == nodata

        pj = ds.GetProjection()
        assert 'WGS 84' in pj
        assert exp_ds.GetProjection() == pj
        for exp, act in zip(exp_ds.GetGeoTransform(), ds.GetGeoTransform()):
            assert exp == pytest.approx(act, abs=0.0001)

    def test_bad_projection(self):
        hdr = self.DEM_HDR.copy()
        hdr[ifc.PYRATE_DATUM] = 'nonexistent projection'
        data_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem')
        self.dest = os.path.join(TEMPDIR, 'tmp_gamma_dem2.tif')
        with pytest.raises(GeotiffException):
            write_fullres_geotiff(hdr, data_path, self.dest, nodata=0)


class TestGammaHeaderParsingTests:
    'Tests conversion of GAMMA headers to Py dicts'

    def test_parse_gamma_epoch_header(self):
        # minimal required headers are:
        # date:      2009  7 13
        # radar_frequency:        5.3310040e+09   Hz
        path = join(GAMMA_TEST_DIR, 'r20090713_VV.slc.par')
        hdrs = gamma.parse_epoch_header(path)

        exp_date = date(2009, 7, 13)
        assert hdrs[ifc.FIRST_DATE] == exp_date

        exp_wavelen = LIGHTSPEED / 5.3310040e+09
        assert hdrs[ifc.PYRATE_WAVELENGTH_METRES] == exp_wavelen

        incidence_angle = 22.9671
        assert hdrs[ifc.PYRATE_INCIDENCE_DEGREES] == incidence_angle

    def test_parse_gamma_dem_header(self):
        path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
        hdrs = gamma.parse_dem_header(path)

        self.assert_equal(hdrs[ifc.PYRATE_NCOLS], 16)
        self.assert_equal(hdrs[ifc.PYRATE_NROWS], 20)
        self.assert_equal(hdrs[ifc.PYRATE_LAT], -33.3831945)
        self.assert_equal(hdrs[ifc.PYRATE_LONG], 150.3870833)
        self.assert_equal(hdrs[ifc.PYRATE_X_STEP], 6.9444445e-05)
        self.assert_equal(hdrs[ifc.PYRATE_Y_STEP], -6.9444445e-05)

    @staticmethod
    def assert_equal(arg1, arg2):
        assert arg1 == arg2


# Test data for the epoch header combination
H0 = {ifc.FIRST_DATE : date(2009, 7, 13),
      ifc.FIRST_DATE : time(12),
      ifc.PYRATE_WAVELENGTH_METRES: 1.8,
      ifc.PYRATE_INCIDENCE_DEGREES: 35.565,
      }

H1 = {ifc.FIRST_DATE : date(2009, 8, 17),
      ifc.FIRST_DATE : time(12, 10, 10),
      ifc.PYRATE_WAVELENGTH_METRES: 1.8,
      ifc.PYRATE_INCIDENCE_DEGREES: 35.56,
      }

H1_ERR1 = {ifc.FIRST_DATE : date(2009, 8, 17),
          ifc.FIRST_DATE : time(12),
          ifc.PYRATE_WAVELENGTH_METRES: 2.4,
          ifc.PYRATE_INCIDENCE_DEGREES: 35.56,
          }

H1_ERR2 = {ifc.FIRST_DATE : date(2009, 8, 17),
          ifc.FIRST_DATE : time(12),
          ifc.PYRATE_WAVELENGTH_METRES: 1.8,
          ifc.PYRATE_INCIDENCE_DEGREES: 35.76,
          }


class TestHeaderCombination:
    'Tests GAMMA epoch and DEM headers can be combined into a single Py dict'

    def setup_method(self):
        self.err = gamma.GammaException
        dem_hdr_path = join(GAMMA_TEST_DIR, 'dem16x20raw.dem.par')
        self.dh = gamma.parse_dem_header(dem_hdr_path)

    @staticmethod
    def assert_equal(arg1, arg2):
        assert arg1 == arg2

    def test_combine_headers(self):
        filenames = ['r20090713_VV.slc.par', 'r20090817_VV.slc.par']
        paths = [join(GAMMA_TEST_DIR, p) for p in filenames]
        hdr0, hdr1 = [gamma.parse_epoch_header(p) for p in paths]

        chdr = gamma.combine_headers(hdr0, hdr1, self.dh)

        exp_timespan = (18 + 17) / ifc.DAYS_PER_YEAR
        self.assert_equal(chdr[ifc.PYRATE_TIME_SPAN], exp_timespan)

        exp_date = date(2009, 7, 13)
        self.assert_equal(chdr[ifc.FIRST_DATE], exp_date)
        exp_date2 = date(2009, 8, 17)
        self.assert_equal(chdr[ifc.SECOND_DATE], exp_date2)

        exp_wavelen = LIGHTSPEED / 5.3310040e+09
        self.assert_equal(chdr[ifc.PYRATE_WAVELENGTH_METRES], exp_wavelen)

    def test_fail_non_dict_header(self):
        self.assertRaises(gamma.combine_headers, H0, '', self.dh)
        self.assertRaises(gamma.combine_headers, '', H0, self.dh)
        self.assertRaises(gamma.combine_headers, H0, H1, None)
        self.assertRaises(gamma.combine_headers, H0, H1, '')

    def test_fail_mismatching_wavelength(self):
        self.assertRaises(gamma.combine_headers, H0, H1_ERR1, self.dh)
        
    def test_fail_mismatching_incidence(self):
        self.assertRaises(gamma.combine_headers, H0, H1_ERR2, self.dh)

    def test_fail_same_date(self):
        self.assertRaises(gamma.combine_headers, H0, H0, self.dh)

    def test_fail_bad_date_order(self):
        with pytest.raises(self.err):
            gamma.combine_headers(H1, H0, self.dh)

    def assertRaises(self, func, * args):
        with pytest.raises(self.err):
            func(* args)


glob_prefix = "*utm_unw_ifg_1rlks_1cr.tif"


@pytest.fixture(scope='module')
def parallel_ifgs(gamma_conf):

    tdir = Path(tempfile.mkdtemp())

    params_p = manipulate_test_conf(gamma_conf, tdir)
    params_p[cf.PARALLEL] = 1

    output_conf_file = 'conf.conf'
    output_conf = tdir.joinpath(output_conf_file)
    cf.write_config_file(params=params_p, output_conf_file=output_conf)

    params_p = Configuration(output_conf).__dict__

    conv2tif.main(params_p)
    prepifg.main(params_p)

    parallel_df = list(Path(tdir).joinpath('out').glob(glob_prefix))

    p_ifgs = small_data_setup(datafiles=parallel_df)
    yield p_ifgs

    shutil.rmtree(params_p[cf.OBS_DIR], ignore_errors=True)


@pytest.fixture(scope='module')
def series_ifgs(gamma_conf):

    print('======================setup series==========================')
    tdir = Path(tempfile.mkdtemp())

    params_s = manipulate_test_conf(gamma_conf, tdir)
    params_s[cf.PARALLEL] = 0

    output_conf_file = 'conf.conf'
    output_conf = tdir.joinpath(output_conf_file)
    cf.write_config_file(params=params_s, output_conf_file=output_conf)

    params_s = Configuration(output_conf).__dict__

    gtif_paths = conv2tif.main(params_s)
    prepifg.main(params_s)

    parallel_df = list(Path(tdir).joinpath('out').glob(glob_prefix))

    s_ifgs = small_data_setup(datafiles=parallel_df)

    yield s_ifgs

    print('======================teardown series==========================')

    shutil.rmtree(params_s[cf.OBS_DIR], ignore_errors=True)


def test_equality(series_ifgs, parallel_ifgs):
    for s, p in zip(series_ifgs, parallel_ifgs):
        np.testing.assert_array_almost_equal(s.phase_data, p.phase_data)


def test_meta_data_exist(series_ifgs, parallel_ifgs):
    for i, (s, p) in enumerate(zip(series_ifgs, parallel_ifgs)):
        print(s, p)
        # all metadata equal
        assert s.meta_data == p.meta_data
        # test that DATA_TYPE exists in metadata
        assert ifc.DATA_TYPE in s.meta_data.keys()
        # test that DATA_TYPE is MULTILOOKED
        assert s.meta_data[ifc.DATA_TYPE] == ifc.MULTILOOKED

    assert i == 16
