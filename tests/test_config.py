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
'''
This Python module contains tests for the config.py PyRate module.
'''
import pytest
import os
import shutil
import tempfile
import unittest
from os.path import join

from .common import SYD_TEST_DIR, SYD_TEST_TIF, SYD_TEST_GAMMA
from pyrate import config
from pyrate.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR,
    LUIGI,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    NO_DATA_AVERAGING_THRESHOLD,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
    APS_METHOD,
    APS_CORRECTION)
from pyrate.tasks.utils import DUMMY_SECTION_NAME
from pyrate.compat import PyAPS_INSTALLED
from tests import common


class ConfigTest(unittest.TestCase):

    @staticmethod
    def test_read_param_file():
        conf_path = join(SYD_TEST_DIR, 'pyrate.conf')
        params = config.get_config_params(conf_path)
        for k in params.keys():
            assert k and len(k) > 1
            assert params[k] != ''
            assert not k.endswith(":")  # are the colons removed?

    @staticmethod
    def test_read_param_file_no_refpixel():
        # ensure the parser can handle empty fields
        conf_path = join(SYD_TEST_DIR, 'pyrate2.conf')
        params = config.get_config_params(conf_path)

        assert params[config.REFX] == -1
        assert params[config.REFY] == -1

    @staticmethod
    def test_parse_namelist():
        nl = join(SYD_TEST_TIF, 'ifms_17')
        result = list(config.parse_namelist(nl))
        assert len(result) == 17
        files = ["geo_060619-061002_unw.tif", "geo_060828-061211_unw.tif",
                    "geo_061002-070430_unw.tif", "geo_070115-070917_unw.tif",
                    "geo_070219-070604_unw.tif"]

        for path in files:
            assert path in result


class ConfigWriteTest(unittest.TestCase):

    def test_write_config_file(self):
        conf_path = join(SYD_TEST_GAMMA, 'pyrate_gamma.conf')
        params = config.get_config_params(conf_path)
        temp_config = tempfile.mktemp(suffix='.conf')
        config.write_config_file(params, temp_config)
        self.assertTrue(os.path.exists(temp_config))
        os.remove(temp_config)

    def test_new_config_file_and_original_match(self):
        conf_path = join(SYD_TEST_GAMMA, 'pyrate_gamma.conf')
        params = config.get_config_params(conf_path)
        temp_config = tempfile.mktemp(suffix='.conf')
        config.write_config_file(params, temp_config)
        new_params = config.get_config_params(temp_config)
        self.maxDiff = None
        self.assertDictEqual(params, new_params)
        os.remove(temp_config)


class ConfigAPSParametersTest(unittest.TestCase):
    def setUp(self):
        self.conf_path = common.SYDNEY_TEST_CONF
        self.params = config.get_config_params(self.conf_path)

    def test_incidence_and_elevation_keys_exist(self):
        self.assertIn(config.APS_INCIDENCE_MAP, self.params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, self.params.keys())

    def test_elevation_ext_should_not_exist(self):
        self.assertIn(config.APS_ELEVATION_EXT, self.params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, self.params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, self.params.keys())
        self.assertIsNone(self.params[config.APS_ELEVATION_MAP])

    def test_impedance_ext_should_exist(self):
        self.assertIn(config.APS_INCIDENCE_EXT, self.params.keys())

    def test_elevation_ext_keys_exist(self):
        self.assertIn(config.APS_INCIDENCE_EXT, self.params.keys())
        self.assertIn(config.APS_ELEVATION_EXT, self.params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, self.params.keys())

    def test_elevation_and_incidence_both_cant_have_values(self):
        self.assertIsNotNone(self.params[config.APS_INCIDENCE_MAP])
        self.assertIsNotNone(self.params[config.APS_INCIDENCE_EXT])
        self.assertIsNone(self.params[config.APS_ELEVATION_MAP])


class TestOneIncidenceOrElevationMap(unittest.TestCase):

    def setUp(self):
        self.base_dir = tempfile.mkdtemp()
        self.conf_file = tempfile.mktemp(suffix='.conf', dir=self.base_dir)
        self.ifgListFile = os.path.join(common.SYD_TEST_GAMMA, 'ifms_17')

    def tearDown(self):
        shutil.rmtree(self.base_dir)

    def make_input_files(self, inc='', ele=''):
        with open(self.conf_file, 'w') as conf:
            conf.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '1'))
            conf.write('{}: {}\n'.format(LUIGI, '0'))
            conf.write('{}: {}\n'.format(
                DEM_HEADER_FILE, os.path.join(
                    common.SYD_TEST_GAMMA, '20060619_utm_dem.par')))
            conf.write('{}: {}\n'.format(IFG_LKSX, '1'))
            conf.write('{}: {}\n'.format(IFG_LKSY, '1'))
            conf.write('{}: {}\n'.format(IFG_CROP_OPT, '1'))
            conf.write('{}: {}\n'.format(NO_DATA_AVERAGING_THRESHOLD, '0.5'))
            conf.write('{}: {}\n'.format(SLC_DIR, ''))
            conf.write('{}: {}\n'.format(DEM_FILE, common.SYD_TEST_DEM_GAMMA))
            conf.write('{}: {}\n'.format(APS_INCIDENCE_MAP, inc))
            conf.write('{}: {}\n'.format(APS_ELEVATION_MAP, ele))
            conf.write('{}: {}\n'.format(APS_CORRECTION, '1'))
            conf.write('{}: {}\n'.format(APS_METHOD, '2'))

    def test_inc_vs_ele_maps_inc_provided(self):
        self.make_input_files(inc=common.SYD_TEST_INCIDENCE)
        assert os.path.exists(self.conf_file)
        params = config.get_config_params(self.conf_file)
        # incidence variables
        self.assertIn(config.APS_INCIDENCE_MAP, params.keys())
        self.assertIn(config.APS_INCIDENCE_EXT, params.keys())
        self.assertIsNotNone(params[config.APS_INCIDENCE_MAP])
        self.assertIsNotNone(params[config.APS_INCIDENCE_EXT])

        # elevation variables
        self.assertIn(config.APS_ELEVATION_MAP, params.keys())
        self.assertIsNone(params[config.APS_ELEVATION_MAP])
        self.assertIn(config.APS_ELEVATION_EXT, params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, params.keys())

    def test_inc_vs_ele_maps_ele_provided(self):
        self.make_input_files(ele=common.SYD_TEST_ELEVATION)
        assert os.path.exists(self.conf_file)
        params = config.get_config_params(self.conf_file)
        # incidence variables
        self.assertIn(config.APS_INCIDENCE_MAP, params.keys())
        self.assertIn(config.APS_INCIDENCE_EXT, params.keys())
        self.assertIsNone(params[config.APS_INCIDENCE_MAP])
        self.assertIsNone(params[config.APS_INCIDENCE_EXT])

        # elevation variables
        self.assertIn(config.APS_ELEVATION_MAP, params.keys())
        self.assertIsNotNone(params[config.APS_ELEVATION_MAP])
        self.assertIn(config.APS_ELEVATION_EXT, params.keys())
        self.assertIn(config.APS_ELEVATION_MAP, params.keys())

    @pytest.mark.skipif(not PyAPS_INSTALLED,
                        reason='grib or pyaps not installed')
    def test_inc_vs_ele_maps_none_provided(self):
        self.make_input_files()
        assert os.path.exists(self.conf_file)
        self.assertRaises(config.ConfigException,
                          config.get_config_params, self.conf_file)


if __name__ == "__main__":
    unittest.main()
