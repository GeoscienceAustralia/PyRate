'''
Tests for the pyrate configuration file parser.

Created on 17/09/2012

.. codeauthor:: Ben Davies
'''

import unittest
from os.path import join
import os
import tempfile
from pyrate import config
from common import SYD_TEST_DIR, SYD_TEST_TIF, SYD_TEST_GAMMA


class ConfigTest(unittest.TestCase):
    def test_read_param_file(self):
        conf_path = join(SYD_TEST_DIR, 'pyrate.conf')
        params = config.get_config_params(conf_path)
        for k in params.keys():
            assert k and len(k) > 1
            assert params[k] != ''
            assert not k.endswith(":")  # are the colons removed?

    def test_read_param_file_no_refpixel(self):
        # ensure the parser can handle empty fields
        conf_path = join(SYD_TEST_DIR, 'pyrate2.conf')
        params = config.get_config_params(conf_path)

        assert params[config.REFX] == -1
        assert params[config.REFY] == -1

    def test_parse_namelist(self):
        nl = join(SYD_TEST_TIF, 'ifms_17')
        result = config.parse_namelist(nl)
        assert len(result) == 17
        files = ["geo_060619-061002.tif", "geo_060828-061211.tif",
                    "geo_061002-070430.tif", "geo_070115-070917.tif",
                    "geo_070219-070604.tif"]

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
        self.assertDictEqual(params, new_params)


if __name__ == "__main__":
    unittest.main()
