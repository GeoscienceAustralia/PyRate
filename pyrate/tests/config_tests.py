'''
Tests for the pyrate configuration file parser.

Created on 17/09/2012
@author: Ben Davies, NCI
'''

import sys, unittest
from os.path import exists, join

from pyrate import config
from common import SYD_TEST_DIR, SYD_TEST_TIF


def test_read_param_file():
	conf_path = join(SYD_TEST_DIR, 'pyrate.conf')
	params = config.get_config_params(conf_path)
	for k in params.keys():
		assert k and len(k) > 1
		assert params[k] != ''
		assert not k.endswith(":") # are the colons removed?

def test_read_param_file_no_refpixel():
	# ensure the parser can handle empty fields
	conf_path = join(SYD_TEST_DIR, 'pyrate2.conf')
	params = config.get_config_params(conf_path)

	assert params[config.REFX] == -1
	assert params[config.REFY] == -1

def test_parse_namelist():
	nl = join(SYD_TEST_TIF, 'ifms_17')
	result = config.parse_namelist(nl)
	assert len(result) == 17
	files = ["geo_060619-061002.tif", "geo_060828-061211.tif",
				"geo_061002-070430.tif", "geo_070115-070917.tif",
				"geo_070219-070604.tif" ]

	for path in files:
		assert path in result
