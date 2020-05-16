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
#
# pylint: disable=trailing-whitespace, missing-docstring
'''
This Python module contains tests for the config.py PyRate module.
'''
import os
import shutil
import tempfile
import unittest
from os.path import join

import pytest

from tests.common import SML_TEST_CONF, SML_TEST_TIF
from tests.common import TEST_CONF_ROIPAC, TEST_CONF_GAMMA
from pyrate.core import config
from pyrate.core.config import (
    validate_parameters, validate_optional_parameters, validate_epochs,
    validate_obs_thresholds, validate_ifgs, validate_gamma_headers,
    validate_coherence_files, validate_tifs_exist, validate_minimum_epochs,
    validate_epoch_thresholds, validate_epoch_cutoff,
    validate_crop_parameters, validate_slpf_cutoff,
    validate_reference_pixel_search_windows,
    validate_multilook_parameters,
    validate_prepifg_tifs_exist,
    _get_temporal_info, _get_prepifg_info, write_config_file)

from pyrate.constants import SIXTEEN_DIGIT_EPOCH_PAIR, TWELVE_DIGIT_EPOCH_PAIR, \
    EIGHT_DIGIT_EPOCH

from pyrate.core.config import (
    _COHERENCE_VALIDATION,
    _ORBITAL_FIT_VALIDATION,
    _APSEST_VALIDATION,
    _TIME_SERIES_VALIDATION,
    _PARAM_VALIDATION,
    _GAMMA_VALIDATION,
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    SLC_DIR,
    HDR_FILE_LIST,
    COH_MASK,
    COH_THRESH,
    COH_FILE_DIR,
    COH_FILE_LIST,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    IFG_XFIRST, IFG_XLAST,
    IFG_YFIRST, IFG_YLAST,
    REFX, REFY,
    REFNX,
    REFNY,
    REF_CHIP_SIZE,
    REF_MIN_FRAC,
    ORBITAL_FIT,
    ORBITAL_FIT_METHOD,
    ORBITAL_FIT_DEGREE,
    ORBITAL_FIT_LOOKS_X,
    ORBITAL_FIT_LOOKS_Y,
    LR_NSIG,
    LR_MAXSIG,
    LR_PTHRESH,
    APSEST,
    TLPF_METHOD,
    TLPF_CUTOFF,
    TLPF_PTHR,
    SLPF_METHOD,
    SLPF_CUTOFF,
    SLPF_ORDER,
    SLPF_NANFILL,
    TIME_SERIES_CAL,
    TIME_SERIES_PTHRESH,
    TIME_SERIES_SM_FACTOR,
    TIME_SERIES_SM_ORDER,
    TIME_SERIES_METHOD,
    PARALLEL,
    PROCESSES,
    NAN_CONVERSION,
    NO_DATA_AVERAGING_THRESHOLD,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
    APS_METHOD,
    APS_CORRECTION,
    ConfigException)
from tests import common
from pyrate.configuration import Configuration
DUMMY_SECTION_NAME = 'pyrate'


class ValidateTestConfig(unittest.TestCase):

    def test_gamma_conf_passes(self):
        config.get_config_params(TEST_CONF_GAMMA)
        
    def test_roipac_conf_passes(self):
        config.get_config_params(TEST_CONF_ROIPAC)


class TestConfigValidation(unittest.TestCase):
    def setUp(self):
        """
        Get a copy of the GAMMA params and also use this to verify that 
        they are correct before we start testing.
        """
        self.params = Configuration(TEST_CONF_GAMMA).__dict__
        self.roipac_params = config.get_config_params(TEST_CONF_ROIPAC)
        self.dummy_dir = '/i/should/not/exist/'
        if os.path.exists(self.dummy_dir):
            raise IOError("'dummy_dir' needs to be non-existant for testing.")
    
    def test_validators(self):
        """
        Test validation functions for 'compulsory' parameters.
        """
        def validate(key, value):
            return _PARAM_VALIDATION[key][0](value)

        self.assertTrue(validate(IFG_FILE_LIST, self.params[IFG_FILE_LIST]))
        self.assertFalse(validate(IFG_FILE_LIST, None))
        self.assertFalse(validate(IFG_FILE_LIST, self.dummy_dir))

        self.assertTrue(validate(DEM_FILE, self.params[DEM_FILE]))
        self.assertFalse(validate(DEM_FILE, None))
        self.assertFalse(validate(DEM_FILE, self.dummy_dir))

        self.assertTrue(validate(DEM_HEADER_FILE, self.params[DEM_HEADER_FILE]))
        self.assertFalse(validate(DEM_HEADER_FILE, None))
        self.assertFalse(validate(DEM_HEADER_FILE, self.dummy_dir))

        self.assertTrue(validate(OUT_DIR, self.params[OUT_DIR]))
        self.assertFalse(validate(OUT_DIR, None))
        # OUT_DIR gets created at runtime
        self.assertTrue(validate(OUT_DIR, self.dummy_dir))

        self.assertTrue(validate(APS_INCIDENCE_MAP, self.params[APS_INCIDENCE_MAP]))
        self.assertFalse(validate(APS_INCIDENCE_MAP, self.dummy_dir))
        self.assertTrue(validate(APS_INCIDENCE_MAP, None))

        self.assertTrue(validate(APS_ELEVATION_MAP, self.params[APS_ELEVATION_MAP]))
        self.assertFalse(validate(APS_ELEVATION_MAP, self.dummy_dir))
        self.assertTrue(validate(APS_ELEVATION_MAP, None))

        self.assertTrue(validate(IFG_CROP_OPT, 1))
        self.assertTrue(validate(IFG_CROP_OPT, 2))
        self.assertTrue(validate(IFG_CROP_OPT, 3))
        self.assertTrue(validate(IFG_CROP_OPT, 4))
        self.assertFalse(validate(IFG_CROP_OPT, 0))
        self.assertFalse(validate(IFG_CROP_OPT, 5))

        self.assertTrue(validate(IFG_LKSX, self.params[IFG_LKSX]))
        self.assertFalse(validate(IFG_LKSX, 0))

        self.assertTrue(validate(IFG_LKSY, self.params[IFG_LKSY]))
        self.assertFalse(validate(IFG_LKSY, 0))

        # TODO: IFG_XFIRST, IFG_XLAST, IFG_YFIRST, IFG_YLAST

        self.assertTrue(validate(NO_DATA_VALUE, self.params[NO_DATA_VALUE]))

        self.assertTrue(validate(COH_MASK, 0))
        self.assertTrue(validate(COH_MASK, 1))
        self.assertFalse(validate(COH_MASK, -1))
        self.assertFalse(validate(COH_MASK, 2))

        self.assertTrue(validate(ORBITAL_FIT, 0))
        self.assertTrue(validate(ORBITAL_FIT, 1))
        self.assertFalse(validate(ORBITAL_FIT, -1))
        self.assertFalse(validate(ORBITAL_FIT, 2))

        self.assertTrue(validate(LR_NSIG, self.params[LR_NSIG]))
        self.assertFalse(validate(LR_NSIG, 0))
        self.assertFalse(validate(LR_NSIG, 11))

        self.assertTrue(validate(LR_PTHRESH, self.params[LR_PTHRESH]))
        self.assertFalse(validate(LR_PTHRESH, 0))

        self.assertTrue(validate(LR_MAXSIG, self.params[LR_MAXSIG]))
        self.assertFalse(validate(LR_MAXSIG, -1))
        self.assertFalse(validate(LR_MAXSIG, 1001))

        self.assertTrue(validate(APSEST, 0))
        self.assertTrue(validate(APSEST, 1))
        self.assertFalse(validate(APSEST, -1))
        self.assertFalse(validate(APSEST, 2))

        self.assertTrue(validate(TIME_SERIES_CAL, 0))
        self.assertTrue(validate(TIME_SERIES_CAL, 1))
        self.assertFalse(validate(TIME_SERIES_CAL, -1))
        self.assertFalse(validate(TIME_SERIES_CAL, 2))

        self.assertTrue(validate(PARALLEL, 0))
        self.assertTrue(validate(PARALLEL, 1))
        self.assertFalse(validate(PARALLEL, 2))
        self.assertFalse(validate(PARALLEL, -1))
        self.assertFalse(validate(PARALLEL, 3))

        self.assertTrue(validate(PROCESSES, 1))
        self.assertFalse(validate(PROCESSES, -1))
        self.assertFalse(validate(PROCESSES, 0))

        self.assertTrue(validate(PROCESSOR, 0))
        self.assertTrue(validate(PROCESSOR, 1))
        self.assertTrue(validate(PROCESSOR, 2))
        self.assertFalse(validate(PROCESSOR, -1))
        self.assertFalse(validate(PROCESSOR, 3))

        self.assertTrue(validate(NAN_CONVERSION, 0))
        self.assertTrue(validate(NAN_CONVERSION, 1))
        self.assertFalse(validate(NAN_CONVERSION, -1))
        self.assertFalse(validate(NAN_CONVERSION, 2))

        self.assertTrue(validate(NO_DATA_AVERAGING_THRESHOLD,
                                 self.params[NO_DATA_AVERAGING_THRESHOLD]))

    def test_gamma_validators(self):
        def validate(key, value):
            return _GAMMA_VALIDATION[key][0](value)

        self.assertFalse(validate(HDR_FILE_LIST, None))
        self.assertFalse(validate(HDR_FILE_LIST, self.dummy_dir))

    def test_coherence_validators(self):
        def validate(key, value):
            return _COHERENCE_VALIDATION[key][0](value)

        self.assertTrue(validate(COH_THRESH, 0.1))
        self.assertFalse(validate(COH_THRESH, -0.1))
        self.assertFalse(validate(COH_THRESH, 1.1))

        self.assertFalse(validate(COH_FILE_LIST, None))

    def test_orbital_validators(self):
        def validate(key, value):
            return _ORBITAL_FIT_VALIDATION[key][0](value)

        self.assertTrue(validate(ORBITAL_FIT_METHOD, 1))
        self.assertTrue(validate(ORBITAL_FIT_METHOD, 2))
        self.assertFalse(validate(ORBITAL_FIT_METHOD, 0))
        self.assertFalse(validate(ORBITAL_FIT_METHOD, 3))

        self.assertTrue(validate(ORBITAL_FIT_DEGREE, 1))
        self.assertTrue(validate(ORBITAL_FIT_DEGREE, 2))
        self.assertTrue(validate(ORBITAL_FIT_DEGREE, 3))
        self.assertFalse(validate(ORBITAL_FIT_DEGREE, 0))
        self.assertFalse(validate(ORBITAL_FIT_DEGREE, 4))

        self.assertFalse(validate(ORBITAL_FIT_LOOKS_X, 0))

        self.assertFalse(validate(ORBITAL_FIT_LOOKS_Y, 0))

    def test_apsest_validators(self):
        def validate(key, value):
            return _APSEST_VALIDATION[key][0](value)

        for i in range(1, 4):
            self.assertTrue(validate(TLPF_METHOD, i))
        self.assertFalse(validate(TLPF_METHOD, 0))
        self.assertFalse(validate(TLPF_METHOD, 4))

        self.assertFalse(validate(TLPF_CUTOFF, 0.0026))
        self.assertTrue(validate(TLPF_CUTOFF, 0.0028))

        self.assertFalse(validate(TLPF_PTHR, 0))
        self.assertTrue(validate(TLPF_PTHR, 1))

        self.assertTrue(validate(SLPF_METHOD, 1))
        self.assertTrue(validate(SLPF_METHOD, 2))

        self.assertTrue(validate(SLPF_CUTOFF, 0.001))
        self.assertFalse(validate(SLPF_CUTOFF, 0.0))

        for i in range(1, 4):
            self.assertTrue(validate(SLPF_ORDER, i))
        self.assertFalse(validate(SLPF_ORDER, 0))
        self.assertFalse(validate(SLPF_ORDER, 4))

        self.assertTrue(validate(SLPF_NANFILL, 0))
        self.assertTrue(validate(SLPF_NANFILL, 1))
        self.assertFalse(validate(SLPF_NANFILL, -1))
        self.assertFalse(validate(SLPF_NANFILL, 2))

    def test_time_series_validators(self):
        def validate(key, value):
            return _TIME_SERIES_VALIDATION[key][0](value)

        self.assertTrue(validate(TIME_SERIES_PTHRESH, 1))
        self.assertFalse(validate(TIME_SERIES_PTHRESH, 0))

        self.assertTrue(validate(TIME_SERIES_SM_FACTOR, -1.0))
        self.assertFalse(validate(TIME_SERIES_SM_FACTOR, 0.1))
        self.assertFalse(validate(TIME_SERIES_SM_FACTOR, -5.1))

        self.assertTrue(validate(TIME_SERIES_SM_ORDER, 1))
        self.assertTrue(validate(TIME_SERIES_SM_ORDER, 2))
        self.assertFalse(validate(TIME_SERIES_SM_ORDER, 0))
        self.assertFalse(validate(TIME_SERIES_SM_ORDER, 3))

        self.assertTrue(validate(TIME_SERIES_METHOD, 1))
        self.assertTrue(validate(TIME_SERIES_METHOD, 2))
        self.assertFalse(validate(TIME_SERIES_METHOD, 0))
        self.assertFalse(validate(TIME_SERIES_METHOD, 3))

    def test_gamma_validation_only_performed_when_on(self):
        self.params[PROCESSOR] = 1
        self.params[HDR_FILE_LIST] = None
        with pytest.raises(ConfigException):
            validate_optional_parameters(self.params)

        self.params[PROCESSOR] = 0
        validate_parameters(self.params)

    def test_coherence_validation_only_performed_when_on(self):
        self.params[COH_MASK] = 1
        self.params[COH_FILE_LIST] = None
        with pytest.raises(ConfigException):
            validate_optional_parameters(self.params)

        self.params[COH_MASK] = 0
        validate_optional_parameters(self.params)

    def test_orbital_validation_only_performed_when_on(self):
        self.params[ORBITAL_FIT] = 1
        self.params[ORBITAL_FIT_METHOD] = 0
        with pytest.raises(ConfigException):
            validate_optional_parameters(self.params)

        self.params[ORBITAL_FIT] = 0
        validate_optional_parameters(self.params)

    def test_apsest_validation_only_performed_when_on(self):
        self.params[APSEST] = 1
        self.params[TLPF_METHOD] = 0
        with pytest.raises(ConfigException):
            validate_optional_parameters(self.params)

        self.params[APSEST] = 0
        validate_optional_parameters(self.params)

    def test_time_series_validation_only_performed_when_on(self):
        self.params[TIME_SERIES_CAL] = 1
        self.params[TIME_SERIES_PTHRESH] = 0
        with pytest.raises(ConfigException):
            validate_optional_parameters(self.params)

        self.params[TIME_SERIES_CAL] = 0
        validate_optional_parameters(self.params)

    def test_epochs_in_gamma_obs(self):
        validate_epochs(self.params[IFG_FILE_LIST], SIXTEEN_DIGIT_EPOCH_PAIR)
        self.params[IFG_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_ifms_17')
        with pytest.raises(ConfigException):
            validate_epochs(self.params[IFG_FILE_LIST], SIXTEEN_DIGIT_EPOCH_PAIR)

    def test_epochs_in_roipac_obs(self):
        validate_epochs(self.roipac_params[IFG_FILE_LIST], TWELVE_DIGIT_EPOCH_PAIR)
        self.roipac_params[IFG_FILE_LIST] = os.path.join(common.SML_TEST_OBS, 'bad_epochs_ifms_17')
        with pytest.raises(ConfigException):
            validate_epochs(self.roipac_params[IFG_FILE_LIST], TWELVE_DIGIT_EPOCH_PAIR)

    def test_epochs_in_gamma_headers(self):
        validate_epochs(self.params[HDR_FILE_LIST], EIGHT_DIGIT_EPOCH)
        self.params[HDR_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_headers')
        with pytest.raises(ConfigException):
            validate_epochs(self.params[HDR_FILE_LIST], EIGHT_DIGIT_EPOCH)

    def test_epochs_in_coherence_files(self):
        self.params[COH_MASK] = 1
        self.params[COH_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'ifms_17')
        validate_epochs(self.params[COH_FILE_LIST], SIXTEEN_DIGIT_EPOCH_PAIR)
        self.params[COH_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_ifms_17')
        with pytest.raises(ConfigException):
            validate_epochs(self.params[COH_FILE_LIST], SIXTEEN_DIGIT_EPOCH_PAIR)

    def test_ifgs_exist(self):
        validate_ifgs(self.params[IFG_FILE_LIST])
        self.params[IFG_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_ifms_17')
        with pytest.raises(ConfigException):
            validate_ifgs(self.params[IFG_FILE_LIST])

    def test_gamma_headers_exist(self):
        validate_gamma_headers(self.params[IFG_FILE_LIST], self.params[HDR_FILE_LIST])
        self.params[HDR_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_headers')
        with pytest.raises(ConfigException):
            validate_gamma_headers(self.params[IFG_FILE_LIST], self.params[HDR_FILE_LIST])

    def test_coherence_files_exist(self):
        self.params[COH_MASK] = 1
        self.params[COH_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'ifms_17')
        validate_coherence_files(self.params[IFG_FILE_LIST], self.params)
        self.params[COH_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'bad_epochs_ifms_17')
        tmp_config = tempfile.mktemp('.conf')
        write_config_file(self.params, output_conf_file=tmp_config)
        with pytest.raises(ConfigException):
            Configuration(tmp_config).__dict__

    def test_obs_thresholds(self):
        self.params[TIME_SERIES_PTHRESH] = 16
        self.params[TLPF_PTHR] = 16
        validate_obs_thresholds(self.params[IFG_FILE_LIST], self.params)

        self.params[TIME_SERIES_PTHRESH] = 18
        self.params[TLPF_PTHR] = 18
        with pytest.raises(ConfigException):
            validate_obs_thresholds(self.params[IFG_FILE_LIST], self.params)


class TestConfigValidationWithFullResGeotiffs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pyrate import conv2tif
        cls.params = Configuration(TEST_CONF_GAMMA).__dict__
        conv2tif.main(cls.params)
    
    @classmethod
    def tearDownClass(cls):
        common.remove_tifs(cls.params[OBS_DIR])
        
    def setUp(self):
        self.params = config.get_config_params(TEST_CONF_GAMMA, step='prepifg')
        self.dummy_dir = '/i/should/not/exist'
        crop_opts = config._crop_opts(self.params)
        self.min_extents, self.n_cols, self.n_rows = \
            config._get_fullres_info(self.params[IFG_FILE_LIST], self.params[OUT_DIR], crop_opts)
        if os.path.exists(self.dummy_dir):
            raise IOError("{dummy_dir} needs to not exist for test purposes.")
    
    def test_validate_tifs_exist(self):
        validate_tifs_exist(self.params[IFG_FILE_LIST], self.params[OUT_DIR])
        with pytest.raises(ConfigException):
            validate_tifs_exist(self.params[IFG_FILE_LIST], self.dummy_dir)

    def test_validate_ifg_crop_coordinates(self):
        self.params[IFG_CROP_OPT] = 3
        self.params[IFG_XFIRST] = 150.0
        with pytest.raises(ConfigException):
            validate_crop_parameters(self.min_extents, self.params)
        self.params[IFG_XFIRST] = 150.92

        self.params[IFG_XLAST] = 150.95
        with pytest.raises(ConfigException):
            validate_crop_parameters(self.min_extents, self.params)
        self.params[IFG_XLAST] = 150.94

        self.params[IFG_YFIRST] = -34.16
        with pytest.raises(ConfigException):
            validate_crop_parameters(self.min_extents, self.params)
        self.params[IFG_YFIRST] = -34.18

        self.params[IFG_YLAST] = -34.24
        with pytest.raises(ConfigException):
            validate_crop_parameters(self.min_extents, self.params)

    def test_validate_ifg_multilook(self):
        self.params[IFG_LKSX] = 48
        self.params[IFG_LKSY] = 73
        with pytest.raises(ConfigException):
            validate_multilook_parameters(self.n_cols, self.n_rows,
                                          IFG_LKSX, IFG_LKSY, self.params)
    

class TestConfigValidationWithPrepifgGeotiffs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pyrate import conv2tif, prepifg
        cls.params = Configuration(TEST_CONF_GAMMA).__dict__
        conv2tif.main(cls.params)
        prepifg.main(cls.params)
    
    @classmethod
    def tearDownClass(cls):
        common.remove_tifs(cls.params[OBS_DIR])
        common.remove_tifs(cls.params[OUT_DIR])
        
    def setUp(self):
        self.params = config.get_config_params(TEST_CONF_GAMMA,
                                               'process')
        self.dummy_dir = '/i/should/not/exist'
        crop_opts = config._crop_opts(self.params)
        self.extents, self.n_cols, self.n_rows, _ = \
            _get_prepifg_info(self.params[IFG_FILE_LIST], self.params[OUT_DIR], self.params)
        self.n_epochs, self.max_span = \
            _get_temporal_info(self.params[IFG_FILE_LIST], self.params[OUT_DIR])
        if os.path.exists(self.dummy_dir):
            raise IOError("{dummy_dir} needs to not exist for test purposes.")

    def test_validate_prepifg_tifs_exist(self):
        validate_prepifg_tifs_exist(self.params[IFG_FILE_LIST], self.params[OUT_DIR], self.params)
        with pytest.raises(ConfigException):
            self.params[IFG_LKSX] = 100
            validate_prepifg_tifs_exist(self.params[IFG_FILE_LIST], self.params[OUT_DIR], self.params)

    def test_validate_epoch_thresholds(self):
        with pytest.raises(ConfigException):
            validate_minimum_epochs(self.n_epochs, 50)
        
        validate_epoch_thresholds(self.n_epochs, self.params)
        self.params[LR_PTHRESH] = 20
        with pytest.raises(ConfigException):
            validate_epoch_thresholds(self.n_epochs, self.params)
        self.params[LR_PTHRESH] = 5
        
        self.params[SLPF_CUTOFF] = 1000
        with pytest.raises(ConfigException):
            validate_epoch_cutoff(self.max_span, SLPF_CUTOFF, self.params)
        


    def test_validate_search_windows(self):
        self.params[REF_CHIP_SIZE] = 21
        self.params[REFNX] = 2
        self.params[REFNY] = 3
        validate_reference_pixel_search_windows(self.n_cols, self.n_rows, self.params)

        self.params[REFNX] = 3
        with pytest.raises(ConfigException):
            validate_reference_pixel_search_windows(self.n_cols, self.n_rows, self.params)

        self.params[REFNX] = 2
        self.params[REFNY] = 4
        with pytest.raises(ConfigException):
            validate_reference_pixel_search_windows(self.n_cols, self.n_rows, self.params)
    
    def test_validate_multilook_parameters(self):
        self.params[ORBITAL_FIT_LOOKS_X] = 48
        self.params[ORBITAL_FIT_LOOKS_Y] = 73
        with pytest.raises(ConfigException):
            validate_multilook_parameters(self.n_cols, self.n_rows,
                ORBITAL_FIT_LOOKS_X, ORBITAL_FIT_LOOKS_Y, self.params)       

    def test_validate_slpf_cutoff(self):
        self.params[SLPF_CUTOFF] = 9999
        with pytest.raises(ConfigException):
            validate_slpf_cutoff(self.extents, self.params)


class ConfigTest(unittest.TestCase):

    @staticmethod
    def test_read_param_file():
        params = config.get_config_params(TEST_CONF_ROIPAC)
        for k in params.keys():
            assert k and len(k) > 1
            assert params[k] != ''
            assert not k.endswith(":")  # are the colons removed?

    @staticmethod
    def test_read_param_file_missing_option():
        # ensure the parser can handle missing option fields
        conf_path = join(SML_TEST_CONF, 'pyrate1.conf')
        params = config.get_config_params(conf_path)

        assert params[REFX] == -1
        assert params[REFY] == -1

    @staticmethod
    def test_read_param_file_missing_value():
        # ensure the parser can handle blank option values
        conf_path = join(SML_TEST_CONF, 'pyrate2.conf')
        params = config.get_config_params(conf_path)

        assert params[REFX] == -1
        assert params[REFY] == -1

    @staticmethod
    def test_parse_namelist():
        nl = join(SML_TEST_TIF, 'ifms_17')
        result = list(config.parse_namelist(nl))
        assert len(result) == 17
        files = ["geo_060619-061002_unw.tif", "geo_060828-061211_unw.tif",
                    "geo_061002-070430_unw.tif", "geo_070115-070917_unw.tif",
                    "geo_070219-070604_unw.tif"]

        for path in files:
            assert path in result


class ConfigWriteTest(unittest.TestCase):

    def test_write_config_file(self):
        params = config.get_config_params(TEST_CONF_GAMMA)
        temp_config = tempfile.mktemp(suffix='.conf')
        config.write_config_file(params, temp_config)
        self.assertTrue(os.path.exists(temp_config))
        os.remove(temp_config)

    def test_new_config_file_and_original_match(self):
        params = config.get_config_params(TEST_CONF_GAMMA)
        temp_config = tempfile.mktemp(suffix='.conf')
        config.write_config_file(params, temp_config)
        new_params = config.get_config_params(temp_config)
        self.maxDiff = None
        self.assertDictEqual(params, new_params)
        os.remove(temp_config)


class ConfigAPSParametersTest(unittest.TestCase):
    def setUp(self):
        self.conf_path = TEST_CONF_ROIPAC
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
        self.ifgListFile = os.path.join(common.SML_TEST_GAMMA, 'ifms_17')

    def tearDown(self):
        shutil.rmtree(self.base_dir)

    def make_input_files(self, inc='', ele=''):
        with open(self.conf_file, 'w') as conf:
            conf.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, common.SML_TEST_GAMMA))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '1'))
            conf.write('{}: {}\n'.format(
                DEM_HEADER_FILE, os.path.join(
                    common.SML_TEST_GAMMA, '20060619_utm_dem.par')))
            conf.write('{}: {}\n'.format(IFG_LKSX, '1'))
            conf.write('{}: {}\n'.format(IFG_LKSY, '1'))
            conf.write('{}: {}\n'.format(IFG_CROP_OPT, '1'))
            conf.write('{}: {}\n'.format(NO_DATA_AVERAGING_THRESHOLD, '0.5'))
            conf.write('{}: {}\n'.format(SLC_DIR, ''))
            conf.write('{}: {}\n'.format(HDR_FILE_LIST, common.SML_TEST_GAMMA_HEADER_LIST))
            conf.write('{}: {}\n'.format(DEM_FILE, common.SML_TEST_DEM_GAMMA))
            conf.write('{}: {}\n'.format(APS_INCIDENCE_MAP, inc))
            conf.write('{}: {}\n'.format(APS_ELEVATION_MAP, ele))
            conf.write('{}: {}\n'.format(APS_CORRECTION, '1'))
            conf.write('{}: {}\n'.format(APS_METHOD, '2'))

    def test_inc_vs_ele_maps_inc_provided(self):
        self.make_input_files(inc=common.SML_TEST_INCIDENCE)
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
        self.make_input_files(ele=common.SML_TEST_ELEVATION)
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


if __name__ == "__main__":
    unittest.main()
