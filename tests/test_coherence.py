import os
import unittest
import pytest
import glob
import copy

import pyrate.core.config as cf
from pyrate import converttogtif, prepifg
from tests import common

class CoherenceMaskingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = cf.get_config_params(common.TEST_CONF_GAMMA)
        cls.params[cf.COH_MASK] = 1
        cls.params[cf.COH_FILE_DIR] = common.SML_TEST_COH_DIR
        cls.params[cf.COH_FILE_LIST] = common.SML_TEST_COH_LIST

    def teardown_method(self, method):
        common.remove_tifs(self.params[cf.OBS_DIR])
        common.remove_tifs(self.params[cf.COH_FILE_DIR])

    def test_find_coherence_files_for_ifgs(self):
        ifg_paths = [os.path.join(self.params[cf.OBS_DIR], ifg) 
                     for ifg in cf.parse_namelist(self.params[cf.IFG_FILE_LIST])]
        coherence_files = [cf.coherence_paths_for(path, self.params, tif=False)
                           for path in ifg_paths]
        self.assertEqual(len(coherence_files), 17)
    
    @unittest.skip("Skip this test until actual coherence test files provided")
    def test_coherence_files_converted(self):
        self.params[cf.COH_MASK] = 1
        gtiff_paths = converttogtif.main(self.params)
        coh_files = [path for path in gtiff_paths if '_cc.tif' in path]
        self.assertEqual(len(coh_files), 17)
        self.params[cf.COH_MASK] = 0

    def test_coherence_files_not_converted(self):
        self.params[cf.COH_MASK] = 0
        gtiff_paths = converttogtif.main(self.params)
        coh_files = [path for path in gtiff_paths if '_cc.tif' in path]
        self.assertEqual(len(coh_files), 0)
        self.params[cf.COH_MASK] = 1
        
    # TOD: Testing actual coherence masking
    # Need coherence file test data
    # Need pre-masked IFGs
    # Compare that output of coherence masking matches pre-masked IFGs   

