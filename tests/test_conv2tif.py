import os
import unittest
import pytest
import glob
import copy

import pyrate.core.config as cf
from pyrate import conv2tif, prepifg
from tests import common


class ConvertToGeotiffTests(unittest.TestCase):
    """
    Test basic functionality of converttogeotif - they get placed in the
    right directory and there's the right amount.
    """
    @classmethod
    def setUpClass(cls):
        cls.gamma_params = cf.get_config_params(common.TEST_CONF_GAMMA) 
        cls.roipac_params = cf.get_config_params(common.TEST_CONF_ROIPAC)

    def test_dem_and_incidence_not_converted(self):
        gp_copy = copy.deepcopy(self.gamma_params)
        gp_copy[cf.DEM_FILE] = None
        gp_copy[cf.APS_INCIDENCE_MAP] = None
        conv2tif.main(gp_copy)
        inc_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*inc.tif'))
        self.assertEqual(len(inc_tif), 0)
        dem_tif = glob.glob(os.path.join(gp_copy[cf.OBS_DIR], '*dem.tif'))
        self.assertEqual(len(dem_tif), 0)

    def test_tifs_placed_in_obs_dir(self):
        # Test no tifs in obs dir
        tifs = glob.glob(os.path.join(self.gamma_params[cf.OBS_DIR], '*.tif'))
        self.assertEqual(len(tifs), 0)
        # Test tifs in obs dir
        conv2tif.main(self.gamma_params)
        tifs = glob.glob(os.path.join(self.gamma_params[cf.OBS_DIR], '*.tif'))
        self.assertEqual(len(tifs), 19)

    def test_num_gamma_tifs_equals_num_unws(self):
        gtifs = conv2tif.main(self.gamma_params)
        # 17 unws + incidence + dem
        self.assertEqual(len(gtifs), 19)
        
    def test_num_roipac_tifs_equals_num_unws(self):
        gtifs = conv2tif.main(self.roipac_params)
        # 17 unws + dem
        self.assertEqual(len(gtifs), 18)

    def test_conversion_not_recomputed(self):
        """
        If a gtif already exists, conv2tif will return None instead
        of path to the gtif.
        """
        conv2tif.main(self.gamma_params)
        gtifs = conv2tif.main(self.gamma_params)
        self.assertTrue(all([gt is None for gt in gtifs]))
        
    def teardown_method(self, method):
        common.remove_tifs(self.gamma_params[cf.OBS_DIR])
        common.remove_tifs(self.roipac_params[cf.OBS_DIR])
 
class PrepifgConversionTests(unittest.TestCase):
    """
    Test that prepifg fails if there are no converted IFGs in the observations
    directory, and that it succeeds if there are.
    """
    @classmethod
    def setUpClass(cls):
        cls.params = cf.get_config_params(common.TEST_CONF_GAMMA) 

    def test_no_tifs_exits(self):
        with pytest.raises(SystemExit):
            prepifg.main(self.params)

    def test_tifs_succeeds(self):
        conv2tif.main(self.params)
        prepifg.main(self.params)
        
    def teardown_method(self, method):
        common.remove_tifs(self.params[cf.OBS_DIR])

