__author__ = 'sudipta'

import unittest
import os
import shutil
import uuid
import glob
import numpy as np
import sys

from pyrate.tests.common import SYD_TEST_DIR, TEMPDIR
from pyrate.tests import common
from pyrate.scripts import run_prepifg, run_pyrate
from pyrate.tests.common import sydney_data_setup
from pyrate import config as cf


from pyrate.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    LUIGI,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    NO_DATA_AVERAGING_THRESHOLD,
    INPUT_IFG_PROJECTION,
    ROIPAC_RESOURCE_HEADER,
    )

DUMMY_SECTION_NAME = 'pyrate'


class TestGammaVsRoipacEquality(unittest.TestCase):

    SYDNEY_GAMMA_TEST = os.path.join(SYD_TEST_DIR, 'gamma_sydney_test')

    @classmethod
    def setUpClass(cls):
        gamma_dir= uuid.uuid4().hex
        cls.gamma_conffile = os.path.join(TEMPDIR,
            '{}/gamma_test.conf'.format(gamma_dir))
        cls.gamma_ifgListFile = os.path.join(TEMPDIR,
            '{}/gamma_ifg.list'.format(gamma_dir))
        cls.gamma_base_dir = os.path.dirname(cls.gamma_conffile)
        common.mkdir_p(cls.gamma_base_dir)

        roipac_dir = uuid.uuid4().hex
        cls.roipac_conffile = os.path.join(
            TEMPDIR, '{}/roipac_test.conf'.format(roipac_dir))
        cls.roipac_ifgListFile = os.path.join(
            TEMPDIR, '{}/roipac_ifg.list'.format(roipac_dir))

        cls.roipac_base_dir = os.path.dirname(cls.roipac_conffile)

        common.mkdir_p(cls.roipac_base_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.gamma_base_dir)
        shutil.rmtree(cls.roipac_base_dir)

    def make_gamma_input_files(self, data):
        with open(self.conf_file, 'w') as conf:
            conf.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '1'))
            conf.write('{}: {}\n'.format(LUIGI, self.LUIGI))
            conf.write('{}: {}\n'.format(
                DEM_HEADER_FILE, os.path.join(
                    self.SYDNEY_GAMMA_TEST, '20060619_utm_dem.par')))
            conf.write('{}: {}\n'.format(IFG_LKSX, '1'))
            conf.write('{}: {}\n'.format(IFG_LKSY, '1'))
            conf.write('{}: {}\n'.format(IFG_CROP_OPT, '1'))
            conf.write('{}: {}\n'.format(NO_DATA_AVERAGING_THRESHOLD, '0.5'))
        with open(self.ifgListFile, 'w') as ifgl:
            ifgl.write('\n'.join(data))

    def test_cmd_ifg_no_gamma_files_created(self):
        self.LUIGI = '0'  # luigi or no luigi
        self.conf_file = self.gamma_conffile
        self.base_dir = self.gamma_base_dir
        self.ifgListFile = self.gamma_ifgListFile
        self.check_gamma(self.gamma_conffile)

    def check_gamma(self, conf_file):
        data_paths = glob.glob(
            os.path.join(self.SYDNEY_GAMMA_TEST, "*_utm.unw"))

        self.make_gamma_input_files(data_paths)
        sys.argv = ['run_pyrate.py', conf_file]

        base_ifg_paths, dest_paths, params = run_pyrate.get_ifg_paths()
        dest_base_ifgs = [os.path.join(
            params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '.tif')
            for q in base_ifg_paths]
        sys.argv = ['run_prepifg.py', conf_file]
        run_prepifg.main()

        for p, q in zip(dest_base_ifgs, dest_paths):
            self.assertTrue(os.path.exists(p),
                            '{} does not exist'.format(p))
            self.assertTrue(os.path.exists(q),
                            '{} does not exist'.format(q))

    def make_roipac_input_files(self, data, projection):
        with open(self.confFile, 'w') as conf:
            conf.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            conf.write('{}: {}\n'.format(INPUT_IFG_PROJECTION, projection))
            conf.write('{}: {}\n'.format(NO_DATA_VALUE, '0.0'))
            conf.write('{}: {}\n'.format(OBS_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(OUT_DIR, self.base_dir))
            conf.write('{}: {}\n'.format(IFG_FILE_LIST, self.ifgListFile))
            conf.write('{}: {}\n'.format(PROCESSOR, '0'))
            conf.write('{}: {}\n'.format(LUIGI, self.luigi))
            conf.write('{}: {}\n'.format(ROIPAC_RESOURCE_HEADER,
                                         common.SYD_TEST_DEM_HDR))
            conf.write('{}: {}\n'.format(IFG_LKSX, '1'))
            conf.write('{}: {}\n'.format(IFG_LKSY, '1'))
            conf.write('{}: {}\n'.format(IFG_CROP_OPT, '1'))
            conf.write('{}: {}\n'.format(NO_DATA_AVERAGING_THRESHOLD, '0.5'))
        with open(self.ifgListFile, 'w') as ifgl:
            ifgl.write('\n'.join(data))

    def test_cmd_ifg_no_roipac_files_created_roipac(self):
        self.dataPaths = common.sydney_data_setup_unw_file_list()
        base_exp = common.sydney_data_setup_ifg_file_list()
        self.expPaths = [os.path.join(self.roipac_base_dir, os.path.basename(i))
                         for i in base_exp]
        self.luigi = '0'
        self.confFile = self.roipac_conffile
        self.ifgListFile = self.roipac_ifgListFile
        self.base_dir = self.roipac_base_dir
        self.check_roipac()

    def check_roipac(self):
        self.make_roipac_input_files(self.dataPaths, 'WGS84')
        sys.argv = ['run_prepifg.py', self.confFile]
        run_prepifg.main()
        for path in self.expPaths:
            self.assertTrue(os.path.exists(path),
                            '{} does not exist'.format(path))

    def test_equality_of_luigi_and_no_luigi_phase_data(self):

        all_gamma_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.gamma_base_dir, "*.tif")))
        all_roipac_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.roipac_base_dir, "*.tif")))

        c = 0
        for c, (r, g) in enumerate(zip(all_roipac_ifgs, all_gamma_ifgs)):
            np.testing.assert_array_equal(r.phase_data, g.phase_data)
        self.assertEquals(c+1, len(all_roipac_ifgs))
        self.assertEquals(c+1, len(all_gamma_ifgs))

    def test_eqality_of_meta_data(self):
        all_gamma_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.gamma_base_dir, "*.tif")))
        all_roipac_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.roipac_base_dir, "*.tif")))
        c = 0
        for c, (i, j) in enumerate(zip(all_gamma_ifgs, all_roipac_ifgs)):
            mdi = i.meta_data
            mdj = j.meta_data
            for k in mdi:  # all key values equal
                if is_number(mdi[k]):
                    self.assertAlmostEqual(
                        float(mdj[k]), float(mdi[k]), places=6)
                else:
                    self.assertEquals(mdj[k], mdi[k])

        self.assertEquals(c + 1, len(all_gamma_ifgs))


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    unittest.main()
