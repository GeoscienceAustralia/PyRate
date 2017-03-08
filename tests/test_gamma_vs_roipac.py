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
This Python module contains tests that compare GAMMA and ROI_PAC
functionality in PyRate.
"""
import glob
import os
import re
import shutil
import sys
import tempfile
import unittest

import numpy as np
from tests import common

from pyrate import config as cf
from pyrate import ifgconstants as ifc
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
    SLC_DIR,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP
    )
from pyrate.scripts import run_prepifg
from tests.common import SYD_TEST_DIR
from tests.common import sydney_data_setup

DUMMY_SECTION_NAME = 'pyrate'


class TestGammaVsRoipacEquality(unittest.TestCase):

    SYDNEY_GAMMA_TEST = os.path.join(SYD_TEST_DIR, 'gamma_sydney_test')

    @classmethod
    def setUpClass(cls):
        cls.gamma_base_dir = tempfile.mkdtemp()
        cls.gamma_conffile = os.path.join(cls.gamma_base_dir, 'gamma_test.conf')
        cls.gamma_ifgListFile = os.path.join(cls.gamma_base_dir,
                                             'gamma_ifg.list')

        cls.roipac_base_dir = tempfile.mkdtemp()
        cls.roipac_conffile = os.path.join(cls.roipac_base_dir,
                                           'roipac_test.conf')
        cls.roipac_ifgListFile = os.path.join(cls.roipac_base_dir,
                                              'roipac_ifg.list')

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
            conf.write('{}: {}\n'.format(SLC_DIR, ''))
            conf.write('{}: {}\n'.format(DEM_FILE, common.SYD_TEST_DEM_GAMMA))
            conf.write('{}: {}\n'.format(APS_INCIDENCE_MAP,
                                         common.SYD_TEST_INCIDENCE))
            conf.write('{}: {}\n'.format(APS_ELEVATION_MAP, ''))
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
        sys.argv = ['pyrate', 'prepifg', conf_file]

        base_ifg_paths, dest_paths, params = cf.get_ifg_paths(conf_file)
        dest_base_ifgs = [os.path.join(
            params[cf.OUT_DIR], os.path.basename(q).split('.')[0] + '_' +
            os.path.basename(q).split('.')[1] + '.tif') 
            for q in base_ifg_paths]
        sys.argv = ['pyrate', 'prepifg', conf_file]
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
            conf.write('{}: {}\n'.format(DEM_FILE, common.SYD_TEST_DEM_ROIPAC))
            conf.write('{}: {}\n'.format(APS_INCIDENCE_MAP, ''))
            conf.write('{}: {}\n'.format(APS_ELEVATION_MAP, ''))
        with open(self.ifgListFile, 'w') as ifgl:
            ifgl.write('\n'.join(data))

    def test_cmd_ifg_no_roipac_files_created_roipac(self):
        self.dataPaths = common.sydney_data_roipac_unws()
        base_exp = common.sydney_ifg_file_list()
        self.expPaths = [os.path.join(self.roipac_base_dir, os.path.basename(i))
                         for i in base_exp]
        self.luigi = '0'
        self.confFile = self.roipac_conffile
        self.ifgListFile = self.roipac_ifgListFile
        self.base_dir = self.roipac_base_dir
        self.check_roipac()

    def check_roipac(self):
        self.make_roipac_input_files(self.dataPaths, 'WGS84')
        sys.argv = ['pyrate', 'prepifg', self.confFile]
        run_prepifg.main()
        for path in self.expPaths:
            self.assertTrue(os.path.exists(path),
                            '{} does not exist'.format(path))

    def test_equality_of_luigi_and_no_luigi_phase_data(self):
        gamma_PTN = re.compile(r'\d{8}')
        gamma_files = []
        for i in glob.glob(os.path.join(self.gamma_base_dir,
                                        "*.tif")):
            if len(gamma_PTN.findall(i)) == 2:
                gamma_files.append(i)
        all_gamma_ifgs = sydney_data_setup(gamma_files)
        all_roipac_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.roipac_base_dir, "geo*.tif")))
        c = 0
        for c, (r, g) in enumerate(zip(all_roipac_ifgs, all_gamma_ifgs)):
            np.testing.assert_array_equal(r.phase_data, g.phase_data)
        self.assertEquals(c+1, len(all_roipac_ifgs))
        self.assertEquals(c+1, len(all_gamma_ifgs))

    def test_equality_of_meta_data(self):
        gamma_PTN = re.compile(r'\d{8}')
        gamma_files = []
        for i in glob.glob(os.path.join(self.gamma_base_dir,
                                        "*.tif")):
            if len(gamma_PTN.findall(i)) == 2:
                gamma_files.append(i)
        all_gamma_ifgs = sydney_data_setup(gamma_files)
        all_roipac_ifgs = sydney_data_setup(
            glob.glob(os.path.join(self.roipac_base_dir, "geo*.tif")))
        c = 0
        for c, (i, j) in enumerate(zip(all_gamma_ifgs, all_roipac_ifgs)):
            mdi = i.meta_data
            mdj = j.meta_data
            for k in mdi:  # all key values equal
                if k == 'INCIDENCE_DEGREES':
                    pass # incidence angle not implemented for roipac
                elif is_number(mdi[k]):
                    self.assertAlmostEqual(
                        float(mdj[k]), float(mdi[k]), places=6)
                elif mdi[k] == 'ROIPAC' or 'GAMMA':
                    pass # INSAR_PROCESSOR can not be equal
                else:
                    self.assertEquals(mdj[k], mdi[k])
            if i.data_path.__contains__(
                    '_{looks}rlks_{crop}cr'.format(looks=1, crop=1)):
                # these are multilooked tifs
                # test that DATA_STEP is MULTILOOKED
                self.assertEqual(mdi[ifc.DATA_TYPE], ifc.MULTILOOKED)
                self.assertEqual(mdj[ifc.DATA_TYPE], ifc.MULTILOOKED)
            else:
                self.assertEqual(mdi[ifc.DATA_TYPE], ifc.ORIG)
                self.assertEqual(mdj[ifc.DATA_TYPE], ifc.ORIG)

        self.assertEquals(c + 1, len(all_gamma_ifgs))


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    unittest.main()
