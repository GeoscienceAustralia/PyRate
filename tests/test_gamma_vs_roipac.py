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
This Python module contains tests that compare GAMMA and ROI_PAC
functionality in PyRate.
"""
import glob
import os
import re
import shutil
import tempfile
import unittest

from . import common
from common import SML_TEST_DIR, small_data_setup, remove_tifs
from core import ifgconstants as ifc, config as cf
from core.config import (
    DEM_HEADER_FILE,
    NO_DATA_VALUE,
    OBS_DIR,
    IFG_FILE_LIST,
    PROCESSOR,
    OUT_DIR,
    IFG_LKSX,
    IFG_LKSY,
    IFG_CROP_OPT,
    NO_DATA_AVERAGING_THRESHOLD,
    INPUT_IFG_PROJECTION,
    SLC_DIR,
    SLC_FILE_LIST,
    DEM_FILE,
    APS_INCIDENCE_MAP,
    APS_ELEVATION_MAP,
)
from core.prepifg_helper import _is_number
from main import prepifg_handler, conv2tif_handler
from . import common

DUMMY_SECTION_NAME = "pyrate"


class TestGammaVsRoipacEquality(unittest.TestCase):
    """ """

    SMLNEY_GAMMA_TEST = os.path.join(SML_TEST_DIR, "gamma_obs")

    @classmethod
    def setUpClass(cls):
        """ """
        cls.gamma_base_dir = tempfile.mkdtemp()
        cls.gamma_conffile = os.path.join(cls.gamma_base_dir, "gamma_test.conf")
        cls.gamma_ifgListFile = os.path.join(cls.gamma_base_dir, "gamma_ifg.list")

        cls.roipac_base_dir = tempfile.mkdtemp()
        cls.roipac_conffile = os.path.join(cls.roipac_base_dir, "roipac_test.conf")
        cls.roipac_ifgListFile = os.path.join(cls.roipac_base_dir, "roipac_ifg.list")

    @classmethod
    def tearDownClass(cls):
        """ """
        shutil.rmtree(cls.gamma_base_dir)
        shutil.rmtree(cls.roipac_base_dir)
        remove_tifs(cls.SMLNEY_GAMMA_TEST)
        remove_tifs(common.SML_TEST_OBS)

    def make_gamma_input_files(self, data):
        """

        Args:
          data: 

        Returns:

        """
        with open(self.conf_file, "w") as conf:
            conf.write("[{}]\n".format(DUMMY_SECTION_NAME))
            conf.write("{}: {}\n".format(NO_DATA_VALUE, "0.0"))
            conf.write("{}: {}\n".format(OBS_DIR, self.SMLNEY_GAMMA_TEST))
            conf.write("{}: {}\n".format(OUT_DIR, self.base_dir))
            conf.write("{}: {}\n".format(IFG_FILE_LIST, self.ifgListFile))
            conf.write("{}: {}\n".format(PROCESSOR, "1"))
            conf.write("{}: {}\n".format(DEM_HEADER_FILE, os.path.join(self.SMLNEY_GAMMA_TEST, "20060619_utm_dem.par")))
            conf.write("{}: {}\n".format(IFG_LKSX, "1"))
            conf.write("{}: {}\n".format(IFG_LKSY, "1"))
            conf.write("{}: {}\n".format(IFG_CROP_OPT, "1"))
            conf.write("{}: {}\n".format(NO_DATA_AVERAGING_THRESHOLD, "0.5"))
            conf.write("{}: {}\n".format(SLC_DIR, ""))
            conf.write("{}: {}\n".format(SLC_FILE_LIST, common.SML_TEST_GAMMA_HEADER_LIST))
            conf.write("{}: {}\n".format(DEM_FILE, common.SML_TEST_DEM_GAMMA))
            conf.write("{}: {}\n".format(APS_INCIDENCE_MAP, common.SML_TEST_INCIDENCE))
            conf.write("{}: {}\n".format(APS_ELEVATION_MAP, ""))
            conf.write("{}: {}\n".format("refny", ""))
            conf.write("{}: {}\n".format("noDataValue", ""))
            conf.write("{}: {}\n".format("tlpfcutoff", ""))
            conf.write("{}: {}\n".format("axsig", ""))
            conf.write("{}: {}\n".format("slpfmethod", ""))
            conf.write("{}: {}\n".format("orbfit", ""))
            conf.write("{}: {}\n".format("slcfilelist", ""))
            conf.write("{}: {}\n".format("obsdir", ""))
            conf.write("{}: {}\n".format("processes", ""))
            conf.write("{}: {}\n".format("orbfitlksy", ""))
            conf.write("{}: {}\n".format("nan_conversion", ""))
            conf.write("{}: {}\n".format("refchipsize", ""))
            conf.write("{}: {}\n".format("slpfcutoff", ""))
            conf.write("{}: {}\n".format("cohfilelist", ""))
            conf.write("{}: {}\n".format("demHeaderFile", ""))
            conf.write("{}: {}\n".format("cohfiledir", ""))
            conf.write("{}: {}\n".format("processor", ""))
            conf.write("{}: {}\n".format("fgxfirst", ""))
            conf.write("{}: {}\n".format("refnx", ""))
            conf.write("{}: {}\n".format("slcFileDir", ""))
            conf.write("{}: {}\n".format("orbfitlksx", ""))
            conf.write("{}: {}\n".format("tlpfpthr", ""))
            conf.write("{}: {}\n".format("slpnanfill_method", ""))
            conf.write("{}: {}\n".format("cohmask", ""))
            conf.write("{}: {}\n".format("ifgylast", ""))
            conf.write("{}: {}\n".format("smorder", ""))
            conf.write("{}: {}\n".format("ifglksx", ""))
            conf.write("{}: {}\n".format("ifgxlast", ""))
            conf.write("{}: {}\n".format("apsest", ""))
            conf.write("{}: {}\n".format("tlpfmethod", ""))
            conf.write("{}: {}\n".format("cohthresh", ""))
            conf.write("{}: {}\n".format("tscal", ""))
            conf.write("{}: {}\n".format("pthr", ""))
            conf.write("{}: {}\n".format("noDataAveragingThreshold", ""))
            conf.write("{}: {}\n".format("tsmethod", ""))
            conf.write("{}: {}\n".format("smfactor", ""))
            conf.write("{}: {}\n".format("parallel", ""))
            conf.write("{}: {}\n".format("orbfitdegrees", ""))
            conf.write("{}: {}\n".format("ifgfilelist", ""))
            conf.write("{}: {}\n".format("ifgcropopt", ""))
            conf.write("{}: {}\n".format("outdir", ""))
            conf.write("{}: {}\n".format("refx", ""))
            conf.write("{}: {}\n".format("orbfitmethod", ""))
            conf.write("{}: {}\n".format("refy", ""))
            conf.write("{}: {}\n".format("ts_pthr", ""))
            conf.write("{}: {}\n".format("nsig", ""))
            conf.write("{}: {}\n".format("ifgyfirst", ""))
            conf.write("{}: {}\n".format("refminfrac", ""))
            conf.write("{}: {}\n".format("refest", ""))
            conf.write("{}: {}\n".format("demfile", ""))
            conf.write("{}: {}\n".format("slpnanfill", ""))
            conf.write("{}: {}\n".format("ifglksy", ""))
            conf.write("{}: {}\n".format("slpforder", ""))

        with open(self.ifgListFile, "w") as ifgl:
            ifgl.write("\n".join(data))

    def test_cmd_ifg_no_gamma_files_created(self):
        self.conf_file = self.gamma_conffile
        self.base_dir = self.gamma_base_dir
        self.ifgListFile = self.gamma_ifgListFile
        self.check_gamma(self.gamma_conffile)

    def check_gamma(self, conf_file):
        """

        Args:
          conf_file: 

        Returns:

        """
        data_paths = glob.glob(os.path.join(self.SMLNEY_GAMMA_TEST, "*_utm.unw"))

        self.make_gamma_input_files(data_paths)

        base_ifg_paths, dest_paths, params = cf.get_ifg_paths(conf_file)
        dest_base_ifgs = [
            os.path.join(params[cf.OBS_DIR],
                         os.path.basename(q).split(".")[0] + "_" + os.path.basename(q).split(".")[1] + ".tif")
            for q in base_ifg_paths
        ]

        conv2tif_handler(conf_file)
        prepifg_handler(conf_file)

        for p, q in zip(dest_base_ifgs, dest_paths):
            self.assertTrue(os.path.exists(p), "{} does not exist".format(p))
            self.assertTrue(os.path.exists(q), "{} does not exist".format(q))

    def make_roipac_input_files(self, data, projection):
        """

        Args:
          data: param projection:
          projection: 

        Returns:

        """
        with open(self.confFile, "w") as conf:
            conf.write("[{}]\n".format(DUMMY_SECTION_NAME))
            conf.write("{}: {}\n".format(INPUT_IFG_PROJECTION, projection))
            conf.write("{}: {}\n".format(NO_DATA_VALUE, "0.0"))
            conf.write("{}: {}\n".format(OBS_DIR, common.SML_TEST_OBS))
            conf.write("{}: {}\n".format(OUT_DIR, self.base_dir))
            conf.write("{}: {}\n".format(IFG_FILE_LIST, self.ifgListFile))
            conf.write("{}: {}\n".format(PROCESSOR, "0"))
            conf.write("{}: {}\n".format(DEM_HEADER_FILE, common.SML_TEST_DEM_HDR))
            conf.write("{}: {}\n".format(IFG_LKSX, "1"))
            conf.write("{}: {}\n".format(IFG_LKSY, "1"))
            conf.write("{}: {}\n".format(IFG_CROP_OPT, "1"))
            conf.write("{}: {}\n".format(NO_DATA_AVERAGING_THRESHOLD, "0.5"))
            conf.write("{}: {}\n".format(DEM_FILE, common.SML_TEST_DEM_ROIPAC))
            conf.write("{}: {}\n".format(APS_INCIDENCE_MAP, ""))
            conf.write("{}: {}\n".format(APS_ELEVATION_MAP, ""))
            conf.write("{}: {}\n".format("ifgfilelist", ""))
            conf.write("{}: {}\n".format("slpforder", ""))
            conf.write("{}: {}\n".format("demHeaderFile", ""))
            conf.write("{}: {}\n".format("ifglksx", ""))
            conf.write("{}: {}\n".format("cohthresh", ""))
            conf.write("{}: {}\n".format("tlpfmethod", ""))
            conf.write("{}: {}\n".format("slcFileDir", ""))
            conf.write("{}: {}\n".format("ifgcropopt", ""))
            conf.write("{}: {}\n".format("ifglksy", ""))
            conf.write("{}: {}\n".format("demfile", ""))
            conf.write("{}: {}\n".format("cohfilelist", ""))
            conf.write("{}: {}\n".format("orbfit", ""))
            conf.write("{}: {}\n".format("ts_pthr", ""))
            conf.write("{}: {}\n".format("slpnanfill_method", ""))
            conf.write("{}: {}\n".format("tscal", ""))
            conf.write("{}: {}\n".format("refy", ""))
            conf.write("{}: {}\n".format("processor", ""))
            conf.write("{}: {}\n".format("slpnanfill", ""))
            conf.write("{}: {}\n".format("orbfitdegrees", ""))
            conf.write("{}: {}\n".format("ifgxfirst", ""))
            conf.write("{}: {}\n".format("cohfiledir", ""))
            conf.write("{}: {}\n".format("pthr", ""))
            conf.write("{}: {}\n".format("refnx", ""))
            conf.write("{}: {}\n".format("tlpfcutoff", ""))
            conf.write("{}: {}\n".format("processes", ""))
            conf.write("{}: {}\n".format("outdir", ""))
            conf.write("{}: {}\n".format("refny", ""))
            conf.write("{}: {}\n".format("nsig", ""))
            conf.write("{}: {}\n".format("tsmethod", ""))
            conf.write("{}: {}\n".format("ifgxlast", ""))
            conf.write("{}: {}\n".format("refchipsize", ""))
            conf.write("{}: {}\n".format("noDataValue", ""))
            conf.write("{}: {}\n".format("noDataAveragingThreshold", ""))
            conf.write("{}: {}\n".format("slcfilelist", ""))
            conf.write("{}: {}\n".format("ifgyfirst", ""))
            conf.write("{}: {}\n".format("cohmask", ""))
            conf.write("{}: {}\n".format("orbfitmethod", ""))
            conf.write("{}: {}\n".format("obsdir", ""))
            conf.write("{}: {}\n".format("refminfrac", ""))
            conf.write("{}: {}\n".format("maxsig", ""))
            conf.write("{}: {}\n".format("smorder", ""))
            conf.write("{}: {}\n".format("apsest", ""))
            conf.write("{}: {}\n".format("orbfitlksy", ""))
            conf.write("{}: {}\n".format("slpfmethod", ""))
            conf.write("{}: {}\n".format("ifgylast", ""))
            conf.write("{}: {}\n".format("smfactor", ""))
            conf.write("{}: {}\n".format("tlpfpthr", ""))
            conf.write("{}: {}\n".format("refest", ""))
            conf.write("{}: {}\n".format("nan_conversion", ""))
            conf.write("{}: {}\n".format("slpfcutoff", ""))
            conf.write("{}: {}\n".format("orbfitlksx", ""))
            conf.write("{}: {}\n".format("parallel", ""))
            conf.write("{}: {}\n".format("refx", ""))


        with open(self.ifgListFile, "w") as ifgl:
            ifgl.write("\n".join(data))

    def test_cmd_ifg_no_roipac_files_created_roipac(self):
        self.dataPaths = common.small_data_roipac_unws()
        base_exp = common.small_ifg_file_list()
        self.expPaths = [os.path.join(common.SML_TEST_OBS, os.path.basename(i)) for i in base_exp]
        self.confFile = self.roipac_conffile
        self.ifgListFile = self.roipac_ifgListFile
        self.base_dir = self.roipac_base_dir
        self.check_roipac()

    def check_roipac(self):
        """ """
        self.make_roipac_input_files(self.dataPaths, "WGS84")

        conv2tif_handler(self.confFile)
        prepifg_handler(self.confFile)
        for path in self.expPaths:
            self.assertTrue(os.path.exists(path), "{} does not exist".format(path))

    def test_equality_of_meta_data(self):
        """ """
        gamma_PTN = re.compile(r"\d{8}")
        gamma_files = []
        for i in glob.glob(os.path.join(self.gamma_base_dir, "*.tif")):
            if len(gamma_PTN.findall(i)) == 2:
                gamma_files.append(i)
        all_gamma_ifgs = small_data_setup(gamma_files)
        all_roipac_ifgs = small_data_setup(glob.glob(os.path.join(self.roipac_base_dir, "geo*.tif")))
        c = 0
        for c, (i, j) in enumerate(zip(all_gamma_ifgs, all_roipac_ifgs)):
            mdi = i.meta_data
            mdj = j.meta_data
            for k in mdi:  # all key values equal
                if k == "INCIDENCE_DEGREES":
                    pass  # incidence angle not implemented for roipac
                elif _is_number(mdi[k]):
                    self.assertAlmostEqual(float(mdj[k]), float(mdi[k]), places=6)
                elif mdi[k] == "ROIPAC" or "GAMMA":
                    pass  # INSAR_PROCESSOR can not be equal
                else:
                    self.assertEqual(mdj[k], mdi[k])
            if i.data_path.__contains__("_{looks}rlks_{crop}cr".format(looks=1, crop=1)):
                # these are multilooked tifs
                # test that DATA_STEP is MULTILOOKED
                self.assertEqual(mdi[ifc.DATA_TYPE], ifc.MULTILOOKED)
                self.assertEqual(mdj[ifc.DATA_TYPE], ifc.MULTILOOKED)
            else:
                self.assertEqual(mdi[ifc.DATA_TYPE], ifc.ORIG)
                self.assertEqual(mdj[ifc.DATA_TYPE], ifc.ORIG)

        self.assertEqual(c + 1, len(all_gamma_ifgs))


if __name__ == "__main__":
    unittest.main()
