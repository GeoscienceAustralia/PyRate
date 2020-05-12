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
# coding: utf-8
"""
This Python module contains tests for the ref_phs_est.py PyRate module.
"""
import glob
import os
import shutil
import sys
import tempfile
import unittest
import numpy as np

import pyrate.core.orbital
from pyrate.core import ifgconstants as ifc, config as cf
from pyrate.core.ref_phs_est import ReferencePhaseError
from pyrate.core.shared import CorrectionStatusError
from pyrate import prepifg, process, conv2tif
from pyrate.configuration import Configuration
from tests import common

legacy_ref_phs_method1 = [-18.2191658020020,
                          27.7119445800781,
                          -18.4944229125977,
                          -2.92210483551025,
                          31.1168708801270,
                          21.2123012542725,
                          9.01810073852539,
                          6.08130645751953,
                          -3.79313516616821,
                          -11.3826837539673,
                          -7.28352737426758,
                          17.6365375518799,
                          -12.8673439025879,
                          5.46325922012329,
                          -35.4149475097656,
                          -13.5371961593628,
                          -12.7864856719971]


legacy_ref_phs_method2 = [-21.4459648132324,
                          27.1714553833008,
                          -20.8264484405518,
                          -3.47468209266663,
                          30.4519863128662,
                          22.3201427459717,
                          9.58487224578857,
                          4.81979084014893,
                          -3.89160847663879,
                          -12.0131330490112,
                          -8.64702987670898,
                          19.2060871124268,
                          -9.92049789428711,
                          4.38952684402466,
                          -34.9590339660645,
                          -14.3167810440063,
                          -11.9066228866577]


class RefPhsTests(unittest.TestCase):
    """Basic reference phase estimation tests"""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.params = dict()
        self.params[cf.REF_EST_METHOD] = 1
        self.params[cf.PARALLEL] = False
        self.params[cf.TMPDIR] = self.tmp_dir
        self.refpx, self.refpy = 38, 58
        common.copytree(common.SML_TEST_TIF, self.tmp_dir)
        small_tifs = glob.glob(os.path.join(self.tmp_dir, "*.tif"))
        for s in small_tifs:
            os.chmod(s, 0o644)
        self.ifgs = common.small_data_setup(self.tmp_dir, is_dir=True)
        for ifg in self.ifgs:
            ifg.close()
       

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except PermissionError:
            print("Files still in use.")

        for ifg in self.ifgs:
            ifg.close()

    def test_need_at_least_two_ifgs(self):
        self.assertRaises(ReferencePhaseError, process._ref_phase_estimation,
                          self.ifgs[:1], self.params, self.refpx, self.refpy)

    def test_metadata(self):
        process._ref_phase_estimation(self.ifgs, self.params, self.refpx, self.refpy)
        for ifg in self.ifgs:
            ifg.open()
            self.assertEqual(ifg.dataset.GetMetadataItem(ifc.PYRATE_REF_PHASE),
                             ifc.REF_PHASE_REMOVED)
    
    def test_mixed_metadata_raises(self):
        # correct reference phase for some of the ifgs
        process._ref_phase_estimation(self.ifgs[:5], self.params, self.refpx, self.refpy)
        for ifg in self.ifgs:
            ifg.open()

        # now it should raise exception if we wnat to correct
        # refernece phase again on all of them
        self.assertRaises(CorrectionStatusError, process._ref_phase_estimation,
                          self.ifgs, self.params, self.refpx, self.refpy)
        

class RefPhsEstimationLegacyTestMethod1Serial(unittest.TestCase):
    """
    Reference phase estimation method 1 is tested vs legacy output
    """

    @classmethod
    def setUpClass(cls):

        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', common.TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = cls.temp_out_dir
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = False

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST],
                                               params[cf.OBS_DIR])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for ifg in ifgs:
            ifg.close()

        cls.ref_phs, cls.ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.temp_out_dir)
        except PermissionError:
            print("File opened by another process.")

        try:
            common.remove_tifs(cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])
        except PermissionError:
            print("File opened by another process.")

        for ifg in cls.ifgs:
            ifg.close()

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method1,
                                             self.ref_phs,
                                             decimal=3)

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            if not ifg.is_open:
                ifg.open()

        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR,
                                                     'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_unw_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data), np.nanmean(self.ifgs[k].phase_data), places=3)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)), np.sum(np.isnan(self.ifgs[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))


class RefPhsEstimationLegacyTestMethod1Parallel(unittest.TestCase):
    """
    Reference phase estimation method 1 is tested vs legacy output
    """
    @classmethod
    def setUpClass(cls):

        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', common.TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = cls.temp_out_dir
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = True

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST],
                                               params[cf.OBS_DIR])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for i in ifgs:
            i.close()

        cls.ref_phs, cls.ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)

        # end run_pyrate copy


    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        shutil.rmtree(cls.temp_out_dir)
        common.remove_tifs(
            cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method1,
                                             self.ref_phs,
                                             decimal=3)

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR,
                                                     'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_unw_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(
                        ifg_data,
                        self.ifgs[k].phase_data,
                        decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data), np.nanmean(self.ifgs[k].phase_data), places=3)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))


class RefPhsEstimationLegacyTestMethod2Serial(unittest.TestCase):
    """
    Reference phase estimation method 2 is tested vs legacy output
    """

    @classmethod
    def setUpClass(cls):

        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', common.TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = cls.temp_out_dir
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = False

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST],
                                               params[cf.OBS_DIR])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop,
                                               params, xlks)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        
        for i in ifgs:
            i.close()

        cls.ref_phs, cls.ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)

    @classmethod
    def tearDownClass(cls):
        for ifg in cls.ifgs:
            ifg.close()
        shutil.rmtree(cls.temp_out_dir)
        common.remove_tifs(
            cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])


    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR,
                                                     'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_unw_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data), np.nanmean(self.ifgs[k].phase_data), places=3)

                    # number of nans must equal
                    self.assertEqual(np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))


        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method2,
                                             self.ref_phs, decimal=3)


class RefPhsEstimationLegacyTestMethod2Parallel(unittest.TestCase):
    """
    Reference phase estimation method 2 is tested vs legacy output

    """
    # TODO: Improve the parallel tests to remove duplication from serial tests
    @classmethod
    def setUpClass(cls):

        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        cls.temp_out_dir = tempfile.mkdtemp()
        sys.argv = ['prepifg.py', common.TEST_CONF_ROIPAC]
        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.TMPDIR] = cls.temp_out_dir
        conv2tif.main(params)
        prepifg.main(params)

        params[cf.OUT_DIR] = cls.temp_out_dir
        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = True

        xlks, ylks, crop = cf.transform_params(params)

        base_ifg_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST],    
                                               params[cf.OBS_DIR])

        dest_paths = cf.get_dest_paths(base_ifg_paths, crop,
                                       params, xlks)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = process._ref_pixel_calc(dest_paths, params)

        # Estimate and remove orbit errors
        pyrate.core.orbital.remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for i in ifgs:
            i.close()

        cls.ref_phs, cls.ifgs = process._ref_phase_estimation(dest_paths, params, refx, refy)

    @classmethod
    def tearDownClass(cls):
        for ifg in cls.ifgs:
            ifg.close()
        shutil.rmtree(cls.temp_out_dir)
        common.remove_tifs(
            cf.get_config_params(common.TEST_CONF_ROIPAC)[cf.OBS_DIR])

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR,
                                            'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                     if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                     and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_unw_1rlks')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(
                        ifg_data, self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    self.assertAlmostEqual(np.nanmean(ifg_data), np.nanmean(self.ifgs[k].phase_data), places=3)

                    # number of nans must equal
                    self.assertEqual(
                        np.sum(np.isnan(ifg_data)),
                        np.sum(np.isnan(self.ifgs[k].phase_data)))

        # ensure we have the correct number of matches
        self.assertEqual(count, len(self.ifgs))

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method2, self.ref_phs, decimal=3)
