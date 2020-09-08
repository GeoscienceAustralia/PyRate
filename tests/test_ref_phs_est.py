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
from pathlib import Path
import shutil
import tempfile
import pytest
import numpy as np

from pyrate.core import ifgconstants as ifc, config as cf
from pyrate.core.ref_phs_est import ReferencePhaseError, ref_phase_est_wrapper
from pyrate.core.refpixel import ref_pixel_calc_wrapper
from pyrate.core.orbital import remove_orbital_error
from pyrate.core.shared import CorrectionStatusError, Ifg
from pyrate import prepifg, correct, conv2tif
from pyrate.configuration import MultiplePaths, Configuration
from tests import common
from tests.common import TEST_CONF_GAMMA

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


class TestRefPhsTests:
    """Basic reference phase estimation tests"""

    def setup_method(self):
        self.params = Configuration(common.TEST_CONF_GAMMA).__dict__
        self.tmp_dir = tempfile.mkdtemp()
        self.params[cf.OUT_DIR] = self.tmp_dir
        self.params[cf.REF_EST_METHOD] = 1
        self.params[cf.PARALLEL] = False
        self.params[cf.TMPDIR] = self.tmp_dir
        common.copytree(common.SML_TEST_TIF, self.tmp_dir)
        self.small_tifs = glob.glob(os.path.join(self.tmp_dir, "*.tif"))
        for s in self.small_tifs:
            os.chmod(s, 0o644)
        self.ifgs = common.small_data_setup(self.tmp_dir, is_dir=True)
        self.params[cf.INTERFEROGRAM_FILES] = [MultiplePaths(p, self.params) for p in self.small_tifs]
        for p in self.params[cf.INTERFEROGRAM_FILES]:
            p.sampled_path = p.converted_path
            p.tmp_sampled_path = p.sampled_path
        for ifg in self.ifgs:
            ifg.close()

        self.params[cf.REFX], self.params[cf.REFY] = -1, -1
        self.params[cf.REFNX], self.params[cf.REFNY] = 10, 10
        self.params[cf.REF_CHIP_SIZE], self.params[cf.REF_MIN_FRAC] = 21, 0.5
        self.params['rows'], self.params['cols'] = 3, 2
        self.params[cf.REF_PIXEL_FILE] = Configuration.ref_pixel_path(self.params)
        correct._update_params_with_tiles(self.params)
        correct.ref_pixel_calc_wrapper(self.params)
       
    def teardown_method(self):
        shutil.rmtree(self.params[cf.OUT_DIR])

    def test_need_at_least_two_ifgs(self):
        self.params[cf.INTERFEROGRAM_FILES] = [MultiplePaths(p, self.params) for p in self.small_tifs[:1]]
        for p in self.params[cf.INTERFEROGRAM_FILES]:
            p.sampled_path = p.converted_path
            p.tmp_sampled_path = p.sampled_path

        with pytest.raises(ReferencePhaseError):
            ref_phase_est_wrapper(self.params)

    def test_metadata(self):
        for ifg in self.ifgs:
            ifg.open()
            assert ifc.PYRATE_REF_PHASE not in ifg.dataset.GetMetadata()
            ifg.close()
        ref_phase_est_wrapper(self.params)
        for ifg in self.ifgs:
            ifg.open()
            assert ifg.dataset.GetMetadataItem(ifc.PYRATE_REF_PHASE) == ifc.REF_PHASE_REMOVED
            ifg.close()
    
    def test_mixed_metadata_raises(self):

        # change config to 5 ifgs
        self.params[cf.INTERFEROGRAM_FILES] = [MultiplePaths(p, self.params) for p in self.small_tifs[:5]]
        for p in self.params[cf.INTERFEROGRAM_FILES]:
            p.sampled_path = p.converted_path
            p.tmp_sampled_path = p.sampled_path

        # correct reference phase for some of the ifgs
        ref_phase_est_wrapper(self.params)
        for ifg in self.ifgs:
            ifg.open()

        # change config to all ifgs
        self.params[cf.INTERFEROGRAM_FILES] = [MultiplePaths(p, self.params) for p in self.small_tifs]
        for p in self.params[cf.INTERFEROGRAM_FILES]:
            p.sampled_path = p.converted_path
            p.tmp_sampled_path = p.sampled_path

        # now it should raise exception if we want to correct refernece phase again on all of them
        with pytest.raises(CorrectionStatusError):
            ref_phase_est_wrapper(self.params)
        

class TestRefPhsEstimationLegacyTestMethod1Serial:
    """
    Reference phase estimation method 1 is tested vs legacy output
    """

    @classmethod
    def setup_class(cls):
        # start with a clean output dir
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        conv2tif.main(params)
        prepifg.main(params)
        for p in params[cf.INTERFEROGRAM_FILES]:  # hack
            p.tmp_sampled_path = p.sampled_path
            Path(p.sampled_path).chmod(0o664)  # assign write permission as conv2tif output is readonly
        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = False
        params[cf.ORBFIT_OFFSET] = True

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)
        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = ref_pixel_calc_wrapper(params)

        # Estimate and remove orbit errors
        remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for ifg in ifgs:
            ifg.close()

        for p in params[cf.INTERFEROGRAM_FILES]:
            p.tmp_sampled_path = p.sampled_path
        params[cf.REFX], params[cf.REFY] = refx, refy
        params['rows'], params['cols'] = 3, 2
        correct._update_params_with_tiles(params)
        cls.ref_phs, cls.ifgs = ref_phase_est_wrapper(params)
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method1, self.ref_phs, decimal=3)

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            if not ifg.is_open:
                ifg.open()

        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR, 'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_correctedgeo_')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_ifg.tif')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    assert np.nanmean(ifg_data) == pytest.approx(np.nanmean(self.ifgs[k].phase_data), abs=0.001)

                    # number of nans must equal
                    assert np.sum(np.isnan(ifg_data)) == np.sum(np.isnan(self.ifgs[k].phase_data))

        # ensure we have the correct number of matches
        assert count == len(self.ifgs)


class TestRefPhsEstimationLegacyTestMethod1Parallel:
    """
    Reference phase estimation method 1 is tested vs legacy output
    """
    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        conv2tif.main(params)
        prepifg.main(params)
        for p in params[cf.INTERFEROGRAM_FILES]:  # hack
            p.tmp_sampled_path = p.sampled_path
            Path(p.sampled_path).chmod(0o664)  # assign write permission as conv2tif output is readonly


        params[cf.REF_EST_METHOD] = 1
        params[cf.PARALLEL] = True
        params[cf.ORBFIT_OFFSET] = True

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = ref_pixel_calc_wrapper(params)

        # Estimate and remove orbit errors
        remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for i in ifgs:
            i.close()
        for p in params[cf.INTERFEROGRAM_FILES]:
            p.tmp_sampled_path = p.sampled_path
        params[cf.REFX], params[cf.REFY] = refx, refy
        params['rows'], params['cols'] = 3, 2
        correct._update_params_with_tiles(params)
        cls.ref_phs, cls.ifgs = ref_phase_est_wrapper(params)
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_estimate_reference_phase(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method1, self.ref_phs, decimal=3)

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR, 'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(
                LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_correctedgeo_')[-1].split('.')[0] == os.path.split(j.data_path)[-1].split(
                        '_ifg.tif')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(
                        ifg_data,
                        self.ifgs[k].phase_data,
                        decimal=3)

                    # means must also be equal
                    assert np.nanmean(ifg_data) == pytest.approx(np.nanmean(self.ifgs[k].phase_data), abs=0.001)

                    # number of nans must equal
                    assert np.sum(np.isnan(ifg_data)) == np.sum(np.isnan(self.ifgs[k].phase_data))

        # ensure we have the correct number of matches
        assert count == len(self.ifgs)


class TestRefPhsEstimationLegacyTestMethod2Serial:
    """
    Reference phase estimation method 2 is tested vs legacy output
    """

    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        conv2tif.main(params)
        prepifg.main(params)
        for p in params[cf.INTERFEROGRAM_FILES]:  # hack
            p.tmp_sampled_path = p.sampled_path
            Path(p.sampled_path).chmod(0o664)  # assign write permission as conv2tif output is readonly

        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = False
        params[cf.ORBFIT_OFFSET] = True

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        mst_grid = common.mst_calculation(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = ref_pixel_calc_wrapper(params)

        # Estimate and remove orbit errors
        remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        
        for i in ifgs:
            i.close()
        for p in params[cf.INTERFEROGRAM_FILES]:
            p.tmp_sampled_path = p.sampled_path
        params[cf.REFX], params[cf.REFY] = refx, refy
        params['rows'], params['cols'] = 3, 2
        correct._update_params_with_tiles(params)

        cls.ref_phs, cls.ifgs = ref_phase_est_wrapper(params)
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR, 'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2geo_')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_ifg.tif')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(ifg_data,
                        self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    assert np.nanmean(ifg_data) == pytest.approx(np.nanmean(self.ifgs[k].phase_data), abs=0.001)

                    # number of nans must equal
                    assert np.sum(np.isnan(ifg_data)) == np.sum(np.isnan(self.ifgs[k].phase_data))


        # ensure we have the correct number of matches
        assert count == len(self.ifgs)

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method2, self.ref_phs, decimal=3)


class TestRefPhsEstimationLegacyTestMethod2Parallel:
    """
    Reference phase estimation method 2 is tested vs legacy output

    """
    # TODO: Improve the parallel tests to remove duplication from serial tests

    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
        conv2tif.main(params)
        prepifg.main(params)
        for p in params[cf.INTERFEROGRAM_FILES]:  # hack
            p.tmp_sampled_path = p.sampled_path
            Path(p.sampled_path).chmod(0o664)  # assign write permission as conv2tif output is readonly

        params[cf.REF_EST_METHOD] = 2
        params[cf.PARALLEL] = 1
        params[cf.ORBFIT_OFFSET] = True

        dest_paths, headers = common.repair_params_for_correct_tests(params[cf.OUT_DIR], params)

        # start run_pyrate copy
        ifgs = common.pre_prepare_ifgs(dest_paths, params)
        # Estimate reference pixel location
        refx, refy = ref_pixel_calc_wrapper(params)

        # Estimate and remove orbit errors
        remove_orbital_error(ifgs, params)

        for i in ifgs:
            i.close()

        ifgs = common.pre_prepare_ifgs(dest_paths, params)

        for i in ifgs:
            i.close()

        for p in params[cf.INTERFEROGRAM_FILES]:
            p.tmp_sampled_path = p.sampled_path
        params[cf.REFX], params[cf.REFY] = refx, refy
        params['rows'], params['cols'] = 3, 2
        correct._update_params_with_tiles(params)
        cls.ref_phs, cls.ifgs = ref_phase_est_wrapper(params)
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_ifgs_after_ref_phs_est(self):
        for ifg in self.ifgs:
            ifg.open()
        LEGACY_REF_PHASE_DIR = os.path.join(common.SML_TEST_DIR, 'ref_phase_est')

        onlyfiles = [f for f in os.listdir(LEGACY_REF_PHASE_DIR)
                     if os.path.isfile(os.path.join(LEGACY_REF_PHASE_DIR, f))
                     and f.endswith('.csv') and f.__contains__('_ref_phase_')
                     and f.__contains__('method2')]

        count = 0
        for i, f in enumerate(onlyfiles):
            ifg_data = np.genfromtxt(os.path.join(LEGACY_REF_PHASE_DIR, f), delimiter=',')
            for k, j in enumerate(self.ifgs):
                if f.split('_corrected_method2geo_')[-1].split('.')[0] == \
                        os.path.split(j.data_path)[-1].split('_ifg.tif')[0]:
                    count += 1
                    # all numbers equal
                    np.testing.assert_array_almost_equal(
                        ifg_data, self.ifgs[k].phase_data, decimal=3)

                    # means must also be equal
                    assert np.nanmean(ifg_data) == pytest.approx(np.nanmean(self.ifgs[k].phase_data), abs=0.001)

                    # number of nans must equal
                    assert np.sum(np.isnan(ifg_data)) == np.sum(np.isnan(self.ifgs[k].phase_data))

        # ensure we have the correct number of matches
        assert count == len(self.ifgs)

    def test_estimate_reference_phase_method2(self):
        np.testing.assert_array_almost_equal(legacy_ref_phs_method2, self.ref_phs, decimal=3)


class TestRefPhsEstReusedFromDisc:

    @classmethod
    def setup_class(cls):
        cls.conf = TEST_CONF_GAMMA
        params = Configuration(cls.conf).__dict__
        conv2tif.main(params)
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[cf.OUT_DIR])

    def test_ref_phase_used_from_disc_on_rerun(self, ref_est_method):
        self.params = Configuration(self.conf).__dict__
        self.params[cf.REF_EST_METHOD] = ref_est_method
        correct._update_params_with_tiles(self.params)

        phase_prev, time_written = self.__run_once()

        # run again
        phase_now, time_written_1 = self.__run_once()

        # and once more
        phase_again, time_written_2 = self.__run_once()

        # assert no new file was written
        assert time_written_1 == time_written
        assert time_written_2 == time_written

        # assert phase data is unchanged after applying ref_ph correction from disc
        np.testing.assert_array_equal(phase_now, phase_prev)
        np.testing.assert_array_equal(phase_now, phase_again)

    def __run_once(self):
        ref_phs_file = Configuration.ref_phs_file(self.params)
        correct._copy_mlooked(self.params)
        multi_paths = self.params[cf.INTERFEROGRAM_FILES]
        ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        ifgs = [Ifg(i) for i in ifg_paths]
        self.params[cf.REFX_FOUND], self.params[cf.REFY_FOUND] = ref_pixel_calc_wrapper(self.params)
        correct._create_ifg_dict(self.params)
        ref_phase_est_wrapper(self.params)
        for i in ifgs:
            i.open()
        phase_prev = [i.phase_data for i in ifgs]
        # assert ref_ph_file present
        assert ref_phs_file.exists()
        time_written = os.stat(ref_phs_file).st_mtime
        for i in ifgs:
            i.close()
        return phase_prev, time_written
