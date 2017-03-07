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
This Python module contains tests for the aps.py PyRate module.
"""
import copy
import glob
import os
import re
import shutil
import subprocess
import tempfile
import unittest
import pytest
import gdal
import numpy as np

from pyrate import config as cf
from pyrate import ifgconstants as ifc
from pyrate.compat import PyAPS_INSTALLED
from pyrate.scripts import run_pyrate, run_prepifg
from tests import common
if PyAPS_INSTALLED:
    from pyrate import aps


@unittest.skipUnless(PyAPS_INSTALLED, 'PyAPS must be available for this test')
@pytest.mark.skipif(not PyAPS_INSTALLED,
                    reason='PyAPS must be available for this test')
class TestMethod1VsMethod2AndMetaData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.tif_dir_method2 = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        file_list = cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                   'ifms_17'))
        cls.params[cf.IFG_FILE_LIST] = tempfile.mktemp(dir=cls.tif_dir)
        # write a short filelist with only 3 gamma unws
        with open(cls.params[cf.IFG_FILE_LIST], 'w') as fp:
            for f in file_list[:2]:
                fp.write(os.path.join(common.SYD_TEST_GAMMA, f) + '\n')
        cls.params[cf.OUT_DIR] = cls.tif_dir
        cls.params[cf.PARALLEL] = True
        cls.params[cf.REF_EST_METHOD] = 1
        cls.params[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
        cls.params[cf.APS_INCIDENCE_MAP] = common.SYD_TEST_INCIDENCE
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])
        # add dem
        base_unw_paths.append(common.SYD_TEST_DEM_GAMMA)
        # add incidence
        base_unw_paths.append(common.SYD_TEST_INCIDENCE)

        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        import copy
        cls.params_method2 = copy.copy(cls.params)
        cls.params_method2[cf.OUT_DIR] = cls.tif_dir_method2
        cls.params_method2[cf.APS_METHOD] = 2

        # dest_paths are tifs that have been geotif converted and multilooked
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params)
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params_method2)

        # removed incidence as we don't want it in ifgs list
        base_unw_paths.pop()
        # removed dem as we don't want it in ifgs list
        base_unw_paths.pop()

        dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params, xlks)
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

        dest_paths_m2 = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params_method2, xlks)
        cls.ifgs_method2 = common.sydney_data_setup(datafiles=dest_paths_m2)
        aps.remove_aps_delay(cls.ifgs, cls.params)
        aps.remove_aps_delay(cls.ifgs_method2, cls.params_method2)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        for i in cls.ifgs_method2:
            i.close()

        shutil.rmtree(cls.tif_dir)
        shutil.rmtree(cls.tif_dir_method2)

    def test_metadata_was_copied(self):
        for i in self.ifgs:
            md = i.meta_data
            i.close()
            self.assertIn(ifc.PYRATE_APS_ERROR, md.keys())
            self.assertIn(aps.APS_STATUS, md.values())

    def test_meta_data_was_written(self):
        for i in self.ifgs:
            md = i.meta_data
            i.close()
            ds = gdal.Open(i.data_path)
            md_w = ds.GetMetadata()
            self.assertDictEqual(md, md_w)
            ds = None

    def test_dem_tifs_present(self):
        # geotiffed dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                                    + '.tif'))

        # multilooked dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                        + '_{looks}rlks_{crop}cr.tif'.format(
                           looks=self.params[cf.IFG_LKSX],
                           crop=self.params[cf.IFG_CROP_OPT])))

    def test_method1_method2_equal_with_uniform_incidence_map(self):
        for i, j in zip(self.ifgs, self.ifgs_method2):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data,
                                                 decimal=4)


@unittest.skipUnless(PyAPS_INSTALLED, 'PyAPS must be available for this test')
@pytest.mark.skipif(not PyAPS_INSTALLED,
                    reason='PyAPS must be available for this test')
class TestOriginalVsEfficientAps(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.tif_dir_original = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        file_list = cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                   'ifms_17'))
        cls.params[cf.IFG_FILE_LIST] = tempfile.mktemp(dir=cls.tif_dir)
        # write a short filelist with only 3 gamma unws
        with open(cls.params[cf.IFG_FILE_LIST], 'w') as fp:
            for f in file_list[:2]:
                fp.write(os.path.join(common.SYD_TEST_GAMMA, f) + '\n')
        cls.params[cf.OUT_DIR] = cls.tif_dir
        cls.params[cf.PARALLEL] = True
        cls.params[cf.REF_EST_METHOD] = 2
        cls.params[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
        cls.params[cf.APS_INCIDENCE_MAP] = common.SYD_TEST_INCIDENCE
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])
        # add dem
        base_unw_paths.append(common.SYD_TEST_DEM_GAMMA)
        # add incidence
        base_unw_paths.append(common.SYD_TEST_INCIDENCE)

        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        import copy
        cls.params_original = copy.copy(cls.params)
        cls.params_original[cf.OUT_DIR] = cls.tif_dir_original
        cls.params_original[cf.APS_METHOD] = 2

        # dest_paths are tifs that have been geotif converted and multilooked
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params)
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params_original)

        # removed incidence as we don't want it in ifgs list
        base_unw_paths.pop()
        # removed dem as we don't want it in ifgs list
        base_unw_paths.pop()

        dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params, xlks)
        cls.ifgs = common.sydney_data_setup(datafiles=dest_paths)

        dest_paths_orig = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params_original, xlks)
        cls.ifgs_orig = common.sydney_data_setup(datafiles=dest_paths_orig)
        aps.remove_aps_delay(cls.ifgs, cls.params)
        aps.remove_aps_delay_original(cls.ifgs_orig, cls.params_original)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs:
            i.close()
        for i in cls.ifgs_orig:
            i.close()

        shutil.rmtree(cls.tif_dir)
        shutil.rmtree(cls.tif_dir_original)

    def test_metadata_was_copied(self):
        for i in self.ifgs:
            md = i.meta_data
            i.close()
            self.assertIn(ifc.PYRATE_APS_ERROR, md.keys())
            self.assertIn(aps.APS_STATUS, md.values())

    def test_meta_data_was_written(self):
        for i in self.ifgs:
            i.close()
            md = i.meta_data
            ds = gdal.Open(i.data_path)
            md_w = ds.GetMetadata()
            self.assertDictEqual(md, md_w)
            ds = None


    def test_dem_tifs_present(self):
        # geotiffed dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                                    + '.tif'))

        # multilooked dem
        os.path.exists(os.path.join(self.params[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                        + '_{looks}rlks_{crop}cr.tif'.format(
                           looks=self.params[cf.IFG_LKSX],
                           crop=self.params[cf.IFG_CROP_OPT])))

    def test_method1_method2_equal_with_uniform_incidence_map(self):
        for i, j in zip(self.ifgs, self.ifgs_orig):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data,
                                                 decimal=4)


@unittest.skipUnless(PyAPS_INSTALLED, 'PyAPS must be available for this test')
@pytest.mark.skipif(not PyAPS_INSTALLED,
                    reason='PyAPS must be available for this test')
class TestAPSIncidenceVsElevationVsParallel(unittest.TestCase):
    """
    This class tests APS method when incidence map is provided vs elevation map
    """

    @classmethod
    def setUpClass(cls):
        cls.tif_dir_inc = tempfile.mkdtemp()
        cls.tif_dir_ele = tempfile.mkdtemp()
        cls.tif_dir_ele_par = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params_inc = cf.get_config_params(cls.test_conf)
        cls.params_inc[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params_inc[cf.PROCESSOR] = 1  # gamma
        file_list = cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                   'ifms_17'))
        # config file
        cls.params_inc[cf.IFG_FILE_LIST] = tempfile.mktemp(dir=cls.tif_dir_inc)

        # write a short filelist with only 3 gamma unws
        with open(cls.params_inc[cf.IFG_FILE_LIST], 'w') as fp:
            for f in file_list[:2]:
                fp.write(os.path.join(common.SYD_TEST_GAMMA, f) + '\n')

        cls.params_inc[cf.OUT_DIR] = cls.tif_dir_inc
        cls.params_inc[cf.PARALLEL] = 0
        cls.params_inc[cf.REF_EST_METHOD] = 1
        cls.params_inc[cf.APS_METHOD] = 2
        cls.params_inc[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
        cls.params_inc[cf.APS_INCIDENCE_MAP] = common.SYD_TEST_INCIDENCE
        run_prepifg.main(cls.params_inc)

        # now create the config for the elevation_map case
        cls.params_ele = copy.copy(cls.params_inc)
        cls.params_ele[cf.OUT_DIR] = cls.tif_dir_ele
        cls.params_ele[cf.APS_METHOD] = 2
        cls.params_ele[cf.APS_INCIDENCE_MAP] = None
        cls.params_ele[cf.APS_INCIDENCE_EXT] = None
        cls.params_ele[cf.APS_ELEVATION_MAP] = common.SYD_TEST_ELEVATION
        cls.params_ele[cf.APS_ELEVATION_EXT] = 'lv_theta'
        run_prepifg.main(cls.params_ele)

        ptn = re.compile(r'\d{8}')
        dest_paths_inc = [f for f in
                          glob.glob(os.path.join(cls.tif_dir_inc, '*.tif'))
                          if ("cr" in f) and ("rlks" in f) and
                          (len(re.findall(ptn, os.path.basename(f))) == 2)]
        cls.ifgs_inc = common.sydney_data_setup(datafiles=dest_paths_inc)

        dest_paths_ele = [f for f in
                          glob.glob(os.path.join(cls.tif_dir_ele, '*.tif'))
                          if "cr" in f and "rlks" in f and
                          (len(re.findall(ptn, os.path.basename(f))) == 2)]

        cls.ifgs_ele = common.sydney_data_setup(datafiles=dest_paths_ele)

        # now create the config for the elevation map parallel case
        cls.params_ele_par = copy.copy(cls.params_ele)
        cls.params_ele_par[cf.OUT_DIR] = cls.tif_dir_ele_par
        cls.params_ele_par[cf.PARALLEL] = True
        run_prepifg.main(cls.params_ele_par)
        dest_paths_ele_par = \
            [f for f in glob.glob(os.path.join(cls.tif_dir_ele_par, '*.tif'))
             if "cr" in f and "rlks" in f and
             (len(re.findall(ptn, os.path.basename(f))) == 2)]

        cls.ifgs_ele_par = common.sydney_data_setup(datafiles=dest_paths_ele_par)

        aps.remove_aps_delay(cls.ifgs_inc, cls.params_inc)
        aps.remove_aps_delay(cls.ifgs_ele, cls.params_ele)
        aps.remove_aps_delay(cls.ifgs_ele_par, cls.params_ele_par)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs_ele:
            i.close()
        for i in cls.ifgs_ele_par:
            i.close()
        for i in cls.ifgs_inc:
            i.close()

        shutil.rmtree(cls.tif_dir_inc)
        shutil.rmtree(cls.tif_dir_ele)
        shutil.rmtree(cls.tif_dir_ele_par)

    def test_inc_vs_ele_equal_with_uniform_incidence_map(self):
        for i, j in zip(self.ifgs_inc, self.ifgs_ele):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data,
                                                 decimal=4)

    def test_inc_serial_vs_parallel(self):
        for i, j in zip(self.ifgs_ele_par, self.ifgs_ele):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data,
                                                 decimal=4)


@unittest.skipUnless(PyAPS_INSTALLED, 'PyAPS must be available for this test')
class MPITests(unittest.TestCase):
    # TODO: add tests for looks > 1
    @classmethod
    def setUpClass(cls):
        cls.tif_dir_serial = tempfile.mkdtemp()
        cls.tif_dir_mpi = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        cls.params = cf.get_config_params(cls.test_conf)
        cls.params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        cls.params[cf.PROCESSOR] = 1  # gamma
        file_list = cf.parse_namelist(os.path.join(common.SYD_TEST_GAMMA,
                                                   'ifms_17'))
        cls.params[cf.IFG_FILE_LIST] = tempfile.mktemp(dir=cls.tif_dir_serial)
        # write a short filelist with only 3 gamma unws
        with open(cls.params[cf.IFG_FILE_LIST], 'w') as fp:
            for f in file_list[:2]:
                fp.write(os.path.join(common.SYD_TEST_GAMMA, f) + '\n')
        cls.params[cf.OUT_DIR] = cls.tif_dir_serial
        cls.params[cf.PARALLEL] = 0
        cls.params[cf.REF_EST_METHOD] = 2
        cls.params[cf.IFG_LKSX] = 1
        cls.params[cf.IFG_LKSY] = 1
        cls.params[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
        cls.params[cf.APS_INCIDENCE_MAP] = common.SYD_TEST_INCIDENCE
        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(
            cls.params[cf.IFG_FILE_LIST])
        # add dem
        base_unw_paths.append(common.SYD_TEST_DEM_GAMMA)
        # add incidence
        base_unw_paths.append(common.SYD_TEST_INCIDENCE)

        xlks, ylks, crop = run_pyrate.transform_params(cls.params)

        cls.params_mpi = copy.copy(cls.params)
        cls.params_mpi[cf.OUT_DIR] = cls.tif_dir_mpi

        # dest_paths are tifs that have been geotif converted and multilooked
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params)
        run_prepifg.gamma_prepifg(base_unw_paths, cls.params_mpi)

        # removed incidence as we don't want it in ifgs list
        base_unw_paths.pop()
        # removed dem as we don't want it in ifgs list
        base_unw_paths.pop()

        dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params, xlks)
        dest_paths_mpi = run_pyrate.get_dest_paths(
            base_unw_paths, crop, cls.params_mpi, xlks)
        run_pyrate.process_ifgs(dest_paths, cls.params)
        cls.ifgs_serial = common.sydney_data_setup(datafiles=dest_paths)
        cls.conf_mpi = tempfile.mktemp('.conf', dir=cls.tif_dir_mpi)
        cf.write_config_file(cls.params_mpi, cls.conf_mpi)
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar.py ' + \
              cls.conf_mpi
        cmd = str.split()
        subprocess.check_call(cmd)
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar_2.py ' + \
              cls.conf_mpi
        cmd = str.split()
        subprocess.check_call(cmd)
        str = 'mpirun -np 4 python pyrate/nci/run_pyrate_pypar_3.py ' + \
              cls.conf_mpi
        cmd = str.split()
        subprocess.check_call(cmd)
        cls.ifgs_mpi = common.sydney_data_setup(datafiles=dest_paths_mpi)

    @classmethod
    def tearDownClass(cls):
        for i in cls.ifgs_serial:
            i.close()
        for i in cls.ifgs_mpi:
            i.close()

        shutil.rmtree(cls.tif_dir_serial)
        shutil.rmtree(cls.tif_dir_mpi)

    def test_metadata_was_copied(self):
        for j, i in zip(self.ifgs_serial, self.ifgs_mpi):
            md = i.meta_data
            md_s = j.meta_data
            self.assertIn(ifc.PYRATE_APS_ERROR, md_s.keys())
            self.assertIn(ifc.PYRATE_APS_ERROR, md.keys())
            self.assertIn(aps.APS_STATUS, md.values())

    def test_meta_data_was_written(self):
        for i in self.ifgs_mpi:
            md = i.meta_data
            ds = gdal.Open(i.data_path)
            md_w = ds.GetMetadata()
            self.assertDictEqual(md, md_w)
            self.assertIn(ifc.PYRATE_APS_ERROR, md_w.keys())
            ds = None

    def test_dem_tifs_present(self):
        # geotiffed dem
        os.path.exists(os.path.join(self.params_mpi[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                                    + '.tif'))

        # multilooked dem
        os.path.exists(os.path.join(self.params_mpi[cf.OUT_DIR],
                       os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                        + '_{looks}rlks_{crop}cr.tif'.format(
                           looks=self.params_mpi[cf.IFG_LKSX],
                           crop=self.params_mpi[cf.IFG_CROP_OPT])))

    def test_serial_vs_mpi_equal(self):
        for i, j in zip(self.ifgs_serial, self.ifgs_mpi):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data,
                                                 decimal=4)


if __name__ == '__main__':
    unittest.main()
