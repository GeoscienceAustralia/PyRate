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
This Python module contains system integration tests for the PyRate workflow.
"""

import glob
import os
import shutil
import tempfile
import unittest
from os.path import join
import numpy as np

import pyrate.core.shared
from pyrate.core import shared, config as cf, config, prepifg_helper, mst
from pyrate import process, prepifg, conv2tif
from pyrate.configuration import MultiplePaths
from tests import common

# taken from
# http://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
if os.name == "nt":
    def symlink_ms(source, link_name):
        import ctypes
        csl = ctypes.windll.kernel32.CreateSymbolicLinkW
        csl.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint32)
        csl.restype = ctypes.c_ubyte
        flags = 1 if os.path.isdir(source) else 0
        try:
            if csl(link_name, source.replace('/', '\\'), flags) == 0:
                raise ctypes.WinError()
        except:
            pass
    os.symlink = symlink_ms

CURRENT_DIR = os.getcwd()


def test_transform_params():
    params = {config.IFG_LKSX: 3, config.IFG_LKSY: 2, config.IFG_CROP_OPT: 1}
    assert cf.transform_params(params) == (3, 2, 1)


def test_warp_required():
    nocrop = prepifg_helper.ALREADY_SAME_SIZE
    assert shared.warp_required(xlooks=2, ylooks=1, crop=nocrop)
    assert shared.warp_required(xlooks=1, ylooks=2, crop=nocrop)
    assert shared.warp_required(xlooks=1, ylooks=1, crop=nocrop)
    assert not shared.warp_required(xlooks=1, ylooks=1, crop=None)

    for c in prepifg_helper.CROP_OPTIONS[:-1]:
        assert shared.warp_required(xlooks=1, ylooks=1, crop=c)


def test_original_ifg_paths():
    ifgdir = common.SML_TEST_TIF
    ifglist_path = join(ifgdir, 'ifms_17')
    paths = cf.original_ifg_paths(ifglist_path, ifgdir)
    assert paths[0] == join(ifgdir, 'geo_060619-061002_unw.tif'), str(paths[0])
    assert paths[-1] == join(ifgdir, 'geo_070709-070813_unw.tif')


def dest_ifg_paths(ifg_paths, outdir):
    """
    Returns paths to out/dest ifgs.
    """

    bases = [os.path.basename(p) for p in ifg_paths]
    return [join(outdir, p) for p in bases]


def test_dest_ifg_paths():
    # given source ifgs to process, get paths of ifgs in out dir
    src_paths = ['tif/ifg0.tif', 'tif/ifg1.tif']
    dest_paths = dest_ifg_paths(src_paths, outdir='out')
    assert dest_paths == [os.path.join('out', i) for i in ['ifg0.tif',
                                                           'ifg1.tif']]


# FIXME: change to read output ifgs
def get_ifgs(out_dir, _open=True):
    paths = glob.glob(join(out_dir, 'geo_*-*_unw.tif'))
    ifgs = [shared.Ifg(p) for p in paths]
    assert len(ifgs) == 17, 'Got %s' % ifgs

    if _open:
        for i in ifgs:
            i.open(readonly=False)
    return ifgs


class PyRateTests(unittest.TestCase):
    # Initialise & run workflow from class setup, ignoring multilooking as it is
    # a separate step. Unit tests verify different steps have completed

    @classmethod
    def setUpClass(cls):

        # testing constants2
        cls.BASE_DIR = tempfile.mkdtemp()
        cls.BASE_OUT_DIR = join(cls.BASE_DIR, 'out')
        cls.BASE_DEM_DIR = join(cls.BASE_DIR, 'dem')
        cls.BASE_DEM_FILE = join(cls.BASE_DEM_DIR, 'roipac_test_trimmed.tif')

        try:
            # copy source data (treat as prepifg already run)
            os.makedirs(cls.BASE_OUT_DIR)
            for path in glob.glob(join(common.SML_TEST_TIF, '*')):
                dest = join(cls.BASE_OUT_DIR, os.path.basename(path))
                shutil.copy(path, dest)
                os.chmod(dest, 0o660)

            os.makedirs(cls.BASE_DEM_DIR)
            orig_dem = common.SML_TEST_DEM_TIF
            os.symlink(orig_dem, cls.BASE_DEM_FILE)
            os.chdir(cls.BASE_DIR)

            # Turn off validation because we're in a different working dir
            #  and relative paths in config won't be work.
            params = config.get_config_params(common.TEST_CONF_ROIPAC)
            params[cf.OUT_DIR] = cls.BASE_OUT_DIR
            params[cf.PROCESSOR] = 0  # roipac
            params[cf.APS_CORRECTION] = 0
            paths = glob.glob(join(cls.BASE_OUT_DIR, 'geo_*-*.tif'))
            params[cf.PARALLEL] = False
            params[cf.INTERFEROGRAM_FILES] = [MultiplePaths(cls.BASE_OUT_DIR, p) for p in paths]
            for p in params[cf.INTERFEROGRAM_FILES]:  # cheat
                p.sampled_path = p.converted_path
            process.process_ifgs(sorted(paths), params, 2, 2)

            if not hasattr(cls, 'ifgs'):
                cls.ifgs = get_ifgs(out_dir=cls.BASE_OUT_DIR)
        except:
            # revert working dir & avoid paths busting other tests
            os.chdir(CURRENT_DIR)
            raise

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.BASE_DIR, ignore_errors=True)
        os.chdir(CURRENT_DIR)

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        self.assertTrue(key in md, 'Missing %s in %s' % (key, ifg.data_path))
        self.assertTrue(md[key], value)

    def test_basic_outputs(self):
        self.assertTrue(os.path.exists(self.BASE_OUT_DIR))

        for i in self.ifgs:
            self.assertFalse(i.is_read_only)

        # log_path = self.get_logfile_path()
        # st = os.stat(log_path)
        # self.assertTrue(st.st_size > 0)

    def test_phase_conversion(self):
        # ensure phase has been converted from radians to millimetres
        key = 'DATA_UNITS'
        value = 'MILLIMETRES'

        for i in self.ifgs:
            self.key_check(i, key, value)

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in self.ifgs:
            self.key_check(i, key, value)


class ParallelPyRateTests(unittest.TestCase):
    """
    parallel vs serial pyrate tests verifying results from all steps equal
    """

    @classmethod
    def setUpClass(cls):
        rate_types = ['stack_rate', 'stack_error', 'stack_samples']
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.TEST_CONF_GAMMA

        from pyrate.configuration import Configuration
        # change the required params
        params = Configuration(cls.test_conf).__dict__
        params[cf.OBS_DIR] = common.SML_TEST_GAMMA
        params[cf.PROCESSES] = 4
        params[cf.PROCESSOR] = 1  # gamma
        params[cf.IFG_FILE_LIST] = os.path.join(common.SML_TEST_GAMMA, 'ifms_17')
        params[cf.OUT_DIR] = cls.tif_dir
        params[cf.PARALLEL] = 1
        params[cf.APS_CORRECTION] = False
        params[cf.TMPDIR] = os.path.join(params[cf.OUT_DIR], cf.TMPDIR)
        rows, cols = params["rows"], params["cols"]

        # xlks, ylks, crop = cf.transform_params(params)

        # base_unw_paths need to be geotiffed by converttogeotif 
        #  and multilooked by run_prepifg
        base_unw_paths = list(cf.parse_namelist(params[cf.IFG_FILE_LIST]))

        multi_paths = [MultiplePaths(params[cf.OUT_DIR], b, ifglksx=params[cf.IFG_LKSX],
                                     ifgcropopt=params[cf.IFG_CROP_OPT]) for b in base_unw_paths]

        # dest_paths are tifs that have been geotif converted and multilooked
        cls.converted_paths = [b.converted_path for b in multi_paths]
        cls.sampled_paths = [b.sampled_path for b in multi_paths]
        from copy import copy
        orig_params = copy(params)
        conv2tif.main(params)
        prepifg.main(orig_params)
        tiles = pyrate.core.shared.get_tiles(cls.sampled_paths[0], rows, cols)
        ifgs = common.small_data_setup()
        params[cf.INTERFEROGRAM_FILES] = multi_paths

        cls.refpixel_p, cls.maxvar_p, cls.vcmt_p = process.process_ifgs(cls.sampled_paths, params, rows, cols)
        cls.mst_p = common.reconstruct_mst(ifgs[0].shape, tiles, params[cf.TMPDIR])
        cls.rate_p, cls.error_p, cls.samples_p = \
            [common.reconstruct_stack_rate(ifgs[0].shape, tiles, params[cf.TMPDIR], t) for t in rate_types]
        
        common.remove_tifs(params[cf.OBS_DIR])

        # now create the non parallel version
        cls.tif_dir_s = tempfile.mkdtemp()
        params[cf.PARALLEL] = 0
        params[cf.PROCESSES] = 1
        params[cf.OUT_DIR] = cls.tif_dir_s
        params[cf.TMPDIR] = os.path.join(params[cf.OUT_DIR], cf.TMPDIR)
        multi_paths = [MultiplePaths(params[cf.OUT_DIR], b, ifglksx=params[cf.IFG_LKSX],
                                     ifgcropopt=params[cf.IFG_CROP_OPT]) for b in base_unw_paths]

        cls.converted_paths_s = [b.converted_path for b in multi_paths]
        cls.sampled_paths_s = [b.sampled_path for b in multi_paths]
        orig_params = copy(params)
        conv2tif.main(params)
        prepifg.main(orig_params)
        params[cf.INTERFEROGRAM_FILES] = multi_paths
        cls.refpixel, cls.maxvar, cls.vcmt = process.process_ifgs(cls.sampled_paths_s, params, rows, cols)
        cls.mst = common.reconstruct_mst(ifgs[0].shape, tiles, params[cf.TMPDIR])
        cls.rate, cls.error, cls.samples = \
            [common.reconstruct_stack_rate(ifgs[0].shape, tiles, params[cf.TMPDIR], t) for t in rate_types]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir, ignore_errors=True)
        shutil.rmtree(cls.tif_dir_s, ignore_errors=True)
        common.remove_tifs(cf.get_config_params(cls.test_conf)[cf.OBS_DIR])

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in common.small_data_setup(datafiles=self.sampled_paths):
            self.key_check(i, key, value)

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        self.assertTrue(key in md, 'Missing %s in %s' % (key, ifg.data_path))
        self.assertTrue(md[key], value)

    def test_phase_conversion(self):
        # ensure phase has been converted from radians to millimetres
        key = 'DATA_UNITS'
        value = 'MILLIMETRES'

        for i in common.small_data_setup(datafiles=self.sampled_paths):
            self.key_check(i, key, value)

    def test_mst_equal(self):
        np.testing.assert_array_equal(self.mst, self.mst_p)

    def test_refpixel_equal(self):
        np.testing.assert_array_equal(self.refpixel, self.refpixel_p)

    def test_maxvar_equal(self):
        np.testing.assert_array_almost_equal(self.maxvar, self.maxvar_p, decimal=4)

    def test_vcmt_equal(self):
        np.testing.assert_array_almost_equal(self.vcmt, self.vcmt_p, decimal=4)

    def test_stack_rate_equal(self):
        np.testing.assert_array_almost_equal(self.rate, self.rate_p, decimal=4)
        np.testing.assert_array_almost_equal(self.error, self.error_p, decimal=4)
        np.testing.assert_array_almost_equal(self.samples, self.samples_p, decimal=4)


class TestPrePrepareIfgs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        params = config.get_config_params(common.TEST_CONF_ROIPAC)
        cls.tmp_dir = tempfile.mkdtemp()
        common.copytree(common.SML_TEST_TIF, cls.tmp_dir)
        tifs = glob.glob(os.path.join(cls.tmp_dir, "*.tif"))
        for t in tifs:
            os.chmod(t, 0o644)
        small_ifgs = common.small_data_setup(datafiles=tifs)
        ifg_paths = [i.data_path for i in small_ifgs]

        cls.ifg_ret = common.pre_prepare_ifgs(ifg_paths, params=params)
        for i in cls.ifg_ret:
            i.close()

        nan_conversion = params[cf.NAN_CONVERSION]

        # prepare a second set
        cls.tmp_dir2 = tempfile.mkdtemp()
        common.copytree(common.SML_TEST_TIF, cls.tmp_dir2)
        tifs = glob.glob(os.path.join(cls.tmp_dir2, "*.tif"))
        for t in tifs:
            os.chmod(t, 0o644)
        small_ifgs = common.small_data_setup(datafiles=tifs)
        ifg_paths = [i.data_path for i in small_ifgs]

        cls.ifgs = [shared.Ifg(p) for p in ifg_paths]

        for i in cls.ifgs:
            if not i.is_open:
                i.open(readonly=False)
            if nan_conversion:  # nan conversion happens here in networkx mst
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_nans()
            if not i.mm_converted:
                i.convert_to_mm()
            i.close()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dir2)
        shutil.rmtree(cls.tmp_dir)

    def test_small_data_prep_phase_equality(self):
        for i, j in zip(self.ifgs, self.ifg_ret):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data)
            self.assertFalse((i.phase_data == 0).any())
            # if there was any 0 still present
            i.phase_data[4, 2] = 0
            self.assertTrue((i.phase_data == 0).any())

    def test_small_data_prep_metadata_equality(self):
        for i, j in zip(self.ifgs, self.ifg_ret):
            self.assertDictEqual(i.meta_data, j.meta_data)
