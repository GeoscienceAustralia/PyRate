"""
Tests end-to-end runs and functions in PyRate workflow module.

Created on 17/09/2012

.. codeauthor:: Ben Davies, Sudipta Basak
"""

import os, glob, shutil, logging, unittest
from os.path import join
import tempfile
import numpy as np
from pyrate import shared, config, prepifg
from pyrate.scripts import run_pyrate, run_prepifg
from pyrate import config as cf
from pyrate.tests import common
from pyrate.shared import Ifg
from pyrate import mst

# taken from http://stackoverflow.com/questions/6260149/os-symlink-support-in-windows
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


TEST_CORE = join(os.environ['PYRATEPATH'], 'tests', 'sydney_test')

CURRENT_DIR = os.getcwd()


def test_transform_params():
    params = {config.IFG_LKSX: 3, config.IFG_LKSY: 2, config.IFG_CROP_OPT: 1}
    assert run_pyrate.transform_params(params) == (3, 2, 1)


def test_warp_required():
    nocrop = prepifg.ALREADY_SAME_SIZE
    assert run_pyrate.warp_required(xlooks=2, ylooks=1, crop=nocrop)
    assert run_pyrate.warp_required(xlooks=1, ylooks=2, crop=nocrop)
    assert run_pyrate.warp_required(xlooks=1, ylooks=1, crop=nocrop)
    assert not run_pyrate.warp_required(xlooks=1, ylooks=1, crop=None)

    for c in prepifg.CROP_OPTIONS[:-1]:
        assert run_pyrate.warp_required(xlooks=1, ylooks=1, crop=c)


def test_original_ifg_paths():
    ifgdir = join(TEST_CORE, 'tif')
    ifglist_path = join(ifgdir, 'ifms_17')
    paths = run_pyrate.original_ifg_paths(ifglist_path)
    assert paths[0] == join(ifgdir, 'geo_060619-061002.tif'), str(paths[0])
    assert paths[-1] == join(ifgdir, 'geo_070709-070813.tif')


# def test_working_ifg_paths():
#     # working ifg paths is not used anymore
#     # base file list changes to the correct mlooked file even though
#     # CROPOUT is the same as the original tif
#     src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
#
#     assert run_pyrate.working_ifg_paths(src_paths, 1, 1, 4) == src_paths
#     assert run_pyrate.working_ifg_paths(src_paths, 0, 0, 4) == src_paths


def test_working_ifg_paths_mlooked_exists_failure():
    src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
    try:
        _ = run_pyrate.working_ifg_paths(src_paths, 2, 2, 4)
        raise Exception("Didn't fail with missing files")
    except IOError:
        return


def test_dest_ifg_paths():
    # given source ifgs to process, get paths of ifgs in out dir
    src_paths = ['tif/ifg0.tif', 'tif/ifg1.tif']
    dest_paths = run_pyrate.dest_ifg_paths(src_paths, outdir='out')
    assert dest_paths == [os.path.join('out', i) for i in ['ifg0.tif', 'ifg1.tif']]


# FIXME: change to read output ifgs
def get_ifgs(out_dir, _open=True):
    paths = glob.glob(join(out_dir, 'geo_*-*.tif'))
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
        cls.BASE_DEM_FILE = join(cls.BASE_DEM_DIR, 'sydney_trimmed.tif')
        from pyrate.tests.common import SYD_TEST_DIR

        try:
            # copy source data (treat as prepifg already run)
            os.makedirs(cls.BASE_OUT_DIR)
            for path in glob.glob(join(TEST_CORE, 'tif/*')):
                dest = join(cls.BASE_OUT_DIR, os.path.basename(path))
                shutil.copy(path, dest)
                os.chmod(dest, 0660)

            os.makedirs(cls.BASE_DEM_DIR)
            orig_dem = join(TEST_CORE, 'dem', 'sydney_trimmed.tif')
            os.symlink(orig_dem, cls.BASE_DEM_FILE)
            os.chdir(cls.BASE_DIR)

            # manually start logging as main() is being bypassed
            run_pyrate.init_logging(logging.DEBUG)

            params = config.get_config_params(
                os.path.join(SYD_TEST_DIR, 'pyrate_system_test.conf'))
            params[cf.SIM_DIR] = cf.PYRATEPATH
            params[cf.OUT_DIR] = cls.BASE_OUT_DIR
            params[cf.PROCESSOR] = 0  # roipac
            params[cf.APS_CORRECTION] = 0
            paths = glob.glob(join(cls.BASE_OUT_DIR, 'geo_*-*.tif'))
            params[cf.PARALLEL] = False
            run_pyrate.process_ifgs(paths, params)

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

    def get_logfile_path(self):
        logpaths = glob.glob(join(self.BASE_DIR, '*.log'))
        if len(logpaths) != 1:
            msg = 'Log not generated. Use --nologcapture if running nosetests'
            self.fail(msg)
        return logpaths[0]

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        self.assertTrue(key in md, 'Missing %s in %s' % (key, ifg.data_path))
        self.assertTrue(md[key], value)

    def test_basic_outputs(self):
        self.assertTrue(os.path.exists(self.BASE_OUT_DIR))

        for i in self.ifgs:
            self.assertFalse(i.is_read_only)

        log_path = self.get_logfile_path()
        st = os.stat(log_path)
        self.assertTrue(st.st_size > 0)

    def test_phase_conversion(self):
        # ensure phase has been converted to millimetres
        key = 'PHASE_UNITS'
        value = 'MILLIMETRES'

        for i in self.ifgs:
            self.key_check(i, key, value)

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in self.ifgs:
            self.key_check(i, key, value)

    def test_refpixel_found(self):
        log_path = self.get_logfile_path()
        for line in file(log_path):
            if 'Reference pixel coordinate:' in line:
                return

        self.fail('No reference pixel found')


class ParallelPyRateTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tif_dir = tempfile.mkdtemp()
        cls.test_conf = common.SYDNEY_TEST_CONF

        # change the required params
        params = cf.get_config_params(cls.test_conf)
        params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        params[cf.PROCESSOR] = 1  # gamma
        params[cf.IFG_FILE_LIST] = os.path.join(
            common.SYD_TEST_GAMMA, 'ifms_17')
        params[cf.OUT_DIR] = cls.tif_dir
        params[cf.PARALLEL] = 1
        params[cf.APS_CORRECTION] = False

        xlks, ylks, crop = run_pyrate.transform_params(params)

        # base_unw_paths need to be geotiffed and multilooked by run_prepifg
        base_unw_paths = run_pyrate.original_ifg_paths(params[cf.IFG_FILE_LIST])

        # dest_paths are tifs that have been geotif converted and multilooked
        cls.dest_paths = run_pyrate.get_dest_paths(
            base_unw_paths, crop, params, xlks)

        run_prepifg.gamma_prepifg(base_unw_paths, params)

        cls.mst_p, cls.refpixel_p, cls.maxvar_p, cls.vcmt_p, cls.rate_p, \
            cls.error_p, cls.samples_p = \
            run_pyrate.process_ifgs(cls.dest_paths, params)
        cls.mst_p_2 = run_pyrate.mst_calculation(cls.dest_paths, params)

        # now create the non parallel version
        cls.tif_dir_s = tempfile.mkdtemp()
        params[cf.PARALLEL] = 0
        params[cf.OUT_DIR] = cls.tif_dir_s
        cls.dest_paths_s = run_pyrate.get_dest_paths(
            base_unw_paths, crop, params, xlks)
        run_prepifg.gamma_prepifg(base_unw_paths, params)
        cls.mst, cls.refpixel, cls.maxvar, cls.vcmt, cls.rate, \
        cls.error, cls.samples = \
            run_pyrate.process_ifgs(cls.dest_paths_s, params)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tif_dir, ignore_errors=True)
        shutil.rmtree(cls.tif_dir_s, ignore_errors=True)

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in common.sydney_data_setup(datafiles=self.dest_paths):
            self.key_check(i, key, value)

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        self.assertTrue(key in md, 'Missing %s in %s' % (key, ifg.data_path))
        self.assertTrue(md[key], value)

    def test_phase_conversion(self):
        # ensure phase has been converted to millimetres
        key = 'PHASE_UNITS'
        value = 'MILLIMETRES'

        for i in common.sydney_data_setup(datafiles=self.dest_paths):
            self.key_check(i, key, value)

    def test_mst_equal(self):
        from pyrate import mst
        ifgs = common.sydney_data_setup(datafiles=self.dest_paths)
        mst_original_p = mst.mst_boolean_array(ifgs)
        ifgs_s = common.sydney_data_setup(datafiles=self.dest_paths_s)
        mst_original_s = mst.mst_boolean_array(ifgs_s)
        np.testing.assert_array_equal(self.mst, mst_original_p)
        np.testing.assert_array_equal(self.mst, mst_original_s)
        np.testing.assert_array_equal(self.mst, self.mst_p_2)
        np.testing.assert_array_equal(self.mst, self.mst_p)

    def test_refpixel_equal(self):
        np.testing.assert_array_equal(self.refpixel, self.refpixel_p)

    def test_maxvar_equal(self):
        np.testing.assert_array_equal(self.maxvar, self.maxvar_p)

    def test_vcmt_equal(self):
        np.testing.assert_array_equal(self.vcmt, self.vcmt_p)

    def test_linear_rate_equal(self):
        np.testing.assert_array_almost_equal(self.rate, self.rate_p)
        np.testing.assert_array_almost_equal(self.error, self.error_p)
        np.testing.assert_array_almost_equal(self.samples, self.samples_p)


class TestPrePrepareIfgs(unittest.TestCase):

    def test_sydney_data_prep(self):
        tmp_dir = tempfile.mkdtemp()
        shared.copytree(common.SYD_TEST_TIF, tmp_dir)
        tifs = glob.glob(os.path.join(tmp_dir, "*.tif"))
        for t in tifs:
            os.chmod(t, 0644)
        sydney_ifgs = common.sydney_data_setup(datafiles=tifs)
        params = config.get_config_params(common.SYDNEY_TEST_CONF)

        ifg_paths = [i.data_path for i in sydney_ifgs]

        ifg_ret = run_pyrate.pre_prepare_ifgs(ifg_paths, params=params)
        for i in ifg_ret:
            i.close()

        nan_conversion = params[cf.NAN_CONVERSION]

        ifgs = [Ifg(p) for p in ifg_paths]

        for i in ifgs:
            if not i.is_open:
                i.open(readonly=False)
            if nan_conversion:  # nan conversion happens here in networkx mst
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_nans()
            if not i.mm_converted:
                i.convert_to_mm()

        for i, j in zip(ifgs, ifg_ret):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data)
            self.assertFalse((i.phase_data == 0).any())
            # if there was any 0 still present
            i.phase_data[4, 2] = 0
            self.assertTrue((i.phase_data == 0).any())

if __name__ == "__main__":
    unittest.main()
