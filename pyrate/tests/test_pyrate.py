"""
Tests end-to-end runs and functions in PyRate workflow module.

Created on 17/09/2012

.. codeauthor:: Ben Davies
"""

import os, glob, shutil, logging, unittest
from os.path import join
from pyrate import shared, config, prepifg
from pyrate.scripts import run_pyrate
from pyrate import config as cf

# testing constants
BASE_DIR = '/tmp/pyrate/workflow'
BASE_OUT_DIR = join(BASE_DIR, 'out')
BASE_DEM_DIR = join(BASE_DIR, 'dem')
BASE_CFG_FILE = join(BASE_DIR, 'pyrate_workflow_test.conf')
BASE_DEM_FILE = join(BASE_DEM_DIR, 'sydney_trimmed.tif')

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
    assert dest_paths == ['out/ifg0.tif', 'out/ifg1.tif']


# FIXME: change to read output ifgs
def get_ifgs(_open=True):
    paths = glob.glob(join(BASE_OUT_DIR, 'geo_*-*.tif'))
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
        from pyrate.tests.common import SYD_TEST_MATLAB_ORBITAL_DIR
        # start each full test run cleanly
        shutil.rmtree(BASE_DIR, ignore_errors=True)

        try:
            # copy source data (treat as prepifg already run)
            os.makedirs(BASE_OUT_DIR)

            for path in glob.glob(join(TEST_CORE, 'tif/*')):
                dest = join(BASE_OUT_DIR, os.path.basename(path))
                shutil.copy(path, dest)
                os.chmod(dest, 0660)

            os.makedirs(BASE_DEM_DIR)
            orig_dem = join(TEST_CORE, 'dem', 'sydney_trimmed.tif')
            os.symlink(orig_dem, BASE_DEM_FILE)
            os.chdir(BASE_DIR)

            # manually start logging as main() is being bypassed
            run_pyrate.init_logging(logging.DEBUG)

            params = config.get_config_params(
                os.path.join(SYD_TEST_MATLAB_ORBITAL_DIR, 'orbital_error.conf'))
            params[cf.SIM_DIR] = cf.PYRATEPATH
            paths = glob.glob(join(BASE_OUT_DIR, 'geo_*-*.tif'))
            run_pyrate.process_ifgs(paths, params)
            os.chdir(CURRENT_DIR)
        except:
            # revert working dir & avoid paths busting other tests
            os.chdir(CURRENT_DIR)
            raise

    def setUp(self):
        # performance: to save constantly opening ifgs
        if not hasattr(self, 'ifgs'):
            self.ifgs = get_ifgs()

    def get_logfile_path(self):
        logpaths = glob.glob(join(BASE_DIR, '*.log'))

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
        self.assertTrue(os.path.exists(BASE_DIR))
        self.assertTrue(os.path.exists(BASE_OUT_DIR))

        for i in self.ifgs:
            self.assertFalse(i.is_read_only)

        log_path = self.get_logfile_path()
        st = os.stat(log_path)
        self.assertTrue(st.st_size > 0)

    def test_wavelength_conversion(self):
        # ensure phase has been converted from metres to millimetres
        key = 'PHASE_UNITS'
        value = 'MILLIMETRES'

        for i in self.ifgs:
            self.key_check(i, key, value)

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in self.ifgs:
            self.key_check(i, key, value)

    #def test_mst_matrix(self):
        # TODO: should matrix be written to disk to be read in?
        #raise NotImplementedError

    def test_refpixel_found(self):
        log_path = self.get_logfile_path()
        for line in file(log_path):
            if 'Reference pixel coordinate:' in line:
                return

        self.fail('No reference pixel found')


if __name__ == "__main__":
    unittest.main()
