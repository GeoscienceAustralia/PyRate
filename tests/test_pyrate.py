#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
import pytest
from os.path import join
from pathlib import Path
import numpy as np

import pyrate.configuration
import pyrate.constants as C
import pyrate.core.prepifg_helper
import pyrate.core.shared
import pyrate.main
import tests.common
from pyrate.core import shared, prepifg_helper
from pyrate.core.shared import dem_or_ifg
from pyrate import correct
from pyrate.configuration import MultiplePaths, Configuration
from tests import common
from tests.common import PY37GDAL304

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
    params = {C.IFG_LKSX: 3, C.IFG_LKSY: 2, C.IFG_CROP_OPT: 1}
    assert pyrate.core.prepifg_helper.transform_params(params) == (3, 2, 1)


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
    paths = tests.common.original_ifg_paths(ifglist_path, ifgdir)
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


class TestPyRate:
    # Initialise & run workflow from class setup, ignoring multilooking as it is
    # a separate step. Unit tests verify different steps have completed

    @classmethod
    def setup_class(cls):

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

            config = Configuration(common.TEST_CONF_ROIPAC)
            params = config.__dict__
            params['correct'] = ['orbfit', 'refphase', 'mst', 'apscorrect', 'maxvar']
            params[C.OUT_DIR] = cls.BASE_OUT_DIR
            params[C.PROCESSOR] = 0  # roipac
            params[C.APS_CORRECTION] = 0
            paths = glob.glob(join(cls.BASE_OUT_DIR, 'geo_*-*.tif'))
            paths = sorted(paths)
            params[C.PARALLEL] = False
            params[C.ORBFIT_OFFSET] = True
            params[C.TEMP_MLOOKED_DIR] = cls.BASE_OUT_DIR.join(C.TEMP_MLOOKED_DIR)
            params[C.INTERFEROGRAM_FILES] = [MultiplePaths(p, params) for p in paths]
            for p in params[C.INTERFEROGRAM_FILES]:  # cheat
                p.sampled_path = p.converted_path
                p.tmp_sampled_path = p.converted_path
            params["rows"], params["cols"] = 2, 2
            params[C.REF_PIXEL_FILE] = Configuration.ref_pixel_path(params)
            Path(params[C.OUT_DIR]).joinpath(C.APS_ERROR_DIR).mkdir(exist_ok=True, parents=True)
            Path(params[C.OUT_DIR]).joinpath(C.MST_DIR).mkdir(exist_ok=True, parents=True)
            correct.correct_ifgs(config)

            if not hasattr(cls, 'ifgs'):
                cls.ifgs = get_ifgs(out_dir=cls.BASE_OUT_DIR)
        except:
            # revert working dir & avoid paths busting other tests
            os.chdir(CURRENT_DIR)
            raise

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.BASE_DIR, ignore_errors=True)
        os.chdir(CURRENT_DIR)

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        assert key in md, 'Missing %s in %s' % (key, ifg.data_path)
        assert md[key] == value

    def test_basic_outputs(self):
        assert os.path.exists(self.BASE_OUT_DIR)

        for i in self.ifgs:
            assert ~i.is_read_only

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


@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL304, reason="Only run in one CI env")
class TestParallelPyRate:
    """
    parallel vs serial pyrate tests verifying results from all steps equal
    """

    @classmethod
    def setup_class(cls):
        gamma_conf = common.TEST_CONF_GAMMA
        from tests.common import manipulate_test_conf
        rate_types = ['stack_rate', 'stack_error', 'stack_samples']
        cls.tif_dir = Path(tempfile.mkdtemp())
        params = manipulate_test_conf(gamma_conf, cls.tif_dir)

        from pyrate.configuration import Configuration
        # change the required params
        params[C.PROCESSES] = 4
        params[C.PROCESSOR] = 1  # gamma
        params[C.IFG_FILE_LIST] = os.path.join(common.GAMMA_SML_TEST_DIR, 'ifms_17')
        params[C.PARALLEL] = 1
        params[C.APS_CORRECTION] = 0
        params[C.REFX], params[C.REFY] = -1, -1
        rows, cols = params["rows"], params["cols"]

        output_conf_file = 'gamma.conf'
        output_conf = cls.tif_dir.joinpath(output_conf_file).as_posix()
        pyrate.configuration.write_config_file(params=params, output_conf_file=output_conf)

        config = Configuration(output_conf)
        params = config.__dict__

        common.sub_process_run(f"pyrate conv2tif -f {output_conf}")
        common.sub_process_run(f"pyrate prepifg -f {output_conf}")

        cls.sampled_paths = [p.tmp_sampled_path for p in params[C.INTERFEROGRAM_FILES]]

        ifgs = common.small_data_setup()
        correct._copy_mlooked(params)
        tiles = pyrate.core.shared.get_tiles(cls.sampled_paths[0], rows, cols)
        correct.correct_ifgs(config)
        pyrate.main.timeseries(config)
        pyrate.main.stack(config)
        cls.refpixel_p, cls.maxvar_p, cls.vcmt_p = \
            (params[C.REFX], params[C.REFY]), params[C.MAXVAR], params[
                C.VCMT]
        cls.mst_p = common.reconstruct_mst(ifgs[0].shape, tiles, params[C.OUT_DIR])
        cls.rate_p, cls.error_p, cls.samples_p = \
            [common.reconstruct_stack_rate(ifgs[0].shape, tiles, params[C.TMPDIR], t) for t in rate_types]
        
        common.remove_tifs(params[C.OUT_DIR])

        # now create the non parallel version
        cls.tif_dir_s = Path(tempfile.mkdtemp())
        params = manipulate_test_conf(gamma_conf, cls.tif_dir_s)
        params[C.PROCESSES] = 4
        params[C.PROCESSOR] = 1  # gamma
        params[C.IFG_FILE_LIST] = os.path.join(common.GAMMA_SML_TEST_DIR, 'ifms_17')
        params[C.PARALLEL] = 0
        params[C.APS_CORRECTION] = 0
        params[C.REFX], params[C.REFY] = -1, -1
        output_conf_file = 'gamma.conf'
        output_conf = cls.tif_dir_s.joinpath(output_conf_file).as_posix()
        pyrate.configuration.write_config_file(params=params, output_conf_file=output_conf)
        config = Configuration(output_conf)
        params = config.__dict__

        common.sub_process_run(f"pyrate conv2tif -f {output_conf}")
        common.sub_process_run(f"pyrate prepifg -f {output_conf}")

        correct._copy_mlooked(params)
        correct.correct_ifgs(config)
        pyrate.main.timeseries(config)
        pyrate.main.stack(config)
        cls.refpixel, cls.maxvar, cls.vcmt = \
            (params[C.REFX], params[C.REFY]), params[C.MAXVAR], params[
                C.VCMT]
        cls.mst = common.reconstruct_mst(ifgs[0].shape, tiles, params[C.OUT_DIR])
        cls.rate, cls.error, cls.samples = \
            [common.reconstruct_stack_rate(ifgs[0].shape, tiles, params[C.TMPDIR], t) for t in rate_types]
        cls.params = params

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[C.OUT_DIR])

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'
        ifgs = [dem_or_ifg(i) for i in self.sampled_paths]
        for i in ifgs:
            i.open()
            i.nodata_value = 0
            self.key_check(i, key, value)

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        assert key in md, 'Missing %s in %s' % (key, ifg.data_path)
        assert md[key] == value

    def test_phase_conversion(self):
        # ensure phase has been converted from radians to millimetres
        key = 'DATA_UNITS'
        value = 'MILLIMETRES'
        ifgs = [dem_or_ifg(i) for i in self.sampled_paths]

        for i in ifgs:
            i.open()
            i.nodata_value = 0
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


@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL304, reason="Only run in one CI env")
class TestPrePrepareIfgs:

    @classmethod
    def setup_class(cls):
        params = Configuration(common.TEST_CONF_ROIPAC).__dict__
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

        nan_conversion = params[C.NAN_CONVERSION]

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
            i.open(readonly=False)
            if nan_conversion:  # nan conversion happens here in networkx mst
                i.nodata_value = params[C.NO_DATA_VALUE]
                i.convert_to_nans()
            if not i.mm_converted:
                i.convert_to_mm()
            i.close()

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.tmp_dir2)
        shutil.rmtree(cls.tmp_dir)

    def test_small_data_prep_phase_equality(self):
        for i, j in zip(self.ifgs, self.ifg_ret):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data)
            assert ~(i.phase_data == 0).any()
            # if there was any 0 still present
            i.phase_data[4, 2] = 0
            assert (i.phase_data == 0).any()

    def test_small_data_prep_metadata_equality(self):
        for i, j in zip(self.ifgs, self.ifg_ret):
            assert i.meta_data == j.meta_data
