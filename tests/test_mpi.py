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
This Python module contains tests for mpi operations in PyRate.
Tun this module as 'mpirun -n 4 pytest tests/test_mpi.py'
"""
from __future__ import print_function

import glob
import shutil
import numpy as np
import pytest
import os
import tempfile
import random
import string
from subprocess import check_output

import pyrate.orbital
import pyrate.shared
import tests.common
from pyrate import ref_phs_est as rpe
from pyrate import shared
from pyrate import vcm
from pyrate import refpixel
from pyrate.scripts import run_pyrate, run_prepifg, postprocessing
from tests.common import sydney_data_setup, reconstruct_mst, \
    reconstruct_linrate
from tests import common
from tests.test_vcm import matlab_maxvar
from pyrate import config as cf
from pyrate import mpiops
from pyrate import algorithm

TRAVIS = True if 'TRAVIS' in os.environ else False
PYTHON3P5 = True if ('TRAVIS_PYTHON_VERSION' in os.environ and
                     os.environ['TRAVIS_PYTHON_VERSION'] == '3.5') else False
GDAL_VERSION = check_output(["gdal-config", "--version"]).decode(
    encoding="utf-8").split('\n')[0]
MPITEST = TRAVIS and GDAL_VERSION == '2.0.0'


@pytest.fixture()
def tempdir():
    """
    tempdir for tests
    """
    def tmpdir():
        return tempfile.mkdtemp()
    return tmpdir


@pytest.fixture
def random_filename(tmpdir_factory):
    def make_random_filename(ext=''):
        dir = str(tmpdir_factory.mktemp('pyrate').realpath())
        fname = ''.join(random.choice(string.ascii_lowercase)
                        for _ in range(10))
        return os.path.join(dir, fname + ext)
    return make_random_filename


@pytest.fixture()
def get_config():
    """
    Parameters
    ----------
    conf_file: str
        config file

    Returns
    -------
    params: dict
        dict of params
    """
    def params(conf_file):
        return cf.get_config_params(conf_file)
    return params


# Make sure all MPI tests use this fixure
@pytest.fixture()
def mpisync(request):
    mpiops.comm.barrier()

    def fin():
        mpiops.comm.barrier()

    request.addfinalizer(fin)
    return mpiops.comm


@pytest.fixture(params=[0, 1])
def roipac_or_gamma(request):
    return request.param


@pytest.fixture(params=[1, 2])
def ref_est_method(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def row_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def col_splits(request):
    return request.param


@pytest.fixture(params=[1, 2, 5])
def modify_config(request, tempdir, get_config):
    test_conf = common.SYDNEY_TEST_CONF
    params_dict = get_config(test_conf)
    params_dict[cf.IFG_LKSX] = request.param
    params_dict[cf.IFG_LKSY] = request.param
    params_dict[cf.OBS_DIR] = tempdir()
    shared.copytree(common.SYD_TEST_GAMMA, params_dict[cf.OBS_DIR])
    params_dict[cf.IFG_FILE_LIST] = os.path.join(
        params_dict[cf.OBS_DIR], 'ifms_17')
    params_dict[cf.PARALLEL] = False
    params_dict[cf.APS_CORRECTION] = 0
    yield params_dict
    # clean up
    shutil.rmtree(params_dict[cf.OBS_DIR])


@pytest.fixture(params=range(1, 6))
def get_lks(request):
    return request.param


@pytest.fixture(params=range(1, 3))
def get_crop(request):
    return request.param


def test_vcm_matlab_vs_mpi(mpisync, tempdir, get_config):
    from tests.common import SYD_TEST_DIR

    params_dict = get_config(os.path.join(SYD_TEST_DIR,
                                          'pyrate_system_test.conf'))

    MATLAB_VCM_DIR = os.path.join(SYD_TEST_DIR, 'matlab_vcm')
    matlab_vcm = np.genfromtxt(os.path.join(MATLAB_VCM_DIR,
                                            'matlab_vcmt.csv'), delimiter=',')
    if mpiops.rank == 0:
        outdir = tempdir()
    else:
        outdir = None
    outdir = mpiops.comm.bcast(outdir, root=0)
    params_dict[cf.OUT_DIR] = outdir
    params_dict[cf.PARALLEL] = False
    xlks, ylks, crop = cf.transform_params(params_dict)
    base_unw_paths = cf.original_ifg_paths(params_dict[cf.IFG_FILE_LIST])
    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params_dict, xlks)

    # run prepifg, create the dest_paths files
    if mpiops.rank == 0:
        run_prepifg.roipac_prepifg(base_unw_paths, params_dict)

    mpiops.comm.barrier()

    tiles = run_pyrate.get_tiles(dest_paths[0], rows=1, cols=1)
    preread_ifgs = run_pyrate.create_ifg_dict(dest_paths,
                                              params=params_dict,
                                              tiles=tiles)
    refpx, refpy = run_pyrate.ref_pixel_calc(dest_paths, params_dict)
    run_pyrate.orb_fit_calc(dest_paths, params_dict)
    run_pyrate.ref_phase_estimation(dest_paths, params_dict, refpx, refpy)

    maxvar, vcmt = run_pyrate.maxvar_vcm_calc(dest_paths, params_dict,
                                              preread_ifgs)
    np.testing.assert_array_almost_equal(maxvar, matlab_maxvar, decimal=4)
    np.testing.assert_array_almost_equal(matlab_vcm, vcmt, decimal=3)
    if mpiops.rank == 0:
        shutil.rmtree(outdir)


@pytest.fixture(params=[1, 2, 5])
def orbfit_lks(request):
    return request.param


@pytest.fixture(params=[1, 2])
def orbfit_method(request):
    return request.param


@pytest.mark.skipif(not MPITEST, reason='skipping mpi tests in travis except '
                                        'python 3.5 and GDAL=2.0.0')
def test_timeseries_linrate_mpi(mpisync, tempdir, modify_config,
                                ref_est_method, row_splits, col_splits,
                                get_crop, orbfit_lks, orbfit_method):
    params = modify_config
    outdir = mpiops.run_once(tempdir)
    params[cf.OUT_DIR] = outdir
    params[cf.REF_EST_METHOD] = ref_est_method
    params[cf.IFG_CROP_OPT] = get_crop
    params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks
    params[cf.ORBITAL_FIT_LOOKS_X] = orbfit_lks
    params[cf.ORBITAL_FIT_METHOD] = orbfit_method
    xlks, ylks, crop = cf.transform_params(params)
    if xlks * col_splits > 45 or ylks * row_splits > 70:
        print('skipping test because lks and col_splits are not compatible')
        return

    # skip some tests in travis to run CI faster
    if TRAVIS and (xlks % 2 or row_splits % 2 or col_splits % 2
                   or orbfit_lks % 2):
        print('Skipping in travis env for faster CI run')
        return
    print("xlks={}, ref_est_method={}, row_splits={}, col_splits={}, "
          "get_crop={}, orbfit_lks={}, orbfit_method={}, "
          "rank={}".format(xlks, ref_est_method, row_splits, col_splits,
                           get_crop, orbfit_lks, orbfit_method, mpiops.rank))

    base_unw_paths = cf.original_ifg_paths(params[cf.IFG_FILE_LIST])
    # dest_paths are tifs that have been geotif converted and multilooked
    dest_paths = cf.get_dest_paths(base_unw_paths, crop, params, xlks)

    # run prepifg, create the dest_paths files
    if mpiops.rank == 0:
        run_prepifg.gamma_prepifg(base_unw_paths, params)

    mpiops.comm.barrier()

    (refpx, refpy), maxvar, vcmt = run_pyrate.process_ifgs(
        ifg_paths=dest_paths, params=params, rows=row_splits, cols=col_splits)

    tiles = mpiops.run_once(run_pyrate.get_tiles, dest_paths[0],
                            rows=row_splits, cols=col_splits)
    postprocessing.postprocess_linrate(row_splits, col_splits, params)
    postprocessing.postprocess_timeseries(row_splits, col_splits, params)
    ifgs_mpi_out_dir = params[cf.OUT_DIR]
    ifgs_mpi = sydney_data_setup(datafiles=dest_paths)

    # single process timeseries/linrate calculation
    if mpiops.rank == 0:
        params_old = modify_config
        params_old[cf.OUT_DIR] = tempdir()
        params_old[cf.REF_EST_METHOD] = ref_est_method
        params_old[cf.IFG_CROP_OPT] = get_crop
        params_old[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks
        params_old[cf.ORBITAL_FIT_LOOKS_X] = orbfit_lks
        params_old[cf.ORBITAL_FIT_METHOD] = orbfit_method
        xlks, ylks, crop = cf.transform_params(params_old)
        base_unw_paths = cf.original_ifg_paths(
            params_old[cf.IFG_FILE_LIST])
        dest_paths = cf.get_dest_paths(
            base_unw_paths, crop, params_old, xlks)
        run_prepifg.gamma_prepifg(base_unw_paths, params_old)

        ifgs = shared.pre_prepare_ifgs(dest_paths, params_old)
        mst_grid = tests.common.mst_calculation(dest_paths, params_old)
        refy, refx = refpixel.ref_pixel(ifgs, params_old)
        assert (refx == refpx) and (refy == refpy)  # both must match
        pyrate.orbital.remove_orbital_error(ifgs, params_old)
        ifgs = common.prepare_ifgs_without_phase(dest_paths, params_old)
        rpe.estimate_ref_phase(ifgs, params_old, refx, refy)
        ifgs = shared.pre_prepare_ifgs(dest_paths, params_old)
        maxvar_s = [vcm.cvd(i, params_old)[0] for i in ifgs]
        vcmt_s = vcm.get_vcmt(ifgs, maxvar)
        tsincr, tscum, _ = tests.common.compute_time_series(
            ifgs, mst_grid, params, vcmt)
        rate, error, samples = tests.common.calculate_linear_rate(
            ifgs, params_old, vcmt, mst_grid)
        mst_mpi = reconstruct_mst(ifgs[0].shape, tiles, ifgs_mpi_out_dir)
        np.testing.assert_array_almost_equal(mst_grid, mst_mpi)
        tsincr_mpi, tscum_mpi = reconstruct_times_series(ifgs[0].shape,
                                                         tiles,
                                                         ifgs_mpi_out_dir)

        rate_mpi, error_mpi, samples_mpi = \
            [reconstruct_linrate(ifgs[0].shape, tiles, ifgs_mpi_out_dir, t)
             for t in ['linrate', 'linerror', 'linsamples']]
        np.testing.assert_array_almost_equal(maxvar, maxvar_s)
        np.testing.assert_array_almost_equal(vcmt, vcmt_s)
        for i, j in zip(ifgs, ifgs_mpi):
            np.testing.assert_array_almost_equal(i.phase_data, j.phase_data)
        np.testing.assert_array_almost_equal(tsincr, tsincr_mpi, decimal=4)
        np.testing.assert_array_almost_equal(tscum, tscum_mpi, decimal=4)
        np.testing.assert_array_almost_equal(rate, rate_mpi, decimal=4)
        np.testing.assert_array_almost_equal(error, error_mpi, decimal=4)
        np.testing.assert_array_almost_equal(samples, samples_mpi, decimal=4)

        # assert linear rate output tifs are same
        _tifs_same(ifgs_mpi_out_dir, params_old[cf.OUT_DIR], 'linrate.tif')
        _tifs_same(ifgs_mpi_out_dir, params_old[cf.OUT_DIR], 'linerror.tif')
        _tifs_same(ifgs_mpi_out_dir, params_old[cf.OUT_DIR], 'linsamples.tif')

        # assert time series output tifs are same
        epochlist = algorithm.get_epochs(ifgs)[0]

        for i in range(tsincr.shape[2]):
            _tifs_same(ifgs_mpi_out_dir, params_old[cf.OUT_DIR],
                       'tsincr' + '_' + str(epochlist.dates[i + 1]) + ".tif")

        # 12 timeseries outputs
        assert i + 1 == tsincr.shape[2]
        shutil.rmtree(ifgs_mpi_out_dir)  # remove mpi out dir
        shutil.rmtree(params_old[cf.OUT_DIR])  # remove serial out dir


def _tifs_same(dir1, dir2, tif):
    linrate_tif_s = os.path.join(dir1, tif)
    linrate_tif_m = os.path.join(dir2, tif)
    common.assert_ifg_phase_equal(linrate_tif_m, linrate_tif_s)


@pytest.mark.skipif(TRAVIS, reason='skipping mpi tests in travis')
def reconstruct_times_series(shape, tiles, output_dir):
    tsincr_file_0 = os.path.join(output_dir, 'tsincr_{}.npy'.format(0))
    shape3 = np.load(tsincr_file_0).shape[2]

    tsincr_mpi = np.empty(shape=(shape + (shape3,)), dtype=np.float32)
    tscum_mpi = np.empty_like(tsincr_mpi, dtype=np.float32)

    for i, t in enumerate(tiles):
        tsincr_file_n = os.path.join(output_dir,
                                     'tsincr_{}.npy'.format(i))
        tsincr_mpi[t.top_left_y:t.bottom_right_y,
                   t.top_left_x: t.bottom_right_x, :] = np.load(tsincr_file_n)

        tscum_file_n = os.path.join(output_dir, 'tscuml_{}.npy'.format(i))

        tscum_mpi[t.top_left_y:t.bottom_right_y,
                  t.top_left_x: t.bottom_right_x, :] = np.load(tscum_file_n)

    return tsincr_mpi, tscum_mpi


def test_prepifg_mpi(mpisync, get_config, tempdir,
                     roipac_or_gamma, get_lks, get_crop):
    from tests.common import SYDNEY_TEST_CONF
    from os.path import join, basename
    params = get_config(SYDNEY_TEST_CONF)
    outdir = mpiops.run_once(tempdir)
    params[cf.OUT_DIR] = outdir
    params[cf.PROCESSOR] = roipac_or_gamma
    params[cf.PARALLEL] = False
    params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
    params[cf.IFG_CROP_OPT] = get_crop
    if roipac_or_gamma == 1:
        params[cf.IFG_FILE_LIST] = join(common.SYD_TEST_GAMMA, 'ifms_17')
        params[cf.OBS_DIR] = common.SYD_TEST_GAMMA
        params[cf.DEM_FILE] = common.SYD_TEST_DEM_GAMMA
    run_prepifg.main(params)

    if mpiops.rank == 0:
        params_s = get_config(SYDNEY_TEST_CONF)
        params_s[cf.OUT_DIR] = tempdir()
        params_s[cf.PARALLEL] = True
        params_s[cf.IFG_LKSX], params_s[cf.IFG_LKSY] = get_lks, get_lks
        params_s[cf.IFG_CROP_OPT] = get_crop
        if roipac_or_gamma == 1:
            base_unw_paths = glob.glob(join(common.SYD_TEST_GAMMA,
                                            "*_utm.unw"))
            run_prepifg.gamma_prepifg(base_unw_paths, params_s)
        else:
            base_unw_paths = glob.glob(join(common.SYD_TEST_OBS, "*.unw"))
            run_prepifg.roipac_prepifg(base_unw_paths, params_s)

        mpi_tifs = glob.glob(join(outdir, "*.tif"))
        serial_tifs = glob.glob(join(params[cf.OUT_DIR], "*.tif"))
        mpi_tifs.sort()
        serial_tifs.sort()
        # 17 geotifs, and 17 mlooked tifs
        assert len(mpi_tifs) == len(serial_tifs)
        for m_f, s_f in zip(mpi_tifs, serial_tifs):
            assert basename(m_f) == basename(s_f)

        shutil.rmtree(outdir)
        shutil.rmtree(params_s[cf.OUT_DIR])
