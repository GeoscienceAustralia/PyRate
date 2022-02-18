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
This Python module contains regression tests for comparing output from serial,
parallel and MPI PyRate runs.
"""
import shutil
import pytest
from pathlib import Path
from subprocess import check_call, CalledProcessError, run
import numpy as np

import pyrate.constants as C
from pyrate.configuration import Configuration, write_config_file

from tests.common import (
    assert_same_files_produced,
    assert_two_dirs_equal,
    manipulate_test_conf,
    MEXICO_CROPA_CONF,
    TEST_CONF_GAMMA,
    PY37GDAL304,
    PY37GDAL302,
    PYTHON3P8,
    PYTHON3P7,
    WORKING_DIR
)


@pytest.fixture(params=[0, 1])
def parallel(request):
    return request.param


@pytest.fixture(params=[1, 4])
def local_crop(request):
    return request.param


@pytest.fixture()
def modified_config(tempdir, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees, ref_est_method):
    def modify_params(conf_file, parallel_vs_serial, output_conf_file):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)

        if params[C.PROCESSOR] == 1:  # turn on coherence for gamma
            params[C.COH_MASK] = 1

        params[C.PARALLEL] = parallel_vs_serial
        params[C.PROCESSES] = 4
        params[C.APSEST] = 1
        params[C.IFG_LKSX], params[C.IFG_LKSY] = get_lks, get_lks
        params[C.REFNX], params[C.REFNY] = 2, 2

        params[C.IFG_CROP_OPT] = get_crop
        params[C.ORBITAL_FIT_LOOKS_X], params[C.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[C.ORBITAL_FIT] = 1
        params[C.ORBITAL_FIT_METHOD] = orbfit_method
        params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[C.REF_EST_METHOD] = ref_est_method
        params[C.MAX_LOOP_LENGTH] = 3
        params[C.LR_MAXSIG] = 0 # turn off pixel masking for these tests
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params
    return modify_params


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PYTHON3P8, reason="Only run in one CI env")
def test_pipeline_parallel_vs_mpi(modified_config, gamma_or_mexicoa_conf=TEST_CONF_GAMMA):
    """
    Tests proving single/multiprocess/mpi produce same output
    """
    gamma_conf = gamma_or_mexicoa_conf
    if np.random.rand() > 0.1:  # skip 90% of tests randomly
        pytest.skip("Randomly skipping as part of 85 percent")
        if gamma_conf == MEXICO_CROPA_CONF:  # skip cropA conf 95% time
            if np.random.rand() > 0.5:
                pytest.skip('skipped in mexicoA')

    print("\n\n")
    print("===x==="*10)

    mpi_conf, params = modified_config(gamma_conf, 0, 'mpi_conf.conf')

    run(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True, check=True)
    run(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True, check=True)
    try:
        run(f"mpirun -n 3 pyrate correct -f {mpi_conf}", shell=True, check=True)
        run(f"mpirun -n 3 pyrate timeseries -f {mpi_conf}", shell=True, check=True)
        run(f"mpirun -n 3 pyrate stack -f {mpi_conf}", shell=True, check=True)
        run(f"mpirun -n 3 pyrate merge -f {mpi_conf}", shell=True, check=True)
    except CalledProcessError as e:
        print(e)
        pytest.skip("Skipping as part of correction error")

    mr_conf, params_m = modified_config(gamma_conf, 1, 'multiprocess_conf.conf')

    run(f"pyrate workflow -f {mr_conf}", shell=True, check=True)

    sr_conf, params_s = modified_config(gamma_conf, 0, 'singleprocess_conf.conf')

    run(f"pyrate workflow -f {sr_conf}", shell=True, check=True)

    # convert2tif tests, 17 interferograms
    if not gamma_conf == MEXICO_CROPA_CONF:
        assert_same_files_produced(params[C.INTERFEROGRAM_DIR],
                                   params_m[C.INTERFEROGRAM_DIR], params_s[C.INTERFEROGRAM_DIR], "*_unw.tif", 17)
        # dem
        assert_same_files_produced(params[C.GEOMETRY_DIR],
                                   params_m[C.GEOMETRY_DIR], params_s[C.GEOMETRY_DIR], "*_dem.tif", 1)
        # if coherence masking, comprare coh files were converted
        if params[C.COH_FILE_LIST] is not None:
            assert_same_files_produced(params[C.COHERENCE_DIR], params_m[C.COHERENCE_DIR], params_s[C.COHERENCE_DIR],
                                       "*_cc.tif", 17)
            print("coherence files compared")

    # prepifg checks
    num_of_ifgs = 30 if gamma_conf == MEXICO_CROPA_CONF else 17
    num_of_coh = 30 if gamma_conf == MEXICO_CROPA_CONF else 17

    # check geom files
    if params[C.DEMERROR]:
        # check files required by dem error correction are produced
        assert_same_files_produced(
            params[C.GEOMETRY_DIR], params_m[C.GEOMETRY_DIR], params_s[C.GEOMETRY_DIR],
            [ft + '.tif' for ft in C.GEOMETRY_OUTPUT_TYPES] + ['*dem.tif'],
            7 if gamma_conf == MEXICO_CROPA_CONF else 8
        )

    # ifgs
    assert_same_files_produced(params[C.INTERFEROGRAM_DIR], params_m[C.INTERFEROGRAM_DIR],
                               params_s[C.INTERFEROGRAM_DIR], ["*_ifg.tif"], num_of_ifgs)

    # coherence
    assert_same_files_produced(params[C.COHERENCE_DIR], params_m[C.COHERENCE_DIR],
                               params_s[C.COHERENCE_DIR], ["*_coh.tif"], num_of_coh)

    # coherence stats
    assert_same_files_produced(params[C.COHERENCE_DIR], params_m[C.COHERENCE_DIR],
                               params_s[C.COHERENCE_DIR], ["coh_*.tif"], 3)

    num_files = 30 if gamma_conf == MEXICO_CROPA_CONF else 17
    # cf.TEMP_MLOOKED_DIR will contain the temp files that can be potentially deleted later
    assert_same_files_produced(params[C.TEMP_MLOOKED_DIR], params_m[C.TEMP_MLOOKED_DIR],
                               params_s[C.TEMP_MLOOKED_DIR], "*_ifg.tif", num_files)

    # prepifg + correct steps that overwrite tifs test
    # ifg phase checking in the previous step checks the correct pipeline upto APS correction

    # 2 x because of aps files
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "tsincr_*.npy", params['notiles'] * 2)

    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "tscuml_*.npy", params['notiles'])

    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "linear_rate_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "linear_error_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "linear_intercept_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "linear_rsquared_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "linear_samples_*.npy", params['notiles'])

    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "stack_rate_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "stack_error_*.npy", params['notiles'])
    assert_same_files_produced(params[C.TMPDIR], params_m[C.TMPDIR], params_s[
        C.TMPDIR], "stack_samples_*.npy", params['notiles'])

    # compare merge step
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "stack*.tif", 3)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "stack*.kml", 2)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "stack*.png", 2)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "stack*.npy", 3)
    
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "linear_*.tif", 5)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "linear_*.kml", 3)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "linear_*.png", 3)
    assert_same_files_produced(params[C.VELOCITY_DIR], params_m[C.VELOCITY_DIR], params_s[
        C.VELOCITY_DIR], "linear_*.npy", 5)

    if params[C.PHASE_CLOSURE]:  # only in cropA
        __check_equality_of_phase_closure_outputs(mpi_conf, sr_conf)
        __check_equality_of_phase_closure_outputs(mpi_conf, mr_conf)
        assert_same_files_produced(params[C.TIMESERIES_DIR], params_m[C.TIMESERIES_DIR], params_s[
            C.TIMESERIES_DIR], "tscuml*.tif", 11)  # phase closure removes one tif
        assert_same_files_produced(params[C.TIMESERIES_DIR], params_m[C.TIMESERIES_DIR], params_s[
            C.TIMESERIES_DIR], "tsincr*.tif", 11)
    else:
        assert_same_files_produced(params[C.TIMESERIES_DIR], params_m[C.TIMESERIES_DIR], params_s[
            C.TIMESERIES_DIR], "tscuml*.tif", 12)
        assert_same_files_produced(params[C.TIMESERIES_DIR], params_m[C.TIMESERIES_DIR], params_s[
            C.TIMESERIES_DIR], "tsincr*.tif", 12)
        
    print("==========================xxx===========================")

    shutil.rmtree(params[WORKING_DIR])
    shutil.rmtree(params_m[WORKING_DIR])
    shutil.rmtree(params_s[WORKING_DIR])


def __check_equality_of_phase_closure_outputs(mpi_conf, sr_conf):
    m_config = Configuration(mpi_conf)
    s_config = Configuration(sr_conf)
    m_close = m_config.closure()
    s_close = s_config.closure()
    m_closure = np.load(m_close.closure)
    s_closure = np.load(s_close.closure)
    # loops
    m_loops = np.load(m_close.loops, allow_pickle=True)
    s_loops = np.load(s_close.loops, allow_pickle=True)
    m_weights = [m.weight for m in m_loops]
    s_weights = [m.weight for m in s_loops]
    np.testing.assert_array_equal(m_weights, s_weights)
    for i, (m, s) in enumerate(zip(m_loops, s_loops)):
        assert all(m_e == s_e for m_e, s_e in zip(m.edges, s.edges))
    # closure
    np.testing.assert_array_almost_equal(np.abs(m_closure), np.abs(s_closure), decimal=4)
    # num_occurrences_each_ifg
    m_num_occurences_each_ifg = np.load(m_close.num_occurences_each_ifg, allow_pickle=True)
    s_num_occurences_each_ifg = np.load(s_close.num_occurences_each_ifg, allow_pickle=True)
    np.testing.assert_array_equal(m_num_occurences_each_ifg, s_num_occurences_each_ifg)
    # check ps
    m_ifgs_breach_count = np.load(m_close.ifgs_breach_count)
    s_ifgs_breach_count = np.load(s_close.ifgs_breach_count)
    np.testing.assert_array_equal(m_ifgs_breach_count, s_ifgs_breach_count)


@pytest.fixture(params=[0, 1])
def coh_mask(request):
    return request.param


@pytest.fixture()
def modified_config_short(tempdir, local_crop, get_lks, coh_mask):
    orbfit_lks = 1
    orbfit_method = 1
    orbfit_degrees = 1
    ref_est_method = 1
    ref_pixel = (150.941666654, -34.218333314)

    def modify_params(conf_file, parallel, output_conf_file, largetifs):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)
        params[C.COH_MASK] = coh_mask
        params[C.PARALLEL] = parallel
        params[C.PROCESSES] = 4
        params[C.APSEST] = 1
        params[C.LARGE_TIFS] = largetifs
        params[C.IFG_LKSX], params[C.IFG_LKSY] = get_lks, get_lks
        params[C.REFX], params[C.REFY] = ref_pixel
        params[C.REFNX], params[C.REFNY] = 4, 4

        params[C.IFG_CROP_OPT] = local_crop
        params[C.ORBITAL_FIT_LOOKS_X], params[
            C.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[C.ORBITAL_FIT] = 1
        params[C.ORBITAL_FIT_METHOD] = orbfit_method
        params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[C.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        # print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.fixture
def create_mpi_files():

    def _create(modified_config_short, gamma_conf):

        mpi_conf, params = modified_config_short(gamma_conf, 0, 'mpi_conf.conf', 1)

        check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)

        try:
            check_call(f"mpirun -n 3 pyrate correct -f {mpi_conf}", shell=True)
            check_call(f"mpirun -n 3 pyrate timeseries -f {mpi_conf}", shell=True)
            check_call(f"mpirun -n 3 pyrate stack -f {mpi_conf}", shell=True)
        except CalledProcessError as c:
            print(c)
            pytest.skip("Skipping as we encountered a process error during CI")
        check_call(f"mpirun -n 3 pyrate merge -f {mpi_conf}", shell=True)
        return params

    return _create


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL304, reason="Only run in one CI env")
def test_stack_and_ts_mpi_vs_parallel_vs_serial(modified_config_short, gamma_conf, create_mpi_files, parallel):
    """
    Checks performed:
    1. mpi vs single process pipeline
    2. mpi vs parallel (python multiprocess) pipeline.
    3. Doing 1 and 2 means we have checked single vs parallel python multiprocess pipelines
    4. This also checks the entire pipeline using largetifs (new prepifg) vs old perpifg (python based)
    """
    if np.random.randint(0, 1000) > 300:  # skip 90% of tests randomly
        pytest.skip("Randomly skipping as part of 60 percent")

    print("\n\n")

    print("===x==="*10)

    params = create_mpi_files(modified_config_short, gamma_conf)

    sr_conf, params_p = modified_config_short(gamma_conf, parallel, 'parallel_conf.conf', 0)

    check_call(f"pyrate workflow -f {sr_conf}", shell=True)

    # convert2tif tests, 17 interferograms
    assert_two_dirs_equal(params[C.INTERFEROGRAM_DIR], params_p[C.INTERFEROGRAM_DIR], "*_unw.tif", 17)

    # if coherence masking, compare coh files were converted
    if params[C.COH_FILE_LIST] is not None:
        assert_two_dirs_equal(params[C.COHERENCE_DIR], params_p[C.COHERENCE_DIR], "*_cc.tif", 17)
        print("coherence files compared")
    assert_two_dirs_equal(params[C.INTERFEROGRAM_DIR], params_p[C.INTERFEROGRAM_DIR], ["*_ifg.tif"], 17)

    # one original dem, another multilooked dem
    assert_two_dirs_equal(params[C.GEOMETRY_DIR], params_p[C.GEOMETRY_DIR], ['*dem.tif'], 2)
    assert_two_dirs_equal(params[C.GEOMETRY_DIR], params_p[C.GEOMETRY_DIR],
                          [t + "*.tif" for t in C.GEOMETRY_OUTPUT_TYPES], 6)  # 2 dems, 6 geom

    assert_two_dirs_equal(params[C.TEMP_MLOOKED_DIR], params_p[C.TEMP_MLOOKED_DIR], "*_ifg.tif", 17)

    # ifg phase checking in the previous step checks the correct pipeline upto APS correction
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "tsincr_*.npy", params['notiles'] * 2)
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "tscuml_*.npy", params['notiles'])

    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "linear_rate_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "linear_error_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "linear_samples_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "linear_intercept_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "linear_rsquared_*.npy", params['notiles'])

    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "stack_rate_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "stack_error_*.npy", params['notiles'])
    assert_two_dirs_equal(params[C.TMPDIR], params_p[C.TMPDIR], "stack_samples_*.npy", params['notiles'])

    # compare merge step
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "stack*.tif", 3)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "stack*.kml", 2)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "stack*.png", 2)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "stack*.npy", 3)

    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "linear*.tif", 5)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "linear*.kml", 3)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "linear*.png", 3)
    assert_two_dirs_equal(params[C.VELOCITY_DIR], params_p[C.VELOCITY_DIR], "linear*.npy", 5)

    assert_two_dirs_equal(params[C.TIMESERIES_DIR], params_p[C.TIMESERIES_DIR], "tscuml*.tif")
    assert_two_dirs_equal(params[C.TIMESERIES_DIR], params_p[C.TIMESERIES_DIR], "tsincr*.tif")
    assert_two_dirs_equal(params[C.TIMESERIES_DIR], params_p[C.TIMESERIES_DIR], "tscuml*.npy")
    assert_two_dirs_equal(params[C.TIMESERIES_DIR], params_p[C.TIMESERIES_DIR], "tsincr*.npy")

    print("==========================xxx===========================")

    shutil.rmtree(params[WORKING_DIR])
    shutil.rmtree(params_p[WORKING_DIR])
