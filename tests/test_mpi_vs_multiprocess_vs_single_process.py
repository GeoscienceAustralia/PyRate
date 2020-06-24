import os
import shutil
import pytest
from pathlib import Path
from subprocess import check_call, check_output, CalledProcessError
import numpy as np
from pyrate.core import config as cf
from tests.common import (
    assert_same_files_produced,
    assert_two_dirs_equal,
    manipulate_test_conf,
    TRAVIS,
    PYTHON3P6,
    PYTHON3P7,
    PYTHON3P8,
    GDAL_VERSION
)

# python3.7 and gdal3.0.4
REGRESSION = PYTHON3P7 and (GDAL_VERSION == '3.0.4')
# python3.7 and gdal3.0.2
REGRESSION2 = PYTHON3P7 and (GDAL_VERSION == '3.0.2')


@pytest.fixture(params=[0, 1])
def parallel(request):
    return request.param


@pytest.fixture(params=[1, 2, 4])
def local_crop(request):
    return request.param


@pytest.fixture()
def modified_config(tempdir, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees, ref_est_method):
    def modify_params(conf_file, parallel_vs_serial, output_conf_file):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)

        if params[cf.PROCESSOR] == 1:  # turn on coherence for gamma
            params[cf.COH_MASK] = 1

        params[cf.PARALLEL] = parallel_vs_serial
        params[cf.PROCESSES] = 4
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.REFNX], params[cf.REFNY] = 2, 2

        params[cf.IFG_CROP_OPT] = get_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["tiles"] = params["rows"] * params["cols"]

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        cf.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params
    return modify_params


@pytest.mark.slow
@pytest.mark.skipif(REGRESSION or PYTHON3P8 or PYTHON3P6, reason="Skip if not python3.7 and gdal=3.0.4")
def test_pipeline_parallel_vs_mpi(modified_config, gamma_conf):

    if np.random.randint(0, 1000) > 149:  # skip 85% of tests randomly
        pytest.skip("Randomly skipping as part of 85 percent")

    print("\n\n")
    print("===x==="*10)

    mpi_conf, params = modified_config(gamma_conf, 0, 'mpi_conf.conf')

    check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)

    try:
        check_call(f"mpirun -n 3 pyrate process -f {mpi_conf}", shell=True)
    except CalledProcessError as c:
        print(c)
        if TRAVIS:
            pytest.skip("Skipping as part of process error")
    check_call(f"mpirun -n 3 pyrate merge -f {mpi_conf}", shell=True)

    mr_conf, params_m = modified_config(gamma_conf, 1, 'multiprocess_conf.conf')

    check_call(f"pyrate workflow -f {mr_conf}", shell=True)

    sr_conf, params_s = modified_config(gamma_conf, 0, 'singleprocess_conf.conf')

    check_call(f"pyrate workflow -f {sr_conf}", shell=True)

    # convert2tif tests, 17 interferograms
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "*_unw_ifg.tif", 17)

    # if coherence masking, comprare coh files were converted
    if params[cf.COH_MASK]:
        assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "*_coh.tif", 17)
        print("coherence files compared")
        # 17 ifgs + 1 dem + 17 mlooked coh files
        assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR],
                                   f"*{params[cf.IFG_CROP_OPT]}cr.tif", 35)
    else:
        # 17 ifgs + 1 dem
        assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR],
                               f"*{params[cf.IFG_CROP_OPT]}cr.tif", 18)

    # prepifg + process steps that overwrite tifs test
    # ifg phase checking in the previous step checks the process pipeline upto APS correction

    # 2 x because of aps files
    assert_same_files_produced(params[cf.TMPDIR], params_m[cf.TMPDIR], params_s[cf.TMPDIR], "tsincr_*.npy",
                               params['tiles'] * 2)

    assert_same_files_produced(params[cf.TMPDIR], params_m[cf.TMPDIR], params_s[cf.TMPDIR], "tscuml_*.npy",
                               params['tiles'])

    assert_same_files_produced(params[cf.TMPDIR], params_m[cf.TMPDIR], params_s[cf.TMPDIR], "stack_rate_*.npy",
                               params['tiles'])
    assert_same_files_produced(params[cf.TMPDIR], params_m[cf.TMPDIR], params_s[cf.TMPDIR], "stack_error_*.npy",
                               params['tiles'])
    assert_same_files_produced(params[cf.TMPDIR], params_m[cf.TMPDIR], params_s[cf.TMPDIR], "stack_samples_*.npy",
                               params['tiles'])

    # compare merge step
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "stack*.tif", 3)
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "stack*.kml", 2)
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "stack*.png", 2)
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "stack*.npy", 3)
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "tscuml*.tif", 12)
    assert_same_files_produced(params[cf.OUT_DIR], params_m[cf.OUT_DIR], params_s[cf.OUT_DIR], "tsincr*.tif", 12)

    print("==========================xxx===========================")

    shutil.rmtree(params[cf.OBS_DIR])
    shutil.rmtree(params_m[cf.OBS_DIR])
    shutil.rmtree(params_s[cf.OBS_DIR])


@pytest.fixture(params=[0, 1])
def coh_mask(request):
    return request.param


@pytest.fixture()
def modified_config_short(tempdir, local_crop, get_lks, coh_mask, ref_pixel):
    orbfit_lks = 1
    orbfit_method = 1
    orbfit_degrees = 1
    ref_est_method = 1

    def modify_params(conf_file, parallel, output_conf_file, largetifs):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)
        params[cf.COH_MASK] = coh_mask
        params[cf.PARALLEL] = parallel
        params[cf.PROCESSES] = 4
        params[cf.APSEST] = 1
        params[cf.LARGE_TIFS] = largetifs
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.REFX], params[cf.REFY] = ref_pixel
        params[cf.REFNX], params[cf.REFNY] = 4, 4

        params[cf.IFG_CROP_OPT] = local_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["savenpy"] = 1
        params["tiles"] = params["rows"] * params["cols"]

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        cf.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.fixture
def create_mpi_files():

    def _create(modified_config_short, gamma_conf):

        mpi_conf, params = modified_config_short(gamma_conf, 0, 'mpi_conf.conf', 1)

        check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)

        try:
            check_call(f"mpirun -n 3 pyrate process -f {mpi_conf}", shell=True)
        except CalledProcessError as c:
            print(c)
            if TRAVIS:
                pytest.skip("Skipping as we encountered a process error")
        check_call(f"mpirun -n 3 pyrate merge -f {mpi_conf}", shell=True)
        return params

    return _create


@pytest.mark.slow
@pytest.mark.skipif(PYTHON3P8 or PYTHON3P7, reason="Skip if not python3.7 and gdal=3.0.4")
def test_stack_and_ts_mpi_vs_parallel_vs_serial(modified_config_short, gamma_conf, create_mpi_files, parallel):
    """
    Checks performed:
    1. mpi vs single process pipeline
    2. mpi vs parallel (python multiprocess) pipeline.
    3. Doing 1 and 2 means we have checked single vs parallel python multiprocess pipelines
    4. This also checks the entire pipeline using largetifs (new prepifg) vs old perpifg (python based)
    """
    if np.random.randint(0, 1000) > 399:  # skip 60% of tests randomly
        pytest.skip("Randomly skipping as part of 60 percent")

    print("\n\n")

    print("===x==="*10)

    params = create_mpi_files(modified_config_short, gamma_conf)

    sr_conf, params_p = modified_config_short(gamma_conf, parallel, 'parallel_conf.conf', 0)

    check_call(f"pyrate workflow -f {sr_conf}", shell=True)

    # convert2tif tests, 17 interferograms
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "*_unw_ifg.tif", 17)

    # if coherence masking, compare coh files were converted
    if params[cf.COH_MASK]:
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "*_coh.tif", 17)
        print("coherence files compared")

    # prepifg + process steps that overwrite tifs test
    # 17 mlooked ifgs + 1 dem + 17 mlooked coherence files
    if params[cf.COH_MASK]:
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], f"*{params[cf.IFG_CROP_OPT]}cr.tif", 35)
    else:
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], f"*{params[cf.IFG_CROP_OPT]}cr.tif", 18)

    # ifg phase checking in the previous step checks the process pipeline upto APS correction
    assert_two_dirs_equal(params[cf.TMPDIR], params_p[cf.TMPDIR], "tsincr_*.npy", params['tiles'] * 2)
    assert_two_dirs_equal(params[cf.TMPDIR], params_p[cf.TMPDIR], "tscuml_*.npy", params['tiles'])
    assert_two_dirs_equal(params[cf.TMPDIR], params_p[cf.TMPDIR], "stack_rate_*.npy", params['tiles'])
    assert_two_dirs_equal(params[cf.TMPDIR], params_p[cf.TMPDIR], "stack_error_*.npy", params['tiles'])
    assert_two_dirs_equal(params[cf.TMPDIR], params_p[cf.TMPDIR], "stack_samples_*.npy", params['tiles'])

    # compare merge step
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "stack*.tif", 3)
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "stack*.kml", 2)
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "stack*.png", 2)
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "stack*.npy", 3)
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "tscuml*.tif")

    print("==========================xxx===========================")

    shutil.rmtree(params[cf.OBS_DIR])
    shutil.rmtree(params_p[cf.OBS_DIR])
