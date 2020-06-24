import shutil
import pytest
from pathlib import Path
from subprocess import check_call
import numpy as np
from pyrate.core import config as cf
from tests.common import (
    assert_two_dirs_equal,
    manipulate_test_conf,
    TRAVIS,
    PYTHON3P6,
    PYTHON3P7,
    GDAL_VERSION
)
# python3.7 and gdal3.0.4
REGRESSION = PYTHON3P7 and (GDAL_VERSION == '3.0.4')
# python3.7 and gdal3.0.2
REGRESSION2 = PYTHON3P7 and (GDAL_VERSION == '3.0.2')


@pytest.fixture(params=[1, 2, 3, 4])
def local_crop(request):
    return request.param


@pytest.fixture()
def modified_config_short(tempdir, local_crop, get_lks, coh_mask):
    orbfit_lks = 1
    orbfit_method = 1
    orbfit_degrees = 1
    ref_est_method = 1

    def modify_params(conf_file, parallel, output_conf_file):

        tdir = Path(tempdir())

        params = manipulate_test_conf(conf_file, tdir)
        params[cf.COH_MASK] = coh_mask
        params[cf.PARALLEL] = parallel
        params[cf.PROCESSES] = 4
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.REFNX], params[cf.REFNY] = 4, 4

        params[cf.IFG_CROP_OPT] = local_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["tiles"] = params["rows"] * params["cols"]

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        cf.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.fixture
def create_mpi_files(modified_config_short):

    def _create(conf):

        mpi_conf, params = modified_config_short(conf, 0, 'mpi_conf.conf')

        check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
        check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)

        return params  # don't need the reamining params
    return _create


@pytest.fixture()
def modified_config_largetifs(tempdir, local_crop, get_lks, coh_mask):
    orbfit_lks = 1
    orbfit_method = 1
    orbfit_degrees = 1
    ref_est_method = 1

    def modify_params(conf_file, parallel, output_conf_file):
        tdir = Path(tempdir())
        params = manipulate_test_conf(conf_file, tdir)
        params[cf.COH_MASK] = coh_mask
        params[cf.LARGE_TIFS] = 1
        params[cf.PARALLEL] = parallel
        params[cf.PROCESSES] = 4
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.REFNX], params[cf.REFNY] = 4, 4

        params[cf.IFG_CROP_OPT] = local_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["tiles"] = params["rows"] * params["cols"]

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        cf.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.mark.slow
@pytest.mark.skipif(PYTHON3P6 or PYTHON3P7, reason="Only run in python 3.8")
def test_prepifg_largetifs_vs_python(modified_config_largetifs, gamma_conf, create_mpi_files):

    print("\n\n")
    print("===x==="*10)
    if TRAVIS and np.random.randint(0, 1000) > 499:  # skip 50% of tests randomly
        pytest.skip("Randomly skipping as part of 50 percent")

    params = create_mpi_files(gamma_conf)
    sr_conf, params_p = modified_config_largetifs(gamma_conf, 1, 'parallel_conf.conf')
    check_call(f"pyrate conv2tif -f {sr_conf}", shell=True)
    check_call(f"pyrate prepifg -f {sr_conf}", shell=True)

    # convert2tif tests, 17 interferograms
    assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "*_unw_ifg.tif", 17)

    # if coherence masking, compare coh files were converted
    if params[cf.COH_MASK]:
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], "*_coh.tif", 17)
        print("coherence files compared")
        # 17 ifgs + 1 dem + 17 mlooked file
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], f"*{params[cf.IFG_CROP_OPT]}cr.tif", 35)
    else:
        # prepifg
        # 17 ifgs + 1 dem
        assert_two_dirs_equal(params[cf.OUT_DIR], params_p[cf.OUT_DIR], f"*{params[cf.IFG_CROP_OPT]}cr.tif", 18)


    print("==========================xxx===========================")

    shutil.rmtree(params[cf.OBS_DIR])
    shutil.rmtree(params_p[cf.OBS_DIR])
