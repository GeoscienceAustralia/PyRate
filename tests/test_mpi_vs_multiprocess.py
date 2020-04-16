import os
import shutil
import pytest
from pathlib import Path
from subprocess import check_call, check_output
import numpy as np
from pyrate.core import mpiops, config as cf
from tests import common
from tests.common import copytree


TRAVIS = True if 'TRAVIS' in os.environ else False
PYTHON3P6 = True if ('TRAVIS_PYTHON_VERSION' in os.environ and os.environ['TRAVIS_PYTHON_VERSION'] == '3.6') else False
GDAL_VERSION = check_output(["gdal-config", "--version"]).decode(encoding="utf-8").split('\n')[0]
REGRESSION = PYTHON3P6 or (TRAVIS and ((GDAL_VERSION == '3.0.4') or (GDAL_VERSION == '2.4.2')))


@pytest.fixture()
def modified_config(tempdir, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees, ref_est_method):
    def modify_params(conf_file, output_conf_file):
        params = cf.get_config_params(conf_file)

        tdir = Path(tempdir())
        copytree(params[cf.OBS_DIR], tdir)

        # manipulate params
        params[cf.OBS_DIR] = tdir.as_posix()
        outdir = tdir.joinpath('out')
        outdir.mkdir(exist_ok=True)
        params[cf.OUT_DIR] = outdir.as_posix()

        params[cf.PARALLEL] = 1
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX], params[cf.IFG_LKSY] = get_lks, get_lks
        params[cf.DEM_FILE] = tdir.joinpath(Path(params[cf.DEM_FILE]).name).as_posix()
        params[cf.DEM_HEADER_FILE] = tdir.joinpath(Path(params[cf.DEM_HEADER_FILE]).name).as_posix()
        params[cf.SLC_FILE_LIST] = tdir.joinpath(Path(params[cf.SLC_FILE_LIST]).name).as_posix()
        params[cf.SLC_DIR] = tdir.as_posix()
        params[cf.IFG_FILE_LIST] = tdir.joinpath(Path(params[cf.IFG_FILE_LIST]).name).as_posix()
        params[cf.COH_FILE_DIR] = tdir.as_posix()
        params[cf.APS_INCIDENCE_MAP] = tdir.joinpath(Path(params[cf.APS_INCIDENCE_MAP]).name).as_posix()
        params[cf.REFNX], params[cf.REFNY] = 4, 4
        params[cf.TMPDIR] = tdir.joinpath(Path(params[cf.TMPDIR]).name).as_posix()

        params[cf.IFG_CROP_OPT] = get_crop
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


@pytest.mark.skipif(REGRESSION, reason='skipping regression tests in travis except GDAL 3.0.2')
def test_conv2tif_prepifg_parallel_vs_mpi(modified_config, roipac_or_gamma_conf):

    BOOL = np.random.randint(0, 10) > 0
    if BOOL:
        pytest.skip("Skipping as part of 90")

    print("\n\n")
    print("===x==="*10)

    mpi_conf, params = modified_config(roipac_or_gamma_conf, 'mpi_conf.conf')

    check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)
    # check_call(f"mpirun -n 3 pyrate process -f {mpi_conf}", shell=True)

    sr_conf, params_s = modified_config(roipac_or_gamma_conf, 'multiprocess_conf.conf')

    check_call(f"pyrate conv2tif -f {sr_conf}", shell=True)
    check_call(f"pyrate prepifg -f {sr_conf}", shell=True)
    # check_call(f"pyrate process -f {sr_conf}", shell=True)
    # TODO: merge step

    # convert2tif tests, 17 interferograms
    __assert_same_files_produced(params[cf.OUT_DIR], params_s[cf.OUT_DIR], "*_unw.tif", 17)

    # prepifg + process steps that overwrite tifs test

    # 17 ifgs + 1 dem
    __assert_same_files_produced(params[cf.OUT_DIR], params_s[cf.OUT_DIR], f"*{params[cf.IFG_CROP_OPT]}cr.tif", 18)

    # ifg phase checking in the previous step checks the process pipeline upto APS correction

    # 2 x because of aps files
    # __assert_same_files_produced(params[cf.TMPDIR], params_s[cf.TMPDIR], "tsincr_*.npy", params['tiles'] * 2)
    #
    # __assert_same_files_produced(params[cf.TMPDIR], params_s[cf.TMPDIR], "tscuml_*.npy", params['tiles'])
    #
    # __assert_same_files_produced(params[cf.TMPDIR], params_s[cf.TMPDIR], "stack_rate_*.npy", params['tiles'])
    # __assert_same_files_produced(params[cf.TMPDIR], params_s[cf.TMPDIR], "stack_error_*.npy", params['tiles'])
    # __assert_same_files_produced(params[cf.TMPDIR], params_s[cf.TMPDIR], "stack_samples_*.npy", params['tiles'])
    print("==========================xxx===========================")

    shutil.rmtree(params[cf.OBS_DIR])
    shutil.rmtree(params_s[cf.OBS_DIR])


def __assert_same_files_produced(dir1, dir2, ext, num_files):
    mpi_files = list(Path(dir1).glob(ext))
    mp_files = list(Path(dir2).glob(ext))  # MultiProcess files
    mpi_files.sort()
    mp_files.sort()
    print(mpi_files)
    print(mp_files)
    # 17 unwrapped geotifs
    # 17 cropped multilooked tifs + 1 dem
    assert len(mpi_files) == num_files
    assert len(mp_files) == num_files
    if mpi_files[0].suffix == '.tif':
        for m_f, s_f in zip(mpi_files, mp_files):
            assert m_f.name == s_f.name
            common.assert_tifs_equal(m_f.as_posix(), s_f.as_posix())
    elif mpi_files[0].suffix == '.npy':
        for m_f, s_f in zip(mpi_files, mp_files):
            print(m_f.name, s_f.name)
            assert m_f.name == s_f.name
            np.testing.assert_array_almost_equal(np.load(m_f), np.load(s_f))
    else:
        raise
