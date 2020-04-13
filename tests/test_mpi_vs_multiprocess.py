import glob
import shutil
import pytest
from os.path import join, basename
from pathlib import Path
from subprocess import check_call
from pyrate.core import mpiops, config as cf
from pyrate.configuration import Configuration
from tests import common
from tests.common import copytree


@pytest.fixture
def coherence_file():
    def create_file():
        return None
    return create_file


@pytest.fixture()
def modified_config(tempdir, get_lks, get_crop, orbfit_lks, orbfit_method, orbfit_degrees, ref_est_method):
    def modify_params(conf_file, output_conf_file):
        params = cf.get_config_params(conf_file)

        tdir = Path(tempdir())
        copytree(params[cf.OBS_DIR], tdir)

        # manipulate params
        params[cf.OBS_DIR] = tdir.as_posix()
        params[cf.OUT_DIR] = tdir.joinpath('out').as_posix()
        params[cf.PARALLEL] = 1
        params[cf.APSEST] = 1
        params[cf.IFG_LKSX] = get_lks
        params[cf.DEM_FILE] = tdir.joinpath(Path(params[cf.DEM_FILE]).name).as_posix()
        params[cf.DEM_HEADER_FILE] = tdir.joinpath(Path(params[cf.DEM_HEADER_FILE]).name).as_posix()
        params[cf.SLC_FILE_LIST] = tdir.joinpath(Path(params[cf.SLC_FILE_LIST]).name).as_posix()
        params[cf.SLC_DIR] = tdir.as_posix()

        params[cf.IFG_CROP_OPT] = get_crop
        params[cf.ORBITAL_FIT_LOOKS_X], params[cf.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[cf.ORBITAL_FIT] = 1
        params[cf.ORBITAL_FIT_METHOD] = orbfit_method
        params[cf.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[cf.REF_EST_METHOD] = ref_est_method

        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        cf.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params
    return modify_params


def test_conv2tif_prepifg_parallel_vs_mpi(modified_config, roipac_or_gamma_conf):
    mpi_conf, params = modified_config(roipac_or_gamma_conf, 'mpi_conf.conf')

    check_call(f"mpirun -n 3 pyrate conv2tif -f {mpi_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {mpi_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate process -f {mpi_conf}", shell=True)

    sr_conf, params_s = modified_config(roipac_or_gamma_conf, 'multiprocess_conf.conf')

    check_call(f"pyrate conv2tif -f {sr_conf}", shell=True)
    check_call(f"pyrate prepifg -f {sr_conf}", shell=True)
    check_call(f"pyrate process -f {sr_conf}", shell=True)

    # convert2tif tests
    __assert_same_files_produced(params[cf.OBS_DIR], params_s[cf.OBS_DIR], "*_unw.tif", 17)

    # prepifg tests
    __assert_same_files_produced(params[cf.OUT_DIR], params_s[cf.OUT_DIR], "*cr.tif", 18)

    # TODO: timeseris and stack asserts

    shutil.rmtree(params[cf.OBS_DIR])
    shutil.rmtree(params_s[cf.OBS_DIR])


def __assert_same_files_produced(dir1, dir2, ext, num_files):
    full_res_mpi_rifs = glob.glob(join(dir1, ext))
    full_res_serial_tifs = glob.glob(join(dir2, ext))
    full_res_mpi_rifs.sort()
    full_res_serial_tifs.sort()
    # 17 unwrapped geotifs
    # 17 cropped multilooked tifs + 1 dem
    assert len(full_res_mpi_rifs) == len(full_res_serial_tifs) == num_files
    for m_f, s_f in zip(full_res_mpi_rifs, full_res_serial_tifs):
        assert basename(m_f) == basename(s_f)
        common.assert_tifs_equal(m_f, s_f)
