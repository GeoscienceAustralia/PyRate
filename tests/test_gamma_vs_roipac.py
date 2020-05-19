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
This Python module contains tests that compare GAMMA and ROI_PAC
functionality in PyRate.
"""
import os
import shutil
import pytest
from pathlib import Path

from pyrate.core.shared import DEM
from pyrate.core import ifgconstants as ifc, config as cf
from pyrate.core.prepifg_helper import _is_number
from pyrate import prepifg, conv2tif, configuration
from tests.common import SML_TEST_DIR, small_data_setup, copytree, TEST_CONF_ROIPAC, TEST_CONF_GAMMA


SMLNEY_GAMMA_TEST = os.path.join(SML_TEST_DIR, "gamma_obs")


def test_files_are_same(tempdir, get_config):
    roipac_params = get_config(TEST_CONF_ROIPAC)
    roipac_tdir = Path(tempdir())
    roipac_params = __workflow(roipac_params, roipac_tdir)

    gamma_params = get_config(TEST_CONF_GAMMA)
    gamma_tdir = Path(tempdir())
    gamma_params = __workflow(gamma_params, gamma_tdir)

    # conv2tif output equal
    __assert_same_files_produced(roipac_params[cf.OUT_DIR], gamma_params[cf.OUT_DIR], "*_unw.tif", 17)

    # prepifg output equal
    __assert_same_files_produced(roipac_params[cf.OUT_DIR], gamma_params[cf.OUT_DIR],
                                 f"*{roipac_params[cf.IFG_CROP_OPT]}cr.tif", 18)

    # clean up
    shutil.rmtree(roipac_params[cf.OBS_DIR])
    shutil.rmtree(gamma_params[cf.OBS_DIR])


def __workflow(params, tdir):
    copytree(params[cf.OBS_DIR], tdir)
    # manipulate params
    params[cf.OBS_DIR] = tdir.as_posix()
    outdir = tdir.joinpath('out')
    outdir.mkdir(exist_ok=True)
    params[cf.OUT_DIR] = outdir.as_posix()

    params[cf.DEM_FILE] = tdir.joinpath(Path(params[cf.DEM_FILE]).name).as_posix()
    params[cf.DEM_HEADER_FILE] = tdir.joinpath(Path(params[cf.DEM_HEADER_FILE]).name).as_posix()
    params[cf.HDR_FILE_LIST] = tdir.joinpath(Path(params[cf.HDR_FILE_LIST]).name).as_posix()
    params[cf.SLC_DIR] = tdir.as_posix()
    params[cf.IFG_FILE_LIST] = tdir.joinpath(Path(params[cf.IFG_FILE_LIST]).name).as_posix()
    params[cf.COH_FILE_DIR] = tdir.as_posix()
    params[cf.APS_INCIDENCE_MAP] = tdir.joinpath(Path(params[cf.APS_INCIDENCE_MAP]).name).as_posix()
    params[cf.TMPDIR] = tdir.joinpath(Path(params[cf.TMPDIR]).name).as_posix()
    output_conf = tdir.joinpath('roipac_temp.conf')
    cf.write_config_file(params=params, output_conf_file=output_conf)
    params = configuration.Configuration(output_conf).__dict__
    conv2tif.main(params)
    prepifg.main(params)
    return params


def __assert_same_files_produced(dir1, dir2, ext, num_files):
    dir1_files = list(Path(dir1).glob(ext))
    dir2_files = list(Path(dir2).glob(ext))
    dir1_files.sort()
    dir2_files.sort()
    # 17 unwrapped geotifs
    # 17 cropped multilooked tifs + 1 dem
    assert len(dir1_files) == num_files
    assert len(dir2_files) == num_files
    c = 0

    all_roipac_ifgs = [f for f in small_data_setup(dir1_files) if not isinstance(f, DEM)]
    all_gamma_ifgs = [f for f in small_data_setup(dir2_files) if not isinstance(f, DEM)]

    for c, (i, j) in enumerate(zip(all_roipac_ifgs, all_gamma_ifgs)):
        mdi = i.meta_data
        mdj = j.meta_data
        for k in mdi:  # all key values equal
            if k == "INCIDENCE_DEGREES":
                pass  # incidence angle not implemented for roipac
            elif _is_number(mdi[k]):
                assert pytest.approx(float(mdj[k]), 0.00001) == float(mdi[k])
            elif mdi[k] == "ROIPAC" or "GAMMA":
                pass  # INSAR_PROCESSOR can not be equal
            else:
                assert mdj[k] == mdi[k]

        if i.data_path.__contains__("_{looks}rlks_{crop}cr".format(looks=1, crop=1)):
            # these are multilooked tifs
            # test that DATA_STEP is MULTILOOKED
            assert mdi[ifc.DATA_TYPE] == ifc.MULTILOOKED
            assert mdj[ifc.DATA_TYPE] == ifc.MULTILOOKED
        else:
            assert mdi[ifc.DATA_TYPE] == ifc.ORIG
            assert mdj[ifc.DATA_TYPE] == ifc.ORIG

    assert c + 1 == len(all_gamma_ifgs)
