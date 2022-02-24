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
This Python module contains tests that compare system vs python prepifg methods.
"""
import shutil
import pytest
from pathlib import Path
import numpy as np

import pyrate.configuration
import pyrate.constants as C
from pyrate import conv2tif, prepifg
from pyrate.configuration import Configuration
from tests.common import (
    assert_two_dirs_equal,
    manipulate_test_conf,
    GITHUB_ACTIONS,
    PY37GDAL302
)


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
        params[C.COH_MASK] = coh_mask
        params[C.PARALLEL] = parallel
        params[C.PROCESSES] = 4
        params[C.APSEST] = 1
        params[C.IFG_LKSX], params[C.IFG_LKSY] = get_lks, get_lks
        params[C.REFNX], params[C.REFNY] = 4, 4

        params[C.IFG_CROP_OPT] = local_crop
        params[C.ORBITAL_FIT_LOOKS_X], params[
            C.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[C.ORBITAL_FIT] = 1
        params[C.ORBITAL_FIT_METHOD] = orbfit_method
        params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[C.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        pyrate.configuration.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.fixture
def create_mpi_files(modified_config_short):

    def _create(conf):

        mpi_conf, params = modified_config_short(conf, 0, 'mpi_conf.conf')
        params = Configuration(mpi_conf).__dict__
        conv2tif.main(params)
        params = Configuration(mpi_conf).__dict__
        prepifg.main(params)

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
        params[C.COH_MASK] = coh_mask
        params[C.LARGE_TIFS] = 1
        params[C.PARALLEL] = parallel
        params[C.PROCESSES] = 4
        params[C.APSEST] = 1
        params[C.IFG_LKSX], params[C.IFG_LKSY] = get_lks, get_lks
        params[C.REFNX], params[C.REFNY] = 4, 4

        params[C.IFG_CROP_OPT] = local_crop
        params[C.ORBITAL_FIT_LOOKS_X], params[
            C.ORBITAL_FIT_LOOKS_Y] = orbfit_lks, orbfit_lks
        params[C.ORBITAL_FIT] = 1
        params[C.ORBITAL_FIT_METHOD] = orbfit_method
        params[C.ORBITAL_FIT_DEGREE] = orbfit_degrees
        params[C.REF_EST_METHOD] = ref_est_method
        params["rows"], params["cols"] = 3, 2
        params["notiles"] = params["rows"] * params["cols"]  # number of tiles

        print(params)
        # write new temp config
        output_conf = tdir.joinpath(output_conf_file)
        pyrate.configuration.write_config_file(params=params, output_conf_file=output_conf)

        return output_conf, params

    return modify_params


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL302, reason="Only run in one CI env")
def test_prepifg_largetifs_vs_python(modified_config_largetifs, gamma_conf, create_mpi_files):

    print("\n\n")
    print("===x==="*10)
    if GITHUB_ACTIONS and np.random.randint(0, 1000) > 899:  # skip 90% of tests randomly
        pytest.skip("Randomly skipping as part of 50 percent")

    params = create_mpi_files(gamma_conf)
    sr_conf, params_p = modified_config_largetifs(gamma_conf, 1, 'parallel_conf.conf')
    params_p = Configuration(sr_conf).__dict__
    conv2tif.main(params_p)
    params_p = Configuration(sr_conf).__dict__
    prepifg.main(params_p)
    params_p = Configuration(sr_conf).__dict__
    # convert2tif tests, 17 interferograms
    assert_two_dirs_equal(params[C.INTERFEROGRAM_DIR], params_p[C.INTERFEROGRAM_DIR], "*_unw.tif", 17)

    # if coherence masking, compare coh files were converted
    if params[C.COH_FILE_LIST] is not None:
        assert_two_dirs_equal(params[C.COHERENCE_DIR], params_p[C.COHERENCE_DIR], "*_cc.tif", 17)
        # 17 ifgs + 1 dem + 17 mlooked file
        assert_two_dirs_equal(params[C.COHERENCE_DIR], params_p[C.COHERENCE_DIR], "*_coh.tif", 17)
        assert_two_dirs_equal(params[C.COHERENCE_DIR], params_p[C.COHERENCE_DIR], ["coh_*.tif"], 3)

    assert_two_dirs_equal(params[C.GEOMETRY_DIR], params_p[C.GEOMETRY_DIR], "*_dem.tif", 1)
    assert_two_dirs_equal(params[C.GEOMETRY_DIR], params_p[C.GEOMETRY_DIR],
                          [t + '.tif' for t in C.GEOMETRY_OUTPUT_TYPES], 6)
    # prepifg
    # 17 ifgs + 1 dem
    assert_two_dirs_equal(params[C.INTERFEROGRAM_DIR], params_p[C.INTERFEROGRAM_DIR], "*_ifg.tif", 17)
    assert_two_dirs_equal(params[C.GEOMETRY_DIR], params_p[C.GEOMETRY_DIR], "dem.tif", 1)

    print("==========================xxx===========================")
    shutil.rmtree(params[C.OUT_DIR])
    shutil.rmtree(params_p[C.OUT_DIR])
