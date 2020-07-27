# coding: utf-8
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
pyrate basic workflow for all supported input datasets

"""
from subprocess import check_call
from pathlib import Path
import pytest
import numpy as np
from pyrate.core import config as cf
from pyrate.configuration import Configuration


@pytest.mark.slow
def test_workflow(system_conf):
    """check the handlers are working as expected"""
    check_call(f"mpirun -n 3 pyrate conv2tif -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate process -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate merge -f {system_conf}", shell=True)

    # assert logs generated in the outdir
    params = Configuration(system_conf).__dict__
    for stage in ['conv2tif', 'prepifg', 'process', 'merge']:
        log_file_name = 'pyrate.log.' + stage
        files = list(Path(params[cf.OUT_DIR]).glob(log_file_name + '.*'))
        assert len(files) == 1


@pytest.mark.slow
def test_single_workflow(gamma_conf):

    check_call(f"mpirun -n 4 pyrate workflow -f {gamma_conf}", shell=True)

    params = Configuration(gamma_conf).__dict__

    log_file_name = 'pyrate.log.' + 'workflow'
    files = list(Path(params[cf.OUT_DIR]).glob(log_file_name + '.*'))
    assert len(files) == 1

    # ref pixel file generated
    ref_pixel_file = Path(params[cf.OUT_DIR]).joinpath(cf.REF_PIXEL_FILE)
    assert ref_pixel_file.exists()
    ref_pixel = np.load(ref_pixel_file)
    np.testing.assert_array_equal(ref_pixel, [38, 58])

    # assert orbfit exists on disc
    from pyrate.core import shared
    looked_files = [p.sampled_path for p in params[cf.INTERFEROGRAM_FILES]]
    ifgs = [shared.Ifg(ifg) for ifg in looked_files]
    orbfits_on_disc = [Path(params[cf.OUT_DIR], cf.ORB_ERROR_DIR,
                          Path(ifg.data_path).stem + '_orbfit.npy')
                       for ifg in ifgs]
    assert all(orbfits_on_disc)
