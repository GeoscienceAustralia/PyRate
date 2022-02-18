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
pyrate basic workflow for all supported input datasets

"""
import shutil
from subprocess import check_call
from pathlib import Path
import pytest
import numpy as np
import pyrate.constants as C
from pyrate.core.mpiops import MPI_INSTALLED
from pyrate.configuration import Configuration
from tests.common import MEXICO_CROPA_CONF, PY37GDAL302


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(not PY37GDAL302, reason="Only run in one CI env")
def test_workflow(system_conf):
    """check the handlers are working as expected"""
    check_call(f"mpirun -n 3 pyrate conv2tif -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate prepifg -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate correct -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate timeseries -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate stack -f {system_conf}", shell=True)
    check_call(f"mpirun -n 3 pyrate merge -f {system_conf}", shell=True)

    # assert logs generated in the outdir
    params = Configuration(system_conf).__dict__
    for stage in ['conv2tif', 'prepifg', 'correct', 'timeseries', 'stack', 'merge']:
        log_file_name = 'pyrate.log.' + stage
        files = list(Path(params[C.OUT_DIR]).glob(log_file_name + '.*'))
        assert len(files) == 1
    shutil.rmtree(params[C.OUT_DIR])


def test_single_workflow(gamma_or_mexicoa_conf):
    if MPI_INSTALLED:
        check_call(f"mpirun -n 4 pyrate workflow -f {gamma_or_mexicoa_conf}", shell=True)
    else:
        check_call(f"pyrate workflow -f {gamma_or_mexicoa_conf}", shell=True)

    params = Configuration(gamma_or_mexicoa_conf).__dict__

    log_file_name = 'pyrate.log.' + 'workflow'
    files = list(Path(params[C.OUT_DIR]).glob(log_file_name + '.*'))
    assert len(files) == 1

    # ref pixel file generated
    ref_pixel_file = params[C.REF_PIXEL_FILE]
    assert Path(ref_pixel_file).exists()
    ref_pixel = np.load(ref_pixel_file)
    if gamma_or_mexicoa_conf == MEXICO_CROPA_CONF:
        np.testing.assert_array_equal(ref_pixel, [42, 2])
        for f in C.GEOMETRY_OUTPUT_TYPES:
            assert Path(params[C.GEOMETRY_DIR]).joinpath(f + '.tif').exists()
    else:
        np.testing.assert_array_equal(ref_pixel, [38, 58])

    # assert orbfit exists on disc
    from pyrate.core import shared
    looked_files = [p.sampled_path for p in params[C.INTERFEROGRAM_FILES]]
    ifgs = [shared.Ifg(ifg) for ifg in looked_files]
    orbfits_on_disc = [Path(params[C.OUT_DIR], C.ORB_ERROR_DIR, Path(ifg.data_path).stem + '_orbfit.npy')
                       for ifg in ifgs]
    assert all(orbfits_on_disc)
    shutil.rmtree(params[C.OUT_DIR])
