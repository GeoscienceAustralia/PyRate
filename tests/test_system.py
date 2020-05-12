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
import pytest
from subprocess import check_call
from pathlib import Path
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
