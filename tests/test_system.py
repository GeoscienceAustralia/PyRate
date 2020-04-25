# coding: utf-8
#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
This Python module contains tests for mpi operations in PyRate.
Tun this module as 'mpirun -n 4 pytest tests/test_mpi.py'
"""
from subprocess import check_call
from pathlib import Path
from pyrate.core import config as cf
from pyrate.configuration import Configuration


def test_workflow(system_conf):
    """check the handlers are working as expected"""
    check_call(f"pyrate conv2tif -f {system_conf}", shell=True)
    check_call(f"pyrate prepifg -f {system_conf}", shell=True)
    check_call(f"pyrate process -f {system_conf}", shell=True)
    check_call(f"pyrate merge -f {system_conf}", shell=True)

    # assert logs generated in the outdir
    params = Configuration(system_conf).__dict__
    for stage in ['conv2tif', 'prepifg', 'process', 'merge']:
        log_file_name = 'pyrate.log.' + stage
        files = list(Path(params[cf.OUT_DIR]).glob(log_file_name + '.*'))
        assert len(files) == 1
