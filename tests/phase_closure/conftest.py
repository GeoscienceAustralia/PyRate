#   This Python module is part of the PyRate software package.
#
#   Copyright 2021 Geoscience Australia
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


from pathlib import Path
import pytest
from pyrate import constants as C
from tests.phase_closure.common import IfgDummy
from tests.common import MEXICO_CROPA_DIR
from pyrate.constants import PYRATEPATH

GEOTIFF = PYRATEPATH.joinpath('tests', 'test_data', 'geotiffs')


@pytest.fixture()
def closure_params(geotiffs):
    ifg_files = [IfgDummy(ifg_path) for ifg_path in geotiffs]
    params = {
        C.INTERFEROGRAM_FILES: ifg_files,
        C.MAX_LOOP_LENGTH: 100,
        'geotiffs': geotiffs
    }
    return params


@pytest.fixture(scope='module')
def cropa_geotifs():
    tifs = [u.as_posix() for u in Path(MEXICO_CROPA_DIR).glob('*_unw.tif')]
    tifs.sort()
    return tifs

@pytest.fixture
def geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs


@pytest.fixture
def ten_geotiffs():
    tifs = [u.as_posix() for u in GEOTIFF.glob('*_unw.tif')]
    tifs.sort()
    return tifs[:10]

