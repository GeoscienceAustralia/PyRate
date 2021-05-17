from pathlib import Path
import pytest
from pyrate import constants as C
from tests.phase_closure.common import IfgDummy
from tests.common import MEXICO_CROPA_DIR


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
