import pytest
from pyrate import constants as C
from tests.phase_closure.common import IfgDummy


@pytest.fixture()
def closure_params(geotiffs):
    ifg_files = [IfgDummy(ifg_path) for ifg_path in geotiffs]
    params = {
        C.INTERFEROGRAM_FILES: ifg_files,
        C.MAX_LOOP_LENGTH: 100,
        'geotiffs': geotiffs
    }
    return params
