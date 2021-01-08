import shutil
from pathlib import Path
from subprocess import check_call
import pytest
from pyrate.core import config as cf
from pyrate import correct
from pyrate.configuration import Configuration
from tests.common import MEXICO_CROPA_CONF


try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cmap = mpl.cm.Spectral
    PLOT = True
except ImportError as e:
    PLOT = False


steps = ['orbfit',  'refphase',  'phase_closure']


@pytest.mark.slow
@pytest.mark.skipif(not PLOT, reason='skipped as plotting packages are missing')
def test_plot_closure(mexico_cropa_params):
    config = Configuration(MEXICO_CROPA_CONF)
    params = config.__dict__
    check_call(f"mpirun -n 3 pyrate prepifg -f {MEXICO_CROPA_CONF}", shell=True)

    correct._copy_mlooked(params)
    correct.__validate_correct_steps(params)
    # house keeping
    correct._update_params_with_tiles(params)
    correct._create_ifg_dict(params)
    params[cf.REFX_FOUND], params[cf.REFY_FOUND] = correct.ref_pixel_calc_wrapper(params)

    # run through the correct steps in user specified sequence
    for step in steps:
        if step == 'phase_closure':
            correct.correct_steps[step](params, config)
        else:
            correct.correct_steps[step](params)

    closure_plot_file = Path(params[cf.OUT_DIR]).joinpath(f'sum_closure.png')
    assert closure_plot_file.exists()
    shutil.rmtree(params[cf.OUT_DIR], ignore_errors=True)
