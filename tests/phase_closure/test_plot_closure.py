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


import shutil
from pathlib import Path
from subprocess import check_call
import pytest

import pyrate.constants as C
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


@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif((not PLOT),
                    reason='skipped as plotting packages are missing')
def test_plot_closure(mexico_cropa_params):
    config = Configuration(MEXICO_CROPA_CONF)
    params = config.__dict__
    check_call(f"mpirun -n 3 pyrate prepifg -f {MEXICO_CROPA_CONF}", shell=True)

    correct._copy_mlooked(params)
    correct.__validate_correct_steps(params)
    # house keeping
    correct._update_params_with_tiles(params)
    correct._create_ifg_dict(params)
    params[C.REFX_FOUND], params[C.REFY_FOUND] = correct.ref_pixel_calc_wrapper(params)

    # run through the correct steps in user specified sequence
    for step in steps:
        if step == 'phase_closure':
            correct.correct_steps[step](params, config)
        else:
            correct.correct_steps[step](params)

    closure_plot_file = Path(config.phase_closure_dir).joinpath(f'closure_loops_iteration_1_fig_0.png')
    assert closure_plot_file.exists()
    shutil.rmtree(params[C.OUT_DIR], ignore_errors=True)
