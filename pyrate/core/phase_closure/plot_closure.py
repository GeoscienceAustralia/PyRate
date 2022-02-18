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


from typing import List
import numpy as np
from pathlib import Path

from pyrate.core.phase_closure.mst_closure import WeightedLoop
from pyrate.core.logger import pyratelogger as log
from pyrate.configuration import Configuration


def plot_closure(closure: np.ndarray, loops: List[WeightedLoop],
                    config: Configuration, thr: float, iteration: int):
    thr = thr * np.pi
    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        cmap = mpl.cm.Spectral
    except ImportError as e:
        log.warn(ImportError(e))
        log.warn("Required plotting packages are not found in environment. "
                 "Closure loop plot will not be generated!!!")
        return

    nrows, ncols, n_loops = closure.shape

    # 49 ifgs per fig
    plt_rows = 7 
    plt_cols = 7
    plots_per_fig = plt_rows * plt_cols
    n_figs = n_loops // plots_per_fig + (n_loops % plots_per_fig > 0)

    all_fig_plots = 1
    for fig_i in range(n_figs):
        fig = plt.figure(figsize=(9*plt_rows, 6*plt_cols))
        this_fig_plots = 0
        for p_r in range(plt_rows):
            for p_c in range(plt_cols):
                if all_fig_plots == n_loops + 1:   
                    break
                    
                ax = fig.add_subplot(plt_rows, plt_cols, plt_cols * p_r + p_c + 1)
                data = closure[:, :, plt_cols * p_r + p_c + fig_i * plots_per_fig]
                loop = loops[plt_cols * p_r + p_c + fig_i * plots_per_fig]
                title = ',\n'.join([repr(l) for l in loop.loop])
                im = ax.imshow(data, vmin=-thr, vmax=thr, cmap=cmap)
                plt.tick_params(axis='both', which='major', labelsize=12)
                text = ax.set_title(title)
                text.set_fontsize(16)
                #text.set_fontsize(min(20, int(n_loops/3)))

                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax)
                this_fig_plots += 1
                all_fig_plots += 1

        fig.tight_layout()

        closure_plot_file = Path(config.phase_closure_dir).joinpath(f'closure_loops_iteration_{iteration}_fig_{fig_i}.png')
        plt.savefig(closure_plot_file, dpi=100)
        plt.close(fig)
        log.info(f'{this_fig_plots} closure loops plotted in {closure_plot_file}')
