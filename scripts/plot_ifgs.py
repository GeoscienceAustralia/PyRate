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
This Python module plots the input interferograms to the PyRate software
"""

import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import math
import numpy as np
import pyrate.constants as C
from pyrate.core.logger import pyratelogger as log, configure_stage_log
from pyrate.core.shared import Ifg, InputTypes, nan_and_mm_convert
from pyrate.main import _params_from_conf


try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cmap = mpl.cm.Spectral_r
    cmap.set_bad(color='grey')
except ImportError as e:
    log.warn(ImportError(e))
    log.warn("Required plotting packages are not found in environment. "
             "Please install matplotlib in your environment to continue plotting!!!")
    raise ImportError(e)


def main():

    parser = argparse.ArgumentParser(prog='plot_ifgs', description="Python script to plot interferograms",
                                     add_help=True,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--verbosity', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help="Increase output verbosity")

    parser.add_argument('-d', '--directory', action="store", type=str, default=None,
                        help="Pass path to directory containing ifgs", required=True)

    parser.add_argument('-f', '--config_file', action="store", type=str, default=None,
                        help="Pass configuration file", required=True)

    parser.add_argument('-n', '--ifgs_per_plot', type=int, default=50, help='number of ifgs per plot', required=False)

    args = parser.parse_args()

    params = _params_from_conf(args.config_file)

    configure_stage_log(args.verbosity, 'plot_ifgs', Path(params[C.OUT_DIR]).joinpath('pyrate.log.').as_posix())

    if args.verbosity:
        log.setLevel(args.verbosity)
        log.info("Verbosity set to " + str(args.verbosity) + ".")

    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    ifgs = sorted(list(Path(args.directory).glob('*_ifg.tif')))

    num_ifgs = len(ifgs)
    if num_ifgs == 0:
        log.warning(f'No interferograms with extension *_ifg.tif were found in {args.directory}')
        return

    log.info(f'Plotting {num_ifgs} interferograms found in {args.directory}')

    ifgs_per_plot = args.ifgs_per_plot

    plt_rows = np.int(np.sqrt(ifgs_per_plot))
    plt_cols = ifgs_per_plot//plt_rows

    if ifgs_per_plot % plt_rows:
        plt_cols += 1

    tot_plots = 0
    num_of_figs = math.ceil(num_ifgs / ifgs_per_plot)

    fig_no = 0
    for i in range(num_of_figs):
        fig_no += 1
        fig = plt.figure(figsize=(12*plt_rows, 8*plt_cols))
        fig_plots = 0
        for p_r in range(plt_rows):
            for p_c in range(plt_cols):
                ax = fig.add_subplot(plt_rows, plt_cols, fig_plots + 1)
                ifg_num = plt_cols * p_r + p_c + ifgs_per_plot * i
                file = ifgs[ifg_num]
                log.info(f'Plotting {file}')
                __plot_ifg(file, cmap, ax, num_ifgs, params)
                tot_plots += 1
                fig_plots += 1
                log.debug(f'Plotted interferogram #{tot_plots}')
                if (fig_plots == ifgs_per_plot) or (tot_plots == num_ifgs):
                    f_name = Path(args.directory).joinpath('ifg-phase-plot-' + str(fig_no) + '.png').as_posix()
                    plt.savefig(f_name, dpi=50)
                    log.info(f'{fig_plots} interferograms plotted in {f_name}')
                    plt.close(fig)
                    break
            if tot_plots == num_ifgs:
                break


def __plot_ifg(file, cmap, ax, num_ifgs, params):
    try:
        ifg = Ifg(file)
        ifg.open()
    except:
        raise AttributeError(f'Cannot open interferogram geotiff: {file}')

    # change nodata values to NaN for display; convert units to mm (if in radians)
    nan_and_mm_convert(ifg, params)

    im = ax.imshow(ifg.phase_data, cmap=cmap)
    text = ax.set_title(Path(ifg.data_path).stem)
    text.set_fontsize(20)
#    text.set_fontsize(min(20, int(num_ifgs / 2)))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, label='mm')
    ifg.close()


if __name__ == "__main__":
    main()
