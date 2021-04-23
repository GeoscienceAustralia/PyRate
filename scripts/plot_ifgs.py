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
"""
This Python module plots the input interferograms to the PyRate software
"""

import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import numpy as np
import pyrate.constants as C
from pyrate.core.logger import pyratelogger as log, configure_stage_log
from pyrate.core.shared import DEM, InputTypes
from pyrate.main import _params_from_conf


try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cmap = mpl.cm.Spectral
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

    parser.add_argument('-f', '--config_file', action="store", type=str, default=None,
                        help="Pass configuration file", required=True)

    args = parser.parse_args()

    params = _params_from_conf(args.config_file)

    configure_stage_log(args.verbosity, 'plot_ifgs', Path(params[C.OUT_DIR]).joinpath('pyrate.log.').as_posix())

    if args.verbosity:
        log.setLevel(args.verbosity)
        log.info("Verbosity set to " + str(args.verbosity) + ".")

    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    ifgs = params[C.INTERFEROGRAM_FILES]
    num_ifgs = len(ifgs)
    log.info(f'Plotting {num_ifgs} interferograms')

    plt_rows = np.int(np.sqrt(num_ifgs))
    plt_cols = num_ifgs//plt_rows
    if num_ifgs % plt_rows:
        plt_cols += 1

    fig = plt.figure(figsize=(12*plt_rows, 8*plt_cols))

    tot_plots = 1
    for p_r in range(plt_rows):
        for p_c in range(plt_cols):
            ax = fig.add_subplot(plt_rows, plt_cols, tot_plots)
            ifg_num = plt_cols * p_r + p_c
            m_path = ifgs[ifg_num]
            __plot_ifg(m_path, cmap, ax, num_ifgs)
            if tot_plots == num_ifgs:
                break
            tot_plots += 1

    f_name = 'ifg-phase-plot.png'
    plt.savefig(f_name)
    log.info(f'Ifg phase data is plotted in {Path(f_name).as_posix()}')


def __plot_ifg(m_path, cmap, ax, num_ifgs):
    if m_path.input_type == InputTypes.IFG:
        ifg = DEM(m_path.converted_path)
    else:
        raise AttributeError("Can only plot tifs")
    ifg.open()
    im = ax.imshow(ifg.data, cmap=cmap)
    text = ax.set_title(Path(ifg.data_path).stem)
    text.set_fontsize(min(20, int(num_ifgs / 2)))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ifg.close()


if __name__ == "__main__":
    main()
