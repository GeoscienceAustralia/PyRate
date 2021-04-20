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
This Python module defines executable run configuration for the PyRate software
"""

import os
import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import numpy as np
import pyrate.constants as C
from pyrate.constants import CLI_DESCRIPTION
from pyrate.core.logger import pyratelogger as log, configure_stage_log
from pyrate.core import mpiops
from pyrate.core.shared import Ifg, InputTypes
from pyrate.configuration import Configuration


def _params_from_conf(config_file):
    config_file = os.path.abspath(config_file)
    config = Configuration(config_file)
    params = config.__dict__
    return params


def main():

    parser = argparse.ArgumentParser(prog='pyrate', description=CLI_DESCRIPTION, add_help=True,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--verbosity', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help="Increase output verbosity")

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    parser_plot = subparsers.add_parser('plot_ifgs', help='Plot interferogram.', add_help=True)
    parser_plot.add_argument('-f', '--config_file', action="store", type=str, default=None,
                                help="Pass configuration file", required=True)

    args = parser.parse_args()

    params = mpiops.run_once(_params_from_conf, args.config_file)

    configure_stage_log(args.verbosity, args.command, Path(params[C.OUT_DIR]).joinpath('pyrate.log.').as_posix())

    log.debug("Plotting")
    log.debug("Arguments supplied at command line: ")
    log.debug(args)

    try:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        cmap = mpl.cm.Spectral
    except ImportError as e:
        log.warn(ImportError(e))
        log.warn("Required plotting packages are not found in environment. "
                 "Please install matplotlib in your environment to continue plotting!!!")
        return

    ifgs = params[C.INTERFEROGRAM_FILES]
    num_ifgs = len(ifgs)

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
            if m_path.input_type == InputTypes.IFG:
                ifg = Ifg(ifgs[ifg_num].converted_path)
            else:
                raise AttributeError("Can only plot tifs")
            ifg.open()
            im = ax.imshow(ifg.phase_data, cmap=cmap)
            text = ax.set_title(Path(ifg.data_path).stem)
            text.set_fontsize(min(20, int(num_ifgs/3)))

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            if tot_plots == num_ifgs:
                break
            tot_plots += 1

    plt.show()


if __name__ == "__main__":
    main()
