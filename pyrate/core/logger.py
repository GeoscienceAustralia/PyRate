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
This Python module contains functions to control PyRate log outputs
"""
import sys
import logging
import warnings
import traceback
from datetime import datetime
from pyrate.core.mpiops import size, rank, run_once


pyratelogger = logging.getLogger(__name__)
formatter = logging.Formatter("%(asctime)s %(module)s:%(lineno)d %(process)d %(levelname)s " + str(rank) + "/" + str(
    size-1)+" %(message)s", "%H:%M:%S")


def configure_stage_log(verbosity, step_name, log_file_name='pyrate.log.'):

    log_file_name = run_once(str.__add__, log_file_name, step_name + '.' + datetime.now().isoformat())

    ch = MPIStreamHandler()
    ch.setLevel(verbosity)
    ch.setFormatter(formatter)

    fh = logging.FileHandler(log_file_name)
    fh.setLevel(verbosity)
    fh.setFormatter(formatter)

    pyratelogger.addHandler(ch)
    pyratelogger.addHandler(fh)


def warn_with_traceback(message, category, filename, lineno, line=None):
    """
    copied from:
    http://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    traceback.print_stack()
    log = sys.stderr
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


class MPIStreamHandler(logging.StreamHandler):
    """
    Only logs messages from Node 0
    """
    def emit(self, record):
        if rank == 0:
            super().emit(record)
