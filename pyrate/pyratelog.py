#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
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
import logging
import sys
import traceback
import warnings

from pyrate import mpiops


def configure(verbosity):
    """
    Function to configure logging properties
    :param verbosity: str
        one of ['DEBUG', 'INFO', 'WARNING', 'ERROR']
    """
    log = logging.getLogger("")
    log.setLevel(verbosity)
    stream = MPIStreamHandler()
    formatter = ElapsedFormatter()
    stream.setFormatter(formatter)
    log.addHandler(stream)


class MPIStreamHandler(logging.StreamHandler):
    """
    Only logs messages from Node 0
    """
    def emit(self, record):
        if mpiops.rank == 0:
            super(MPIStreamHandler, self).emit(record)


class ElapsedFormatter:
    """
    Convenience class for used in showing timestamps
    """
    # pylint: disable=too-few-public-methods
    # pylint: disable=no-self-use
    def format(self, record):
        """ time formatter """
        lvl = record.levelname
        name = record.name
        t = int(round(record.relativeCreated/1000.0))
        msg = record.getMessage()
        logstr = "+{}s {}:{} {}".format(t, name, lvl, msg)
        return logstr


def warn_with_traceback(message, category, filename, lineno, line=None):
    """
    copied from:
    http://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    traceback.print_stack()
    log = sys.stderr
    log.write(warnings.formatwarning(
        message, category, filename, lineno, line))
