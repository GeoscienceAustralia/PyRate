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
This Python module contains functions to control PyRate log outputs
"""
import sys
import logging
import warnings
import traceback
from datetime import datetime
from mpi4py import MPI
from os.path import abspath
from pyrate.core.mpiops import size, rank, run_once


pyratelogger = logging.getLogger(__name__)
formatter = logging.Formatter("%(asctime)s %(module)s %(lineno)d %(process)d %(levelname)s " + str(rank) + "/" + str(
    size-1)+" %(message)s", "%H:%M:%S")


class MPIFileHandler(logging.FileHandler):
    """
    stolen from: https://gist.github.com/chengdi123000/42ec8ed2cbef09ee050766c2f25498cb
    """
    def __init__(self,
                 filename,
                 mode=MPI.MODE_WRONLY | MPI.MODE_CREATE | MPI.MODE_APPEND,
                 encoding='utf-8',
                 delay=False,
                 comm=MPI.COMM_WORLD):
        self.baseFilename = abspath(filename)
        self.mode = mode
        self.encoding = encoding
        self.comm = comm
        if delay:
            # We don't open the stream, but we still need to call the
            # Handler constructor to set level, formatter, lock etc.
            logging.Handler.__init__(self)
            self.stream = None
        else:
            logging.StreamHandler.__init__(self, self._open())

    def _open(self):
        stream = MPI.File.Open(self.comm, self.baseFilename, self.mode)
        stream.Set_atomicity(True)
        return stream

    def emit(self, record):
        """
        Emit a record.

        If a formatter is specified, it is used to format the record.
        The record is then written to the stream with a trailing newline.  If
        exception information is present, it is formatted using
        traceback.print_exception and appended to the stream.  If the stream
        has an 'encoding' attribute, it is used to determine how to do the
        output to the stream.

        Modification:
            stream is MPI.File, so it must use `Write_shared` method rather
            than `write` method. And `Write_shared` method only accept
            bytestring, so `encode` is used. `Write_shared` should be invoked
            only once in each all of this emit function to keep atomicity.
        """
        try:
            msg = self.format(record)
            stream = self.stream
            stream.Write_shared((msg + self.terminator).encode(self.encoding))
            # self.flush()
        except Exception:
            self.handleError(record)

    def close(self):
        if self.stream:
            self.stream.Sync()
            self.stream.Close()
            self.stream = None


def configure_stage_log(verbosity, step_name, log_file_name='pyrate.log.'):

    log_file_name = run_once(str.__add__, log_file_name, step_name + '.' + datetime.now().isoformat())

    ch = logging.StreamHandler()
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
    log.write(warnings.formatwarning(
        message, category, filename, lineno, line))
