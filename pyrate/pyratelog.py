import logging
import sys
import traceback
import warnings

from pyrate import mpiops


def configure(verbosity):
    log = logging.getLogger("")
    log.setLevel(verbosity)
    ch = MPIStreamHandler()
    formatter = ElapsedFormatter()
    ch.setFormatter(formatter)
    log.addHandler(ch)


class MPIStreamHandler(logging.StreamHandler):
    """
    Only logs messages from Node 0
    """
    def emit(self, record):
        if mpiops.chunk_index == 0:
            super().emit(record)


class ElapsedFormatter():

    def format(self, record):
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