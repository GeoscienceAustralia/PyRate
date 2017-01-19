# coding: utf-8
# pylint: disable= invalid-name,  unused-import
"""For compatibility"""

from __future__ import absolute_import
import sys


PY3 = (sys.version_info[0] == 3)

if PY3:
    from io import StringIO
    import pickle
else:
    from StringIO import StringIO
    import cPickle as pickle

try:
    import pygrib
    import PyAPS
    PyAPS_INSTALLED = True
except ImportError:
    PyAPS_INSTALLED = False


class PyPASException(Exception):
    pass


def validate_sklearn():
    if not PyAPS_INSTALLED:
        raise PyPASException('PyAPS needs to be installed in order '
                             'to use this module')
