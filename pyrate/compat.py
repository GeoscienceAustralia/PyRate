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
This Python module tests compatibilities
"""
# coding: utf-8
# pylint: disable= invalid-name,  unused-import

from __future__ import absolute_import
import sys


PY3 = (sys.version_info[0] == 3)


try:
    import pygrib
    import PyAPS
    PyAPS_INSTALLED = True
except ImportError:
    PyAPS_INSTALLED = False


class PyAPSException(Exception):
    """
    Convenience class for PyAPS status
    """


def validate_sklearn():
    """
    Convenience function validating PyAPS status
    """
    if not PyAPS_INSTALLED:
        raise PyAPSException('PyAPS needs to be installed in order '
                             'to use this module')
