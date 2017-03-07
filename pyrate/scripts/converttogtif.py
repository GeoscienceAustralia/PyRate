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
This Python module uses Luigi to convert interferograms to geotiff format
"""
from __future__ import print_function
import sys
import luigi
from pyrate.tasks.utils import pythonify_config
from pyrate.tasks import ConvertToGeotiff


def main():
    """
    Wrapper function to convert unwrapped interferrograms into geotifs
    for all data types
    """
    usage = 'Usage: gamma.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' \
            or sys.argv[1] == '--help':  # pragma: no cover
        print(usage)
        return
    raw_config_file = sys.argv[1]     # this does '~' expansion automatically
    luigi.configuration.LuigiConfigParser.add_config_path(
        pythonify_config(raw_config_file))
    luigi.build([ConvertToGeotiff()], local_scheduler=True)


if __name__ == "__main__":
    main()
