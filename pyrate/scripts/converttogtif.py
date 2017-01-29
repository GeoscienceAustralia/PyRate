from __future__ import print_function
"""
Script to convert Gamma headers to ESRI's BIL format.
"""

import sys
import luigi
from pyrate.tasks.utils import pythonify_config
from pyrate.tasks import ConvertToGeotiff


def main():
    usage = 'Usage: gamma.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(usage)
        return
    rawConfigFile = sys.argv[1]     # this does '~' expansion automatically
    luigi.configuration.LuigiConfigParser.add_config_path(
        pythonify_config(rawConfigFile))
    luigi.build([ConvertToGeotiff()], local_scheduler=True)


if __name__ == "__main__":
    main()
