'''
Script to convert Gamma headers to ESRI's BIL format.

.. codeauthor:: Ben Davies, NCI
'''

import sys, luigi
from pyrate.tasks.utils import pythonifyConfig
from pyrate.tasks import ConvertToGeotiff


def main():
    usage = 'Usage: gamma.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print usage
        return
    rawConfigFile = sys.argv[1]
    luigi.configuration.LuigiConfigParser.add_config_path(pythonifyConfig(rawConfigFile))
    luigi.build([ConvertToGeotiff()], local_scheduler=True)


if __name__ == "__main__":
    main()
