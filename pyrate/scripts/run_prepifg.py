# -*- coding: utf-8 -*-
import sys, luigi
from pyrate.tasks.utils import pythonifyConfig
from pyrate.tasks.prepifg import PrepareInterferograms


def main():
    """
    :param config_file: config file to use. This provides a convenient way to
     use run_prepifg from within the module.
    :return:
    """
    usage = 'Usage: python run_prepifg.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print usage
        return
    raw_config_file = sys.argv[1]

    luigi.configuration.LuigiConfigParser.add_config_path(
        pythonifyConfig(raw_config_file))
    luigi.build([PrepareInterferograms()], local_scheduler=True)


if __name__ == '__main__':
    main()
