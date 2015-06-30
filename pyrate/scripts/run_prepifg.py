# -*- coding: utf-8 -*-





import sys, luigi
from pyrate.tasks.utils import pythonifyConfig
from pyrate.tasks.prepifg import PrepareInterferograms



def main():
    usage = 'Usage: python run_prepifg.py <config file>'
    if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print usage
        return
    rawConfigFile = sys.argv[1]
    luigi.configuration.LuigiConfigParser.add_config_path(pythonifyConfig(rawConfigFile))
    luigi.build([PrepareInterferograms()], local_scheduler=True)



if __name__ == '__main__':
    main()
