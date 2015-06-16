import os, luigi
from pyrate import config

DUMMY_SECTION_NAME = 'pyrate'



class InputParam(dict):
    def __init__(self, name):
        self['section'] = DUMMY_SECTION_NAME
        self['name'] = name



class IfgListMixin(object):
    """
    Mixin to aid access to commonly used 'complex' valuea.
    """

    ifgListFile = luigi.Parameter(config_path=InputParam(config.IFG_FILE_LIST))
    obsDir = luigi.Parameter(config_path=InputParam(config.OBS_DIR))

    def ifgList(self, tif=True):
        fileNames = config.parse_namelist(self.ifgListFile)

        if tif:
            fileNames = ['%s.tif' % os.path.splitext(os.path.basename(fn))[0] for fn in fileNames]

        obsDir = self.obsDir
        if obsDir:
            fileNames = [os.path.join(obsDir, fn) for fn in fileNames]

        return fileNames



def pythonifyConfig(configFile):
    outputFilename = '{}.python'.format(configFile)
    with open(configFile, 'r') as inputFile:
        with open(outputFilename, 'w') as outputFile:
            outputFile.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            for line in inputFile:
                outputFile.write(line)
    return outputFilename

