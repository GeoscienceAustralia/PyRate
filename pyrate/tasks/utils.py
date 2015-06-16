import os, pickle, luigi
from StringIO import StringIO
from pyrate import config
from pyrate.shared import Ifg, DEM, Incidence

DUMMY_SECTION_NAME = 'pyrate'
EXTENTS_FILE_NAME = 'region_extents.pkl'



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

    @property
    def extentsFileName(self):
        return os.path.join(self.obsDir, EXTENTS_FILE_NAME)



class DictParam(luigi.Parameter):
    def parse(self, string):
        sio = StringIO(string)
        return pickle.load(sio)

    def serialize(self, dct):
        sio = StringIO()
        pickle.dump(dct, sio)
        return sio.getvalue()



class RasterParam(DictParam):
    def parse(self, string):
        dct = super(RasterParam, self).parse(string)
        rasterType = dct['type']
        path = dct['path']

        if rasterType == 'DEM':
            return DEM(path)
        elif rasterType == 'Ifg':
            return Ifg(path)
        elif rasterType == 'Incidence':
            return Incidence(path)
        else:
            raise luigi.parameter.UnknownParameterException(
                'rasterBase must be an inscance DEM, Ifg or Incidence is valid')

    def serialize(self, rasterBase):
        path = rasterBase.data_path

        if isinstance(rasterBase, DEM):
            d = {'type':'DEM', 'path':path}
        elif isinstance(rasterBase, Ifg):
            d = {'type':'Ifg', 'path':path}
        elif isinstance(rasterBase, Incidence):
            d = {'type':'Incidence', 'path':path}
        else:
            raise luigi.parameter.UnknownParameterException(
                'rasterBase must be an inscance DEM, Ifg or Incidence is valid')

        return super(RasterParam, self).serialize(d)



def pythonifyConfig(configFile):
    outputFilename = '{}.python'.format(configFile)
    with open(configFile, 'r') as inputFile:
        with open(outputFilename, 'w') as outputFile:
            outputFile.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            for line in inputFile:
                outputFile.write(line)
    return outputFilename

