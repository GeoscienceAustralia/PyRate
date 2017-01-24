import os
import pickle
import luigi
from io import StringIO, BytesIO
from pyrate import config
from pyrate.config import (OBS_DIR,
                           IFG_FILE_LIST,
                           DEM_FILE,
                           DEM_HEADER_FILE,
                           OUT_DIR,
                           ROIPAC_RESOURCE_HEADER)
from pyrate.shared import Ifg, DEM, Incidence

DUMMY_SECTION_NAME = 'pyrate'
EXTENTS_FILE_NAME = 'region_extents.pkl'


class InputParam(dict):
    '''
    Convenience class for specifying the parameter in a PyRate configuration
    file, intended for use in specifying which value in a PyRate configuration
    file should be used for a :py:class:`luigi.Parameter` (which is specified
    by the parameter *conf_path* to its constructor.)
    '''

    def __init__(self, name):
        self['section'] = DUMMY_SECTION_NAME
        self['name'] = name


class IfgListMixin(object):
    """
    Mixin to aid access to commonly used computed values from the PyRate
    config file.

    .. todo:: This should perhaps be renamed to something like *ConfigMixin*
        for clarity, as it is ued for accessing more than the list of
        interferograms.
    """

    ifgListFile = luigi.Parameter(config_path=InputParam(config.IFG_FILE_LIST))
    obsDir = luigi.Parameter(config_path=InputParam(config.OBS_DIR))
    out_dir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))

    def ifgList(self, tif=True):
        """
        Get the list of interferograms to process.

        :param tif: Should the tif files be returned (*True*) or the raw
            interferograms (*False*). The latter will probably only be required
            before conversion to geotif files occurs.
        """

        file_names = config.parse_namelist(self.ifgListFile)

        if tif:
            file_names = ['%s.tif' %
                          os.path.splitext(os.path.basename(fn))[0]
                          for fn in file_names]

        obs_dir = self.obsDir
        if obs_dir:
            file_names = [os.path.join(obs_dir, fn) for fn in file_names]

        return file_names

    def ifgTiffList(self, tif=True):
        """
        Get the list of interferograms to process.

        :param tif: Should the tif files be returned (*True*) or the raw
            interferograms (*False*). The latter will probably only be required
            before conversion to geotif files occurs.
        """

        file_names = config.parse_namelist(self.ifgListFile)

        if tif:
            file_names = ['%s.tif' %
                          os.path.splitext(os.path.basename(fn))[0]
                          for fn in file_names]

        out_dir = self.out_dir
        if out_dir:
            file_names = [os.path.join(out_dir, fn) for fn in file_names]

        return file_names

    @property
    def extentsFileName(self):
        return os.path.join(self.out_dir, EXTENTS_FILE_NAME)


class DictParam(luigi.Parameter):
    """
    Parameter for dictionaries.

    The parameter is serialised to a string using :py:mod:`pickle`.
    """

    def parse(self, string):
        """
        override of :py:meth:`luigi.Parameter.parse`.
        """

        sio = StringIO(string)
        return pickle.load(sio)

    def serialize(self, dct):
        """
        override of :py:meth:`luigi.Parameter.serialize`.
        """
        sio = BytesIO()
        pickle.dump(dct, sio)
        return sio.getvalue()


class RasterParam(DictParam):
    """
    Parameter representing a :py:class:`pyrate.shared.RasterBase` sub class.
    """

    def parse(self, string):
        """
        override of :py:meth:`DictParam.parse`.
        """
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
                'rasterBase must be an inscance DEM, '
                'Ifg or Incidence is valid')

    def serialize(self, rasterBase):
        """
        override of :py:meth:`DictParam.serialize`.
        """

        path = rasterBase.data_path

        if isinstance(rasterBase, DEM):
            d = {'type': 'DEM', 'path': path}
        elif isinstance(rasterBase, Ifg):
            d = {'type': 'Ifg', 'path': path}
        elif isinstance(rasterBase, Incidence):
            d = {'type': 'Incidence', 'path': path}
        else:
            raise luigi.parameter.UnknownParameterException(
                'rasterBase must be an inscance DEM, '
                'Ifg or Incidence is valid')

        return super(RasterParam, self).serialize(d)


def pythonify_config(config_file):
    """
    Make a copy of a PyRate config file which can be used by Python.

    Python and luigi config files *require* section headings. These are not
    present in the PyRate config file, so we make a copy of a PyRate config
    file and add a section heading to it so it can be parsed.
    """

    out_file = '{}.python'.format(config_file)

    with open(config_file, 'r') as inputFile:
        with open(out_file, 'w') as f:
            f.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            for line in inputFile:
                if any(x in line for x in [OBS_DIR, IFG_FILE_LIST,
                                           DEM_FILE, DEM_HEADER_FILE,
                                           OUT_DIR, ROIPAC_RESOURCE_HEADER]):
                    pos = line.find('~')
                    if pos != -1:
                        line = line[:pos] + os.environ['HOME'] + \
                               line[(pos+1):]    # create expanded line
                f.write(line)
    return out_file
