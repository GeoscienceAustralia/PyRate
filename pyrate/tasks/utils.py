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
This Python module contains utilities for prepifg Luigi tasks.
"""
import os
import pickle
from io import StringIO, BytesIO
import luigi

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
    """
    Convenience class for specifying the parameter in a PyRate configuration
    file, intended for use in specifying which value in a PyRate configuration
    file should be used for a :py:class:`luigi.Parameter` (which is specified
    by the parameter *conf_path* to its constructor.)
    """

    def __init__(self, name):
        self['section'] = DUMMY_SECTION_NAME
        self['name'] = name
        super(InputParam, self).__init__()


class IfgListMixin(object):
    """
    Mixin to aid access to commonly used computed values from the PyRate
    config file.

    .. todo:: This should perhaps be renamed to something like *ConfigMixin*
        for clarity, as it is ued for accessing more than the list of
        interferograms.
    """

    ifg_list_file = luigi.Parameter(config_path=InputParam(
        config.IFG_FILE_LIST))
    obs_dir = luigi.Parameter(config_path=InputParam(config.OBS_DIR))
    out_dir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))

    def ifg_list(self, tif=True):
        """
        Get the list of interferograms to process.

        :param tif: Should the tif files be returned (*True*) or the raw
            interferograms (*False*). The latter will probably only be required
            before conversion to geotif files occurs.
        """

        file_names = config.parse_namelist(self.ifg_list_file)

        if tif:
            file_names = ['%s.tif' %
                          os.path.splitext(os.path.basename(fn))[0]
                          for fn in file_names]

        obs_dir = self.obs_dir
        if obs_dir:
            file_names = [os.path.join(obs_dir, fn) for fn in file_names]

        return file_names

    def ifg_tiff_list(self, tif=True):
        """
        Get the list of interferograms to process.

        :param tif: Should the tif files be returned (*True*) or the raw
            interferograms (*False*). The latter will probably only be required
            before conversion to geotif files occurs.
        """

        file_names = config.parse_namelist(self.ifg_list_file)

        if tif:
            file_names = [os.path.join(os.path.basename(fn).split('.')[0] +
                                       '_' + os.path.basename(fn).split('.')[1]
                                       + '.tif') for fn in file_names]

        out_dir = self.out_dir
        if out_dir:
            file_names = [os.path.join(out_dir, fn) for fn in file_names]

        return file_names

    @property
    def extents_file_name(self):
        """
        :return: extents file name
        """
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


class RasterParam(DictParam):  # pragma: no cover
    """
    Parameter representing a :py:class:`pyrate.shared.RasterBase` sub class.
    """

    def parse(self, string):
        """
        override of :py:meth:`DictParam.parse`.
        """
        dct = super(RasterParam, self).parse(string)
        raster_type = dct['type']
        path = dct['path']

        if raster_type == 'DEM':
            return DEM(path)
        elif raster_type == 'Ifg':
            return Ifg(path)
        elif raster_type == 'Incidence':
            return Incidence(path)
        else:  # pragma: no cover
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
        else:  # pragma: no cover
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

    with open(config_file, 'r') as input_file:
        with open(out_file, 'w') as f:
            f.write('[{}]\n'.format(DUMMY_SECTION_NAME))
            for line in input_file:
                if any(x in line for x in [OBS_DIR, IFG_FILE_LIST,
                                           DEM_FILE, DEM_HEADER_FILE,
                                           OUT_DIR, ROIPAC_RESOURCE_HEADER]):
                    pos = line.find('~')
                    if pos != -1:
                        line = line[:pos] + os.environ['HOME'] + \
                               line[(pos+1):]    # create expanded line
                f.write(line)
    return out_file
