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
This Python module is a Luigi wrapper for the prepifg module.
"""
import os
import pickle
import luigi

import pyrate.config as cf
from pyrate.prepifg import (
    Ifg,
    get_analysis_extent,
    prepare_ifg,
    PreprocessError)
from pyrate.tasks.converttogeotif import ConvertToGeotiff
from pyrate.tasks.utils import (
    IfgListMixin,
    InputParam,
    RasterParam)
from pyrate.shared import warp_required


class GetAnalysisExtents(IfgListMixin, luigi.Task):
    """
    Dummy analysis extents gather class used during luigi tasks
    """
    crop_opt = luigi.IntParameter(config_path=InputParam(cf.IFG_CROP_OPT))
    ifgx_first = luigi.FloatParameter(default=None,
                                      config_path=InputParam(cf.IFG_XFIRST))
    ifgy_first = luigi.FloatParameter(default=None,
                                      config_path=InputParam(cf.IFG_YFIRST))
    ifgx_last = luigi.FloatParameter(default=None,
                                     config_path=InputParam(cf.IFG_XLAST))
    ifgy_last = luigi.FloatParameter(default=None,
                                     config_path=InputParam(cf.IFG_YLAST))
    xlooks = luigi.IntParameter(config_path=InputParam(cf.IFG_LKSX))
    ylooks = luigi.IntParameter(config_path=InputParam(cf.IFG_LKSY))

    def requires(self):
        return [ConvertToGeotiff()]

    def run(self):
        user_exts = (self.ifgx_first, self.ifgy_first,
                     self.ifgx_last, self.ifgy_last)

        if not all(user_exts):
            if self.crop_opt == 3:
                raise PreprocessError('No custom cropping extents specified')
            user_exts = None

        ifgs = [Ifg(path) for path in self.ifg_tiff_list()]

        extents = get_analysis_extent(
            self.crop_opt,
            ifgs,
            self.xlooks,
            self.ylooks,
            user_exts)

        with open(self.extents_file_name, 'wb') as ext_file:
            pickle.dump(extents, ext_file)

    def output(self):
        return luigi.LocalTarget(self.extents_file_name)


class PrepareInterferogram(IfgListMixin, luigi.WrapperTask):
    """
    Produces multilooked/resampled data files for PyRate analysis.

    :param ifgs: sequence of Ifg objs (DEM obj may be included for processing)
    :param crop_opt: integer cropping type option (see config)
    :param xlooks: multilooking factor for the X axis
    :param ylooks: Y axis multilooking factor
    :param float thresh: (0.0, 1.0). Controls NaN handling when resampling to
    coarser grids. Value is the proportion above which the number of NaNs in
    an area is considered invalid. thresh=0 resamples to NaN if 1 or more
    contributing cells are NaNs. At 0.25, it resamples to NaN if 1/4 or
    more contributing cells are NaNs. At 1.0, areas are resampled to NaN
    only if all contributing cells are NaNs.
    :param user_exts: CustomExts tuple with user sepcified lat long corners
    :param verbose: Controls level of gdalwarp output
    """
    # pylint: disable=bad-super-call, no-member

    ifg = RasterParam()
    thresh = luigi.FloatParameter(config_path=InputParam(
        cf.NO_DATA_AVERAGING_THRESHOLD))
    crop_opt = luigi.IntParameter(config_path=InputParam(cf.IFG_CROP_OPT))
    xlooks = luigi.IntParameter(config_path=InputParam(cf.IFG_LKSX))
    ylooks = luigi.IntParameter(config_path=InputParam(cf.IFG_LKSY))
    # verbose = luigi.BooleanParameter(default=True, significant=False)

    def requires(self):
        return [GetAnalysisExtents()]

    def run(self):
        with open(self.extents_file_name, 'rb') as ext_file:
            extents = pickle.load(ext_file)
        prepare_ifg(
            self.ifg.data_path,
            self.xlooks,
            self.ylooks,
            extents,
            self.thresh,
            self.crop_opt
            )
        self.ifg.close()

    def output(self):
        if warp_required(self.xlooks, self.ylooks, self.crop_opt):
            return luigi.LocalTarget(
                cf.mlooked_path(self.ifg.data_path,
                                self.ylooks,
                                self.crop_opt))
        else:
            return []

    def complete(self):
        if self.output():
            return super(luigi.WrapperTask, self).complete()
        else:
            # then this didn't do anything, so check that
            # the requres are complete
            # ... which is exactly what luigi.WrapperTask does.
            # TODO: This requires knowledge of prepare_ifg,
            # that is opaque. Address with refactoring.
            return super(PrepareInterferogram, self).complete()


class PrepareInterferograms(IfgListMixin, luigi.WrapperTask):
    """ Luigi wrapper class """

    def __init__(self, *args, **kwargs):
        super(PrepareInterferograms, self).__init__(*args, **kwargs)
        self.extents_removed = False

    def requires(self):
        return [PrepareInterferogram(ifg=Ifg(path))
                for path in self.ifg_tiff_list()]

    def run(self):
        try:
            if os.path.exists(self.extents_file_name):
                os.remove(self.extents_file_name)
        except:
            raise PrepifgException(
                'Extents file was not found in the desired '
                'location: {}'.format(self.extents_file_name),
                'Make sure your paths are setup correctly in config file')

        self.extents_removed = True

    def complete(self):
        return self.extents_removed and \
               super(PrepareInterferograms, self).complete()


class PrepifgException(Exception):
    """
    Prepifg exception class
    """
