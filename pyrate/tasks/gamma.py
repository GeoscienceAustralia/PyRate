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
This Python module is a Luigi wrapper for converting GAMMA format input data.
"""
# pylint: disable=attribute-defined-outside-init
import luigi
from pyrate import config
from pyrate.gamma import manage_headers, get_header_paths
from pyrate.shared import write_geotiff, output_tiff_filename
from pyrate.tasks.utils import IfgListMixin, InputParam


class GammaHasRun(luigi.task.ExternalTask):
    """
    Phantom task used to ensure that the required outputs from GAMMA exist.
    """
    fileName = luigi.Parameter()
    masterHeader = luigi.Parameter(default=None)
    slaveHeader = luigi.Parameter(default=None)

    def output(self):
        targets = [luigi.LocalTarget(self.fileName)]
        if self.masterHeader is not None:
            targets.append(luigi.LocalTarget(self.masterHeader))
        if self.slaveHeader is not None:
            targets.append(luigi.LocalTarget(self.slaveHeader))
        return targets


class ConvertFileToGeotiff(luigi.Task):
    """
    Task responsible for converting a GAMMA file to GeoTiff.
    """
    input_file = luigi.Parameter()
    demHeaderFile = luigi.Parameter(
        config_path=InputParam(config.DEM_HEADER_FILE))
    out_dir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))
    no_data_value = luigi.FloatParameter(
        config_path=InputParam(config.NO_DATA_VALUE))
    slc_dir = luigi.Parameter(config_path=InputParam(config.SLC_DIR))

    def requires(self):
        """
        Overload of :py:meth:`luigi.Task.requires`. Ensures that the required
        input exists.
        """
        self.header_paths = get_header_paths(self.input_file, self.slc_dir)

        if len(self.header_paths) == 2:
            tasks = [GammaHasRun(
                fileName=self.input_file,
                masterHeader=self.header_paths[0],
                slaveHeader=self.header_paths[1])]
        else:
            tasks = [GammaHasRun(fileName=self.input_file)]

        return tasks

    def output(self):
        """
        Overload of :py:meth:`luigi.Task.output`.
        """
        self.out_file = output_tiff_filename(self.input_file, self.out_dir)
        return [luigi.LocalTarget(self.out_file)]

    def run(self):
        """
        Overload of :py:meth:`luigi.Task.run`.
        """
        combined_header = manage_headers(self.demHeaderFile, self.header_paths)
        write_geotiff(combined_header, self.input_file,
                      self.out_file, self.no_data_value)


class ConvertToGeotiff(IfgListMixin, luigi.WrapperTask):
    """ Wrapper class for GAMMA conversion to geotiff"""
    def requires(self):
        return [ConvertFileToGeotiff(input_file=fn)
                for fn in self.ifg_list(tif=False)]
