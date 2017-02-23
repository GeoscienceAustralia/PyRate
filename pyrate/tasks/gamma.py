"""
Luigi tasks to convert GAMMA headers to ESRI's BIL format.

GAMMA headers need to be translated into a GDAL recognisable format for use in
PyRate. This module translates GAMMA headers into into ESRI's BIL format,
allowing GDAL to access the raster data

The types of GAMMA files converted by PyRate are:

- DEM: with a .unw float32 binary data file (MSB order), & '.par' header. There
  is only a single height band in the binary data.

- Interferograms: have a .unw float32 bit binary data file (MSB order), with a
  single band for phase data. Two .par resource/header files, each containing
  details of the epochs used to create the interferogram. No geographic date is
  stored in these, so the DEM header is required for raster sizes/location etc.

The interferograms are geocoded/orthorectified to the DEM geometry, so all
datasets will share the same pixel size and dimensions.

.. todo:: describe incidence files (and any others (for later versions).
"""
# pylint: disable=attribute-defined-outside-init
import os
from os.path import join
import re
import glob2
import luigi
from pyrate import config
from pyrate.gamma import manage_headers
from pyrate.shared import write_geotiff
from pyrate.tasks.utils import IfgListMixin, InputParam

PTN = re.compile(r'\d{8}')  # match 8 digits for the dates


class GammaHasRun(luigi.task.ExternalTask):
    """
    Phaux task used to ensure that the required outputs from GAMMA exist.
    """

    file_name = luigi.Parameter()
    master_header = luigi.Parameter(default=None)
    slave_header = luigi.Parameter(default=None)

    def output(self):
        targets = [luigi.LocalTarget(self.file_name)]
        if self.master_header is not None:
            targets.append(luigi.LocalTarget(self.master_header))
        if self.slave_header is not None:
            targets.append(luigi.LocalTarget(self.slave_header))
        return targets


def get_header_paths(input_file, slc_dir=None):
    """
    function that matches input file names with header file names
    :param input_file: input gamma .unw file
    :return: corresponding header files that matches, or empty list if no match
    found
    """
    if slc_dir:
        dir_name = slc_dir
        _, file_name = os.path.split(input_file)
    else:  # header file must exist in the same dir as that of .unw
        dir_name, file_name = os.path.split(input_file)
    matches = PTN.findall(file_name)
    return [glob2.glob(join(dir_name, '**/*%s*slc.par' % m))[0]
            for m in matches]


class ConvertFileToGeotiff(luigi.Task):
    """
    Task responsible for converting a GAMMA file to GeoTif.
    """

    input_file = luigi.Parameter()
    demHeader_file = luigi.Parameter(
        config_path=InputParam(config.DEM_HEADER_FILE))
    out_dir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))
    no_data_value = luigi.FloatParameter(
        config_path=InputParam(config.NO_DATA_VALUE))
    slc_dir = luigi.Parameter(config_path=InputParam(config.SLC_DIR))

    def requires(self):
        """
        Overload of :py:meth:`luigi.Task.requires`.

        Ensures that the required input exists.
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

        self.out_file = os.path.join(
            self.out_dir,
            '%s.tif' % os.path.splitext(os.path.basename(self.input_file))[0])
        return [luigi.LocalTarget(self.out_file)]

    def run(self):
        """
        Overload of :py:meth:`luigi.Task.run`.
        """
        combined_header = manage_headers(self.demHeader_file,
                                         self.header_paths)
        write_geotiff(combined_header, self.input_file,
                      self.out_file, self.no_data_value)


class ConvertToGeotiff(IfgListMixin, luigi.WrapperTask):
    """ Wrapper class for gamma convert to geotiff"""
    def requires(self):
        return [ConvertFileToGeotiff(inputFile=fn)
                for fn in self.ifg_list(tif=False)]
