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

.. codeauthor:: Ben Davies, NCI
"""

import os, re, luigi
import glob2
from os.path import join
from pyrate import config
from pyrate.gamma import *
from pyrate.tasks.utils import IfgListMixin, InputParam

PTN = re.compile(r'\d{8}')  # match 8 digits for the dates

class GammaHasRun(luigi.task.ExternalTask):
    """
    Phaux task used to ensure that the required outputs from GAMMA exist.
    """

    fileName     = luigi.Parameter()
    masterHeader = luigi.Parameter(default=None)
    slaveHeader  = luigi.Parameter(default=None)

    def output(self):
        targets = [luigi.LocalTarget(self.fileName)]
        if self.masterHeader is not None:
            targets.append(luigi.LocalTarget(self.masterHeader))
        if self.slaveHeader is not None:
            targets.append(luigi.LocalTarget(self.slaveHeader))
        return targets


def get_header_paths(input_file, slc_dir=None):
    """
    function that matches input file names with header file names
    :param input_file: input gamma .unw file
    :return: corresponding header files that matches, or empty list if no match
    found
    """
    if slc_dir:
        dirName = slc_dir
        _, fileName = os.path.split(input_file)
    else:  # header file must exist in the same dir as that of .unw
        dirName, fileName = os.path.split(input_file)
    matches = PTN.findall(fileName)
    return [glob2.glob(join(dirName, '*%s*slc.par' % m))[0] for m in matches]


class ConvertFileToGeotiff(luigi.Task):
    """
    Task responsible for converting a GAMMA file to GeoTif.
    """

    inputFile = luigi.Parameter()
    demHeaderFile = luigi.Parameter(
        config_path=InputParam(config.DEM_HEADER_FILE))
    outputDir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))
    noDataValue = luigi.FloatParameter(
        config_path=InputParam(config.NO_DATA_VALUE))
    slc_dir = luigi.Parameter(config_path=InputParam(config.SLC_DIR))

    def requires(self):
        """
        Overload of :py:meth:`luigi.Task.requires`.

        Ensures that the required input exists.
        """
        self.headerPaths = get_header_paths(self.inputFile, self.slc_dir)

        if len(self.headerPaths) == 2:
            tasks = [GammaHasRun(
                fileName=self.inputFile,
                masterHeader=self.headerPaths[0],
                slaveHeader=self.headerPaths[1])]
        else:
            tasks = [GammaHasRun(fileName=self.inputFile)]

        return tasks

    def output(self):
        """
        Overload of :py:meth:`luigi.Task.output`.
        """

        self.outputFile = os.path.join(
            self.outputDir,
            '%s.tif' % os.path.splitext(os.path.basename(self.inputFile))[0])
        return [luigi.file.LocalTarget(self.outputFile)]

    def run(self):
        """
        Overload of :py:meth:`luigi.Task.run`.
        """
        demHeader = parse_dem_header(self.demHeaderFile)

        # find param files containing filename dates
        if len(self.headerPaths) == 2:
            headers = [parse_epoch_header(hp) for hp in self.headerPaths]
            combinedHeader = combine_headers(headers[0], headers[1], demHeader)
        else:
            # probably have DEM or incidence file
            combinedHeader = demHeader
        to_geotiff(combinedHeader, self.inputFile,
                   self.outputFile, self.noDataValue)


class ConvertToGeotiff(IfgListMixin, luigi.WrapperTask):
    def requires(self):
        return [ConvertFileToGeotiff(inputFile=fn)
                for fn in self.ifgList(tif=False)]
