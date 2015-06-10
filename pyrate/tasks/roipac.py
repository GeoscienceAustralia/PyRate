"""
Library/script to convert ROIPAC headers to ESRI's BIL format.

GDAL lacks a driver to parse ROIPAC headers. This module translates ROIPAC
headers into ESRI's BIL format, which is supported by GDAL. A basic command line
interface is provided for testing purposes.

The types of ROIPAC files/data used in PyRate are:

- *Interferograms*: a .unw 32 bit float data file, with a .rsc resource/header.
  The binary data is assumed to contain 2 bands, amplitude and phase.

- *DEM*: with a .unw 16 bit signed int binary data file, and a .rsc header
  There is only a single height band for the binary data.

.. todo:: Implement & describe incidence files, and any others (for later version)


There may be differences with the .rsc file content, with short and long forms.
The short form has 7 fields, covering raster size, location and wavelength. The
longer form can have up to 40 fields (see the test data for examples). PyRate
attempts to handle both forms of header.

Created on 12/09/2012

.. codeauthor:: Ben Davies, NCI <ben.davies@anu.edu.au>
"""

import os, luigi
import pyrate.ifgconstants as ifc
from pyrate.roipac import *



class RoipacHasRun(luigi.task.ExternalTask):
    '''
    Phaux task used to ensure that the required outputs from GAMMA exist.
    '''

    fileName = luigi.Parameter()
    headerFile = luigi.Parameter()

    def output(self):
        targets = [
            luigi.LocalTarget(self.fileName),
            luigi.LocalTarget(self.headerFile)]
        return targets



class ConvertToGeotiff(luigi.Task):
    '''
    Task responsible for converting a GAMMA file to GeoTif.
    '''

    inputFile = luigi.Parameter()
    projection = luigi.Parameter()
    outputDir = luigi.Parameter(default=None, significant=False)
    noDataValue = luigi.Parameter(significant=False)

    def requires(self):
        '''
        Overload of :py:meth:`luigi.Task.requires`.

        Ensures that the required input exists.
        '''

        self.headerFile = "%s.%s" % (self.inputFile, ROI_PAC_HEADER_FILE_EXT)
        tasks = [RoipacHasRun(
            fileName = self.inputFile,
            headerFile = self.headerFile)]

        return tasks

    def run(self):
        header = parse_header(self.headerFile)

        if ifc.PYRATE_DATUM not in header:  # DEM already has DATUM
            header[ifc.PYRATE_DATUM] = self.projection

        to_geotiff(header, self.inputFile, self.outputFile, self.noDataValue)

    def output(self):
        '''
        Overload of :py:meth:`luigi.Task.output`.

        .. todo:: This is the same as for gamma... refactor.
        '''

        self.outputFile = '%s.tif' % os.path.splitext(os.path.basename(self.inputFile))[0]
        if self.outputDir is not None:
            self.outputFile = os.path.join(self.outputDir, self.outputFile)
        return [luigi.file.LocalTarget(self.outputFile)]

