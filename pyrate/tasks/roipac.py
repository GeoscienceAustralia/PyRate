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
"""
import os
import luigi
import pyrate.config as config
from pyrate.roipac import manage_header, ROI_PAC_HEADER_FILE_EXT, parse_header
from pyrate import ifgconstants as ifc
from pyrate.shared import write_geotiff, output_tiff_filename
from pyrate.tasks.utils import InputParam, IfgListMixin


class RoipacHasRun(luigi.task.ExternalTask):
    """
    Phaux task used to ensure that the required outputs from ROIPAC exist.
    """

    fileName = luigi.Parameter()
    headerFile = luigi.Parameter()

    def output(self):
        targets = [
            luigi.LocalTarget(self.fileName),
            luigi.LocalTarget(self.headerFile)]
        return targets


class ResourceHeaderExists(luigi.ExternalTask):
    """
    Ensure that the resource header exists.
    """

    resourceHeader = luigi.Parameter()
    property = 100

    def output(self):
        return [luigi.LocalTarget(self.resourceHeader)]


class ConvertFileToGeotiff(luigi.Task):
    """
    Task responsible for converting a ROIPAC file to GeoTif.
    """
    # pylint: disable=attribute-defined-outside-init
    inputFile = luigi.Parameter()
    projection = luigi.Parameter()
    outputDir = luigi.Parameter(config_path=InputParam(config.OUT_DIR))
    no_data_value = luigi.FloatParameter(
        config_path=InputParam(config.NO_DATA_VALUE))

    def requires(self):
        """
        Overload of :py:meth:`luigi.Task.requires`.

        Ensures that the required input exists.
        """
        self.header_file = "%s.%s" % (self.inputFile, ROI_PAC_HEADER_FILE_EXT)
        tasks = [RoipacHasRun(
            fileName=self.inputFile,
            headerFile=self.header_file)]

        return tasks

    def run(self):
        header = manage_header(self.header_file, self.projection)
        write_geotiff(header, self.inputFile, self.output_file,
                      self.no_data_value)

    def output(self):
        """
        Overload of :py:meth:`luigi.Task.output`.
        """
        self.output_file = output_tiff_filename(self.inputFile, self.outputDir)
        return [luigi.LocalTarget(self.output_file)]


class _DoConvertToGeotiffRoipac(IfgListMixin, luigi.WrapperTask):
    projection = luigi.Parameter(
        default=None,
        config_path=InputParam(config.INPUT_IFG_PROJECTION))

    resourceHeader = luigi.Parameter(
        default=None,
        config_path=InputParam(config.ROIPAC_RESOURCE_HEADER))

    @property
    def priority(self):
        """
        The requires method of this Task *may* reqire the existence of a header
        file... so that needs to be checked first.
        """
        return ResourceHeaderExists.priority - 1

    def requires(self):
        if self.resourceHeader is not None:
            header = parse_header(self.resourceHeader)
            if ifc.PYRATE_DATUM not in header:
                raise Exception('Error: header/resource file does not '
                                'include DATUM and -p option not given')
            projection = header[ifc.PYRATE_DATUM]
        else:
            if self.projection:
                projection = self.projection
            else:
                raise Exception('Error: no header/resource file given '
                                'and -p option not specified')

        ifg_files = self.ifg_list(tif=False)
        tasks = [ConvertFileToGeotiff(
            inputFile=path,
            projection=projection) for path in ifg_files]
        return tasks


class ConvertToGeotiff(luigi.WrapperTask):
    """
    Convert ROIPAC files to geotifs.

    This delegates the actual conversions tasks which operate on individual
    files and is purely a convenience wrapper.
    """

    resourceHeader = luigi.Parameter(
        default=None,
        config_path=InputParam(config.ROIPAC_RESOURCE_HEADER))

    def requires(self):
        tasks = [_DoConvertToGeotiffRoipac()]
        if self.resourceHeader is not None:
            tasks.append(ResourceHeaderExists(self.resourceHeader))

        return tasks
