'''
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
'''

import os, re, pickle, luigi
from glob import glob
from os.path import join
from StringIO import StringIO
from pyrate.gamma import *



class GammaHasRun(luigi.task.ExternalTask):
    '''
    Phaux task used to ensure that the required outputs from GAMMA exist.
    '''

    fileName = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget(self.fileName)]



class ConvertToGeotiff(luigi.Task):
    '''
    Task responsible for converting a GAMMA file to GeoTif.
    '''

    inputFile = luigi.Parameter()
    outputDir = luigi.Parameter(default=None, significant=False)
    demHeaderPkl = luigi.Parameter(significant=False)
    noDataValue = luigi.Parameter(significant=False)

    def requires(self):
        '''
        Overload of :py:meth:`luigi.Task.requires`.

        Ensures that the required input exists.
        '''

        return [GammaHasRun(fileName = self.inputFile)]

    def output(self):
        '''
        Overload of :py:meth:`luigi.Task.output`.
        '''

        self.outputFile = '%s.tif' % os.path.splitext(os.path.basename(self.inputFile))[0]
        if self.outputDir is not None:
            self.outputFile = os.path.join(self.outputDir, self.outputFile)
        return [luigi.file.LocalTarget(self.outputFile)]

    def run(self):
        '''
        Overload of :py:meth:`luigi.Task.run`.
        '''

        sio = StringIO(self.demHeaderPkl)
        demHeader = pickle.load(sio)
        # find param files containing filename dates
        ptn = re.compile(r'\d{8}')  # match 8 digits for the datea
        dirName, fileName = os.path.split(self.inputFile)
        matches = ptn.findall(fileName)

        if len(matches) == 2:
            hpaths = [glob(join(dirName, '*%s*slc.par' % m))[0] for m in matches]
            tmp = [parse_epoch_header(hp) for hp in hpaths]
            hdrs = combine_headers(tmp[0], tmp[1], demHeader)
        else:
            # probably have DEM or incidence file
            hdrs = demHeader

        to_geotiff(hdrs, self.inputFile, self.outputFile, self.noDataValue)

