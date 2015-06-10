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

import sys, luigi
from optparse import OptionParser
import pyrate.ifgconstants as ifc
from pyrate.tasks.roipac import ConvertToGeotiff



def main():
    usage = "Usage: %prog [options] ROIPAC_FILE [ROIPAC_FILE...]"
    parser = OptionParser(usage=usage)

    proj_help = 'GDAL well known projection (eg. "WGS84")'
    parser.add_option('-p', '--projection', help=proj_help, type='str')
    res_help = 'Resource/header file with projection data (usually DEM header)'
    parser.add_option('-r', '--resource-header', help=res_help, metavar='FILE')
    parser.add_option('-n', '--nodata', help='NODATA value', type='float', default=0.0)
    parser.add_option('-d', '--dest-dir', help='Write to DIR', type='string')
    options, args = parser.parse_args()

    if len(args) == 0:
        parser.error(usage)

    if options.resource_header:
        header = parse_header(options.resource_header)
        if ifc.PYRATE_DATUM not in header:
            if options.projection:
                projection = options.projection
            else:
                sys.exit('Error: header/resource file does not include DATUM and -p option not given')
        projection = header[ifc.PYRATE_DATUM]
    else:
        if options.projection:
            projection = options.projection
        else:
           sys.exit('Error: no header/resource file given and -p option not specified')

    # translate files
    tasks = [ConvertToGeotiff(
        inputFile = path,
        projection = projection,
        outputDir = options.dest_dir,
        noDataValue = options.nodata) for path in args]

    luigi.build(tasks, local_scheduler=True)



if __name__ == '__main__':
    main()

