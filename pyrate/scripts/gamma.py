'''
Script to convert Gamma headers to ESRI's BIL format.

.. codeauthor:: Ben Davies, NCI
'''

import pickle, luigi
from StringIO import StringIO
from optparse import OptionParser
from pyrate.gamma import parse_dem_header
from pyrate.tasks.gamma import ConvertToGeotiff



def main():
    usage = 'Usage: %prog [options] DEM-HEADER GAMMA_FILE [GAMMA_FILE...]'
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--nodata', help='NODATA value', type='float', default=0.0)
    parser.add_option('-d', '--dest-dir', help='Write to DIR', type='string')
    options, args = parser.parse_args()

    if len(args) < 2:
        parser.error(usage)

    demHeaderFile = args[0]
    noData = options.nodata

    demHeader = parse_dem_header(demHeaderFile)
    sio = StringIO()
    pickle.dump(demHeader, sio)
    demHeaderPkl = sio.getvalue()

    tasks = [ConvertToGeotiff(
        inputFile=fn,
        outputDir=options.dest_dir,
        demHeaderPkl=demHeaderPkl,
        noDataValue=options.nodata) for fn in args[1:]]
    luigi.build(tasks, local_scheduler=True)



if __name__ == "__main__":
    main()

