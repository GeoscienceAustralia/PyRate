"""
python utility to crop a larger ifg into smaller ifgs

example usage:
python pyrate/utils/crop_ifgs.py -i tests/sydney_test/tif/geo_060619-061002.tif
-o out.tif -e '150.91 -34.229999976 150.949166651  -34.17'

"""
from optparse import OptionParser
import subprocess


def crop_using_gdalwarp(input_file, output_file, extents):
    # TODO: add extents checking between input_file and extents
    extents_str = [str(e) for e in extents]
    cmd = ['gdalwarp', '-overwrite', '-srcnodata', 'None', '-q', '-te'] \
          + extents_str
    cmd += [input_file, output_file]
    subprocess.check_call(cmd)


if __name__ == '__main__':
    parser = OptionParser(usage='%prog -i input_file -o output_file'
                                ' -e extents\n'
                                'Crop a larger interferrogram into '
                                'smaller ones')
    parser.add_option('-i', '--input', type=str, dest='input_file',
                      help='name of input interferrogram')

    parser.add_option('-o', '--out', type=str, dest='output_file',
                      help='name of cropped output interferrogram')

    parser.add_option('-e', '--extents', type=str, dest='extents',
                      help='extents to be used for the cropped file.\n'
                           'needs to be a list of 4 floats with spaces\n'
                           'example: '
                           "-e '150.91 -34.229999976 150.949166651 -34.17'")

    options, args = parser.parse_args()

    if not options.input_file:  # if filename is not given
        parser.error('Input filename not given.')

    if not options.output_file:  # if filename is not given
        parser.error('Output filename not given.')

    if not options.extents:  # if filename is not given
        parser.error('Crop extents must be provided')

    extents = [float(t) for t in options.extents.split()]

    if len(extents) != 4:
        raise AttributeError('extents to be used for the cropped file.\n'
                             'needs to be a list or tuples of 4 floats\n'
                             "example:"
                             "--extents "
                             "'150.91 -34.229999976 150.949166651 -34.17'")

    crop_using_gdalwarp(input_file=options.input_file,
                        output_file=options.output_file,
                        extents=extents)
