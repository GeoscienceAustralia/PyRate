'''
* a script to show prepifg and pyrate output with matplotlib
* use as python inspect_output.py <prepifg or pyrate output file>
    * e.g. >>python inspect_output.py linrate.tif
* requires file be a TIFF and contain pr_output geotiff tag
* would be good to have this in path variable so could call at command line
'''

import sys
from PIL import Image, TiffImagePlugin
from osgeo import gdal
import numpy as np

import matplotlib.pyplot as pp
from matplotlib import figure

import pyrate.ifgconstants as ifc

if len(sys.argv[1:]) != 1:
    print 'USAGE: python inspect_output.py <prepifg or pyrate output file>'
    sys.exit(0)

# TODO check if it is a TIFF file

try:
    ds = gdal.Open(sys.argv[1])
except IOError:     # assuming this is the exception gdal raises (documentation is bad)
    print 'gdal probably could not open the file'
    sys.exit(0)

'''
md = ds.GetMetadata()
if not (ifc.PR_OUT_TYPE in md.keys()):
    print 'file is not a pyrate output file'
    sys.exit(0)
'''

try:
    tif = Image.open(fp=sys.argv[1])
except IOError as msg:
    print 'PIL could not open the file'
    print msg
    sys.exit(0)

tif_dat = np.array(tif.getdata()).reshape(tif.size[::-1])

'''
tif_dat_flat = np.array(tif.getdata())
for dat in tif_dat_flat:
    if np.isnan(dat):
        print 'got a nan'
print 'out of loop'
while True: pass
'''

plt = pp.imshow(tif_dat)
#pp.title('my plot')
#pp.title('my plot')
fig = pp.gcf()
fig.canvas.set_window_title(sys.argv[1])
pp.axis('off')
pp.show()