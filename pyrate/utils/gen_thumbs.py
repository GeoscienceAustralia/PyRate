from PIL import Image
import os
import numpy as np
import matplotlib.pyplot as pp
import sys
from osgeo import gdal

import pyrate.ifgconstants as ifc

direc = r'C:\Users\gap\grk_viz\tif'
direc = sys.argv[1]

#print direc
#while True: pass

lev = 0
for trav in os.walk(direc):
    #print trav
    break

#print trav
root = trav[0]
fns = trav[2]
print root
print fns
#while True: pass

ifg_1s, ts_cums, ts_vels, ts_incrs = [], [], [], []
for fn in fns:
    ds = gdal.Open(os.path.join(root,fn))
    md = ds.GetMetadata()
    if md[ifc.PR_OUT_TYPE] == 'ifg_1':
        ifg_1s.append(fn)   # hopefully this just creates a reference and not a copy. don't want string copies
    if md[ifc.PR_OUT_TYPE] == 'ts_cum':
        pass
    if md[ifc.PR_OUT_TYPE] == 'ts_vel':
        pass
    if md[ifc.PR_OUT_TYPE] == 'ts_incr':
        pass


print tifs

# get actual data
tif_dats = [np.array(tif.getdata()).reshape(tif.size[::-1]) for tif in tifs]    # reverse shape. works. don't know why

#print tif_dats

# matplotlibification
fig = pp.figure()
for it, tif_dat in enumerate(tif_dats):
    plt = fig.add_subplot(4, 5, it+1)
    pp.imshow(tif_dat)


#print dir(plt)

pp.show()



#img = Image.open()