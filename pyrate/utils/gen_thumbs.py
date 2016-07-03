from PIL import Image
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from osgeo import gdal
from os.path import join
import math

import getopt
opts, args = getopt.getopt(sys.argv[1:], '2icv')        # 1 is one. not l
# usage checking
if opts == []:
    print 'read usage'
    sys.exit(0)
if not (len(args) == 1 or args ==[]):
    print 'read usage'
    sys.exit(0)

# gettting command line options

if args == []:
    direc = os.getcwd()
else:
    direc = args[0]

print direc

#while True: pass

do_ifg1 = False
do_ts_cum = False
do_ts_vel = False
do_ts_incr = False

for opt in opts:
    if opt[0] == '-2':
        # check if is a valid directory and ensure ends in '/' or '\'
        do_ifg1 = True
    if opt[0] == '-i':
        do_ts_incr = True
    if opt[0] == '-c':
        do_ts_cum = True
    if opt[0] == '-v':
        do_ts_vel = True
import pyrate.ifgconstants as ifc

# will get below from argument
#direc = r'C:\Users\gap\grk_viz\tif'
#direc = sys.argv[1]
#direc = '/home/gap/pr_testing/tests/orbitalfit/syd_g/1/python'

#print direc
#while True: pass

#lev = 0
for trav in os.walk(direc):
    #print trav
    break

#print trav
root = trav[0]
fns = trav[2]
#print root
#print fns

#while True: pass

ifg_1s, ifg_2s, ts_cums, ts_vels, ts_incrs = [], [], [], [], []
'''
# this would be the better way of doing it but i'm too lazy to write tests for these metadata additions
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
'''
# can't actually identify ifg_1s from file name. so can't do
for fn in fns:
    if fn[-6:] == 'cr.tif':
        ifg_2s.append(fn)
    elif fn[:6] == 'tscuml':
        ts_cums.append(fn)
    elif fn[:5] == 'tsvel':
        ts_vels.append(fn)
    elif fn[:6] == 'tsincr':
        ts_incrs.append(fn)

#print ifg_2s
#print ts_cums
#print ts_vels
#print ts_incrs

if do_ifg1:
    # create ifg_2 thumbnails...
    ifg_2s.sort()
    # have to assume they are of the same size...
    #img_dat = []        # will resuse this for each image type
    fig = plt.figure()
    fig.clf()

    # figure out how many rows we need
    img_ratio = Image.open(join(root, ifg_2s[0])).size      # gives (47, 72)
    #print img_ratio
    n_rows = int(math.ceil(float(len(ifg_2s))/5))
    #print n_rows
    #print len(ifg_2s)

    pix_per_row = float(192)*img_ratio[1]/float(img_ratio[0])           # 192 is 960/5... comes from 192/height = 47/72
    #print pix_per_row

    #while True: pass

    fig.set_size_inches(9.6, n_rows*pix_per_row/100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    for it, fn in enumerate(ifg_2s):
        img = Image.open(join(root, fn))
        img_dat = np.array(img.getdata()).reshape(img.size[::-1])
        plt.subplot(n_rows, 5, it+1)
        plt.imshow(img_dat)
        plt.axis('off')
        plt.tight_layout()
        plt.title(fn.split('_')[0], fontsize=10)

    #plt.show()
    #plt.set_size_inches(9.6, 9.6)

    plt.savefig(join(direc, 'thumbs_ifg2'), dpi=100)   # pixels per data point...
    plt.close()

if do_ts_vel:
    # create ifg_2 thumbnails...
    ts_vels.sort()
    # have to assume they are of the same size...
    #img_dat = []        # will resuse this for each image type
    fig = plt.figure()
    fig.clf()

    # figure out how many rows we need
    img_ratio = Image.open(join(root, ts_vels[0])).size      # gives (47, 72)
    #print img_ratio
    n_rows = int(math.ceil(float(len(ts_vels))/5))
    #print n_rows
    #print len(ts_vels)

    pix_per_row = float(192)*img_ratio[1]/float(img_ratio[0])           # 192 is 960/5... comes from 192/height = 47/72
    #print pix_per_row

    #while True: pass

    fig.set_size_inches(9.6, n_rows*pix_per_row/100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    for it, fn in enumerate(ts_vels):
        img = Image.open(join(root, fn))
        img_dat = np.array(img.getdata()).reshape(img.size[::-1])
        plt.subplot(n_rows, 5, it+1)
        plt.imshow(img_dat)
        plt.axis('off')
        plt.tight_layout()
        plt.title(fn, fontsize=10)

    #plt.show()
    #plt.set_size_inches(9.6, 9.6)

    plt.savefig(join(direc, 'thumbs_ts_vel'), dpi=100)   # pixels per data point...
    plt.close()

if do_ts_cum:
    # create ifg_2 thumbnails...
    ts_cums.sort()
    # have to assume they are of the same size...
    #img_dat = []        # will resuse this for each image type
    fig = plt.figure()
    fig.clf()

    # figure out how many rows we need
    img_ratio = Image.open(join(root, ts_cums[0])).size      # gives (47, 72)
    #print img_ratio
    n_rows = int(math.ceil(float(len(ts_cums))/5))
    #print n_rows
    #print len(ts_cums)

    pix_per_row = float(192)*img_ratio[1]/float(img_ratio[0])           # 192 is 960/5... comes from 192/height = 47/72
    #print pix_per_row

    #while True: pass

    fig.set_size_inches(9.6, n_rows*pix_per_row/100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    for it, fn in enumerate(ts_cums):
        img = Image.open(join(root, fn))
        img_dat = np.array(img.getdata()).reshape(img.size[::-1])
        plt.subplot(n_rows, 5, it+1)
        plt.imshow(img_dat)
        plt.axis('off')
        plt.tight_layout()
        plt.title(fn, fontsize=10)

    #plt.show()
    #plt.set_size_inches(9.6, 9.6)

    plt.savefig(join(direc, 'thumbs_ts_cuml'), dpi=100)   # pixels per data point...
    plt.close()

if do_ts_incr:
    # create ifg_2 thumbnails...
    ts_incrs.sort()
    # have to assume they are of the same size...
    #img_dat = []        # will resuse this for each image type
    fig = plt.figure()
    fig.clf()

    # figure out how many rows we need
    img_ratio = Image.open(join(root, ts_incrs[0])).size      # gives (47, 72)
    #print img_ratio
    n_rows = int(math.ceil(float(len(ts_incrs))/5))
    #print n_rows
    #print len(ts_cums)

    pix_per_row = float(192)*img_ratio[1]/float(img_ratio[0])           # 192 is 960/5... comes from 192/height = 47/72
    #print pix_per_row

    #while True: pass

    fig.set_size_inches(9.6, n_rows*pix_per_row/100)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    for it, fn in enumerate(ts_cums):
        img = Image.open(join(root, fn))
        img_dat = np.array(img.getdata()).reshape(img.size[::-1])
        plt.subplot(n_rows, 5, it+1)
        plt.imshow(img_dat)
        plt.axis('off')
        plt.tight_layout()
        plt.title(fn, fontsize=10)

    #plt.show()
    #plt.set_size_inches(9.6, 9.6)

    plt.savefig(join(direc, 'thumbs_ts_incr'), dpi=100)   # pixels per data point...
    plt.close()