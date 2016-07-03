'''
v1
==
* compare the outputs of run_cfgs.py

usage
-----
* $ python cmp_out.py
* e.g. $ python pyrate/utils/cmp_out.py --od ~/pr_testing/output/luigi_networkx -d ~/pr_testing/tests/luigi_networkx/syd_g/ -htl
'''

import sys
import os
from os.path import join
import time
import numpy as np
import pickle
import scipy.io     # scipy is the package. need must always import the module (io)
from decimal import Decimal

from PIL import Image
from PIL import ImageDraw

import matplotlib.pyplot as plt
from matplotlib import figure

MPL_SCALING = 11        # matplotlib scaling

# =================================================================
# MARKER DEFINITIONS...
# =================================================================
marker_tol = Image.new('RGBA', (MPL_SCALING, MPL_SCALING), (0, 0, 0, 0))      # fully transparent
marker_nan1 = Image.new('RGBA', (MPL_SCALING, MPL_SCALING), (0, 0, 0, 0))
marker_nan2 = Image.new('RGBA', (MPL_SCALING, MPL_SCALING), (0, 0, 0, 0))

# blue square box for tolerance errors
drawer = ImageDraw.ImageDraw(marker_tol)
drawer.line((0, 0, 0, MPL_SCALING-1), (0, 0, 255, 255))     # fully opaque
drawer.line((0, 0, MPL_SCALING-1, 0), (0, 0, 255, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, MPL_SCALING-1, 0), (0, 0, 255, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, 0, MPL_SCALING-1), (0, 0, 255, 255))

# green '+' for nan type 1
drawer = ImageDraw.ImageDraw(marker_nan1)
drawer.line((0, 0, 0, MPL_SCALING-1), (0, 255, 0, 255))     # fully opaque
drawer.line((0, 0, MPL_SCALING-1, 0), (0, 255, 0, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, MPL_SCALING-1, 0), (0, 255, 0, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, 0, MPL_SCALING-1), (0, 255, 0, 255))
# ------------------
mid = MPL_SCALING/2 # don't need to +1 because indexing starts @ 0
drawer.line((0, mid, MPL_SCALING-1, mid), (0, 255, 0, 255))
drawer.line((mid, 0, mid, MPL_SCALING-1), (0, 255, 0, 255))

# green 'x' for nan type 2
drawer = ImageDraw.ImageDraw(marker_nan2)
drawer.line((0, 0, 0, MPL_SCALING-1), (0, 255, 0, 255))     # fully opaque
drawer.line((0, 0, MPL_SCALING-1, 0), (0, 255, 0, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, MPL_SCALING-1, 0), (0, 255, 0, 255))
drawer.line((MPL_SCALING-1, MPL_SCALING-1, 0, MPL_SCALING-1), (0, 255, 0, 255))
# ------------------
drawer.line((0, 0, MPL_SCALING-1, MPL_SCALING-1), (0, 255, 0, 255))
drawer.line((0, MPL_SCALING-1, MPL_SCALING-1, 0), (0, 255, 0, 255))

# =================================================================
# ARRAY CHECKING FUNCTIONS...
# these could really be made into just 1 function
# =================================================================
# create whether actually do it or not. change later
fig = plt.figure()
#fig2 = plt.figure()

# histogram data (percentages)
# get the length of these to get the total amount of errors
hist_rate_tol = []
hist_rate_nan1 = []
hist_rate_nan2 = []
hist_rate_nans = []
hist_rate_tot = []
# ------------------------------
hist_ts_tol = []
hist_ts_nan1 = []
hist_ts_nan2 = []
hist_ts_nans = []
hist_ts_tot = []

n_rate_ref = 0
n_ts_ref = 0
n_no_ref_fail = 0   # could probably get this some better way. but watever...

n_rate_total = 0
n_ts_total = 0

M1_N = 'py'  # matrix name. should be passed as a function
M2_N = 'mt'
def chk_out_mats(m1, m2):
    global fig, fig2 # these might not need to be global
    global n_rate_ref, n_ts_ref, n_no_ref_fail
    global n_rate_total, n_ts_total
    # appending to histograms shouldn't need require them being global?
    # checks pirate/pyrate (m1/m2) output arrays
    # will determine if 2d or 3d... can only deal with either of those two
    # returns an error string that you just plonk into the error file
    n_rate_tol, n_rate_nan1, n_rate_nan2 = (0, 0, 0)
    n_ts_tol, n_ts_nan1, n_ts_nan2 = (0, 0, 0)

    fun_mode = len(m1.shape)              # function mode. 2d or 3d matrices. don't bother checking otherwise. assume can't happen

    relative_tolerance=0.1      # todo: have this as function and program argument

    er_str = ''     # the error string we build up

    # check reference points here. not the right place but whatever
    global rp_mt, rp_py
    if rp_mt != rp_py:
        er_str += '\t'*4+'* REFERENCE PIXELS DO NOT MATCH...\n'
        er_str += '\t'*5+'* '+M1_N+' is '+repr(rp_py)+'\n'
        er_str += '\t'*5+'* '+M2_N+' is '+repr(rp_mt)+'\n'
        # don't include in histogram or others... only keep track of in n_ref_fail
        if fun_mode == 2: n_rate_ref += 1
        if fun_mode == 3: n_ts_ref += 1
        return er_str           # exit here. won't be anything meaningful

    n_no_ref_fail += 1

    er_array_shapes = ''
    if m1.shape != m2.shape:
        if verbosity == 1:
            er_array_shapes += '\t'*4+'* ARRAY SHAPES DO NOT MATCH. CHECKING WITH STRIPPED DATA...\n'
            er_array_shapes += '\t'*5+'* '+M1_N+' is '+repr(m1.shape)+'\n'
            er_array_shapes += '\t'*5+'* '+M2_N+' is '+repr(m2.shape)+'\n'

        # strip rows or cols to make have the same shape
        # don't worry about epochs. you'd hope they're equal...
        shape_diff = tuple(np.subtract(m1.shape, m2.shape))
        for elem in [0, 1]:
            if shape_diff[elem] > 0:
                # strip diff off m1
                if elem == 0:
                    if len(shape_diff) == 2:
                        m1 = m1[:m1.shape[elem]-shape_diff[elem],:]
                    else:
                        # must be 3d array
                        m1 = m1[:m1.shape[elem]-shape_diff[elem],:,:]
                else:
                    # must be elem == 1
                    if len(shape_diff) == 2:
                        m1 = m1[:,:m1.shape[elem]-shape_diff[elem]]
                    else:
                        m1 = m1[:,:m1.shape[elem]-shape_diff[elem],:]
            else:
                # must be negative. strip abs off m2
                if elem == 0:
                    if len(shape_diff) == 2:
                        m2 = m2[:m2.shape[elem]-abs(shape_diff[elem]),:]
                    else:
                        m2 = m2[:m2.shape[elem]-abs(shape_diff[elem]),:,:]
                else:
                    # must be elem == 1
                    if len(shape_diff) == 2:
                        m2 = m2[:,:m2.shape[elem]-abs(shape_diff[elem])]
                    else:
                        m2 = m2[:,:m2.shape[elem]-abs(shape_diff[elem]),:]

    # sanity checking
    assert m1.size == m2.size, 'array shapes do not match up even after stripping'
    '''
    print len(m1.shape)
    while True: pass
    '''

    if fun_mode == 2:
        n_r, n_c = m1.shape
        n_e = 1
    else:   # gauranteed to be 3. else is faster than testing for 3
        n_r, n_c, n_e = m1.shape        # number of rows, cols, epochs

    it_r, it_c, it_e = (0,0,0)      # iterators. not pythonic, but whatever

    # open up an image to put pixels to
    #img = Image.new('RGB', (n_r, n_c), (255, 255, 255))
    if fun_mode == 2 and do_lin:
        #fig = plt.figure()     # this creates a new figure. use below instead.
        fig.clf()
        fig.set_size_inches(m2.shape[1], m2.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        plt.set_cmap('hot')
        ax.imshow(m2, aspect = 'normal')
        img_fn = join(image_direct, 'rate_error_'+cont_fld+'.png')
        plt.savefig(img_fn, dpi=MPL_SCALING)   # pixels per data point...
        # open it up in PIL now so we can add our hit markers *pttt* *pttt* <-- mw2. MLG!
        img = Image.open(img_fn)

    while it_e < n_e:
        if fun_mode == 3 and do_ts:
            #fig = plt.figure()     # this creates a new figure. use below instead.
            fig.clf()
            fig.set_size_inches(m2.shape[1], m2.shape[0])
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            plt.set_cmap('hot')
            ax.imshow(m2[:,:,it_e], aspect = 'normal')      # background as matlab epoch
            this_d = join(image_direct2, cont_fld)
            if not os.path.isdir(this_d):
                # create the directory
                os.mkdir(this_d)
            img_fn = join(this_d, 'ts_error_epoch'+str(it_e)+'.png')
            plt.savefig(img_fn, dpi=MPL_SCALING)   # pixels per data point...
            # open it up in PIL now so we can add our hit markers *pttt* *pttt* <-- mw2. MLG!
            img = Image.open(img_fn)

        it_r = 0
        while it_r < n_r:
            it_c = 0
            while it_c < n_c:
                if fun_mode == 2:
                    it = (it_r, it_c)
                else:
                    it = (it_r, it_c, it_e)     # this can be used to index ndarray
                if (not np.isnan(m1[it])) and (not np.isnan(m2[it])):  # both floats
                    # Duncan's way of checking if floats are equal...
                    if abs(m1[it].item() - m2[it].item()) > relative_tolerance:
                        if fun_mode == 2:
                            n_rate_tol += 1
                            if do_lin: img.paste(marker_tol, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_tol)
                        if fun_mode == 3:
                            n_ts_tol += 1
                            if do_ts: img.paste(marker_tol, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_tol)
                        if verbosity == 1 or verbosity == 2:
                            er_str += '\t'*4+'* do not match to tolerance @ '+str(it)+'\n'
                            er_str += '\t'*5+'* '+M1_N+' = '+str(Decimal(m1[it].item()))+'\n'
                            er_str += '\t'*5+'* '+M2_N+' = '+str(Decimal(m2[it].item()))+'\n'
                else:   # atleast one of py[it], m2[it] is NaN
                    if (not np.isnan(m1[it])) and (np.isnan(m2[it])):
                        if fun_mode == 2:
                            n_rate_nan1 += 1
                            if do_lin: img.paste(marker_nan1, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_nan1)
                        if fun_mode == 3:
                            n_ts_nan1 += 1
                            if do_ts: img.paste(marker_nan1, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_nan1)
                        if verbosity == 1 or verbosity == 2:
                            er_str += '\t'*4+'* '+M1_N+' should be NaN (but is not) @ '+str(it)+'\n'
                    if (np.isnan(m1[it])) and (not np.isnan(m2[it])):
                        if fun_mode == 2:
                            n_rate_nan2 += 1
                            if do_lin: img.paste(marker_nan2, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_nan2)
                        if fun_mode == 3:
                            n_ts_nan2 += 1
                            if do_ts: img.paste(marker_nan2, (it_c*MPL_SCALING, it_r*MPL_SCALING), marker_nan2)
                        if verbosity == 1 or verbosity == 2:
                            er_str += '\t'*4+'* '+M1_N+' should not be NaN (but is) @ '+str(it)+'\n'
                it_c += 1
            it_r += 1
        # duplication not needed. fix later
        if fun_mode == 2 and do_lin:
            img.save(img_fn, 'PNG')
        if fun_mode == 3 and do_ts:
            img.save(img_fn, 'PNG')
        it_e += 1

    er_str_prep = ''
    '''
    print n_tol_mis
    print n_nan_mis
    while True: pass
    '''

    tot_pix = m1.shape[0]*m1.shape[1]*n_e       # will be the same for either function mode
    if fun_mode == 2:
        n_rate_total += 1   # add to number of tests...
        percent_tol = (float(n_rate_tol)/tot_pix)*100
        if percent_tol != 0.0: hist_rate_tol.append(percent_tol)       # append only if not zero...
        percent_nan1 = (float(n_rate_nan1)/tot_pix)*100
        if percent_nan1 != 0.0: hist_rate_nan1.append(percent_nan1)
        percent_nan2 = (float(n_rate_nan2)/tot_pix)*100
        if percent_nan2 != 0.0: hist_rate_nan2.append(percent_nan2)
        percent_nans = percent_nan1 + percent_nan2
        if percent_nans != 0.0: hist_rate_nans.append(percent_nans)
        percent_tot = percent_nans + percent_tol
        if percent_tot != 0.0: hist_rate_tot.append(percent_tot)

    if fun_mode == 3:
        #print n_ts_tol
        #print tot_pix
        n_ts_total += 1
        percent_tol = (float(n_ts_tol)/tot_pix)*100
        if percent_tol != 0.0: hist_ts_tol.append(percent_tol)
        percent_nan1 = (float(n_ts_nan1)/tot_pix)*100
        if percent_nan1 != 0.0: hist_ts_nan1.append(percent_nan1)
        percent_nan2 = (float(n_ts_nan2)/tot_pix)*100
        if percent_nan2 != 0.0: hist_ts_nan2.append(percent_nan2)
        percent_nans = percent_nan1 + percent_nan2
        if percent_nans != 0.0: hist_ts_nans.append(percent_nans)
        percent_tot = percent_nans + percent_tol
        if percent_tot != 0.0: hist_ts_tot.append(percent_tot)

    if percent_tot != 0.0:        # only print to log file if have errors...
        er_str_prep += '\t'*4+'* tolerance errors = '+str(percent_tol)+'%\n'
        er_str_prep += '\t'*4+'* NaN errors       = '+str(percent_nans)+'%\n'
        er_str_prep += '\t'*4+'* total fail       = '+str(percent_tot)+'%\n'

    er_str = er_array_shapes+er_str_prep+er_str
    return er_str
# =================================================================
# =================================================================

import getopt
opts, args = getopt.getopt(sys.argv[1:], 'd:i:j:v:htl', ['od=','of='])
# usage checking
if opts == [] and args != []:
    print 'read usage'
    sys.exit(0)
# gettting command line options
root_direct = False
ign_list = []
jst_list = []
out_direct = False
out_fn = False
verbosity = 3       # default verbosity (don't tell about stripping rows)
# flags for creating linear rate & ts error images, histograms
do_lin = False
do_ts = False
do_hist = False
# assiging and checking valid options
for opt in opts:
    if opt[0] == '-d':
        # check if is a valid directory and ensure ends in '/' or '\'
        root_direct = opt[1]
    if opt[0] == '-i':
        ign_list = opt[1]
    if opt[0] == '-j':
        jst_list = opt[1]
    if opt[0] == '--od':
        out_direct = opt[1]
    if opt[0] == '--of':
        out_fn = opt[1]
    if opt[0] == '-v':
        verbosity = int(opt[1])
    if opt[0] == '-h':
        do_hist = True
    if opt[0] == '-t':
        do_ts = True
    if opt[0] == '-l':
        do_lin = True
# checking for valid output options
'''
print opts
print do_lin, do_ts, do_hist
while True: pass
'''
if out_direct != False and out_fn != False:
    print 'can only have one or the other or none of --od and --of'
    sys.exit(0)
# check verbosity options are correct
if not (verbosity in [1, 2, 3]):
    print 'verbositiy can only be 1 or 2. default is 1 (most verbose)'
    sys.exit(0)
# configuring global program variables
if root_direct == False:
    root_direct = os.getcwd()
# figuring out where to write output
if out_direct == False and out_fn == False:
    out_fn = join(os.getcwd(),'pr_testing_'+time.ctime().replace(' ','_').replace(':','_')+'.log')
elif out_direct != False:
    out_fn = join(out_direct,'pr_testing_'+time.ctime().replace(' ','_').replace(':','_')+'.log')
'''
else:
    # out_direct == False and out_fn != False
    # out_fn should be what we want
    pass
'''
if not os.path.isdir(out_direct):
    # create the directory
    os.mkdir(out_direct)

if do_lin:
    # check if output image directory exists... if not create it...
    image_direct = join(out_direct, 'error_maps_rate')
    if not os.path.isdir(image_direct):
        # create the directory
        os.mkdir(image_direct)

if do_ts:
    # check if output image directory exists... if not create it...
    image_direct2 = join(out_direct, 'error_maps_ts')
    if not os.path.isdir(image_direct2):
        # create the directory
        os.mkdir(image_direct2)

print 'writing output to '+out_fn
out_fp = open(name=out_fn, mode='w')
out_fp.write('cmp_out.py log file\n')
out_fp.write('===================\n')

# these used in below loop to keep track of what indent writing to in log file
in_tst = ''
in_dat = ''
write_tst = False
write_dat = False

sorted_nums = sorted(filter(lambda x: os.path.isdir(os.path.join(root_direct, x)), os.listdir(root_direct)))
#print sorted_nums
#while True: pass
dat_fld = os.path.split(os.path.split(root_direct)[0])[1]
tst_fld = os.path.split(os.path.split(os.path.split(root_direct)[0])[0])[1]
out_fp.write('* '+tst_fld+'\n')
out_fp.write('\t'*1+'* '+dat_fld+'\n')
for cont_fld in sorted_nums:
    path = join(root_direct, cont_fld)

    pick_f = open(join(path, 'python', 'out_py.pkl'), 'r')      # todo: don't create pickle. read in tiff files (this requires changing run_pyrate.py though)
    out_py = pickle.load(pick_f)
    pick_f.close()
    out_mt = scipy.io.loadmat(join(path, 'matlab', 'out_mt.mat'))

    # load reference pixel values
    refpix_f = open(join(path, 'python', 'refpt.pkl'), 'r')
    rp_py = pickle.load(refpix_f)
    refpix_f.close()
    rp_mt = scipy.io.loadmat(join(path, 'matlab', 'refpt.mat'))
    rp_mt = rp_mt['refpt']
    rp_mt = rp_mt[0,0]
    rp_mt = (rp_mt[0].item()-1, rp_mt[1].item()-1)  # matlab reference point seems to be greater by 1 in each dimension

    # paths where write out error image file
    # overlay ontop of matlab image
    #er_img_tsincr_fn = join(path, 'er_tsincr.png')     # do it just for rate for now
    #er_img_rate_fn = join(path, 'er_rate.png')

    ers = []
    if ('tsincr' in out_py) and ('tsincr' in out_mt):
        py = out_py['tsincr']   # data is double
        mt = out_mt['tsincr']   # data is single
        #py = py[:py.shape[0]-1,:,:]     # remove last row... have to check this happens for all inputs
        ers.append(('tsincr', chk_out_mats(m1=py, m2=mt)))

    if ('rate' in out_py) and ('stackmap' in out_mt):
        py = out_py['rate']   # data is double
        mt = out_mt['stackmap']   # data is single
        #py = py[:py.shape[0]-1,:]     # remove last row... have to check this happens for all inputs
        ers.append(('rate/stackmap', chk_out_mats(m1=py, m2=mt)))

    # check if actually have an error to write
    if ['', ''] != [er[1] for er in ers]:
        # write indents if need to
        #if write_tst:
        #    out_fp.write('* '+tst_fld+'\n')
        #if write_dat:
        #    out_fp.write('\t'*1+'* '+dat_fld+'\n')
        out_fp.write('\t'*2+'* '+cont_fld+'\n')
        # write errors
        for er in ers:
            if er[1] != '':
                out_fp.write('\t'*3+'* '+er[0]+'\n')
                out_fp.write(er[1])     # <-- pass an indent level parameter

# write out end of run metrics (n_ref_fail, etc.)
out_fp.write('\nSET METRICS')
out_fp.write('\n======================\n')
# -----------------------------------------------------
met_string = ''
tot_tested = n_no_ref_fail + n_rate_ref + n_ts_ref
met_string += '* '+str(tot_tested/2)+' equivalent configuration files tested = '+str(tot_tested)+' rate & ts'+'\n'
temp_val = (float(n_rate_ref)/(tot_tested/2))*100
met_string += '* reference fails / # equiv. configurations = '+str(n_rate_ref)+'/'+str(tot_tested/2)+' = '+str(temp_val)+'%\n'
# ------------------------------
met_string += '======================\n'
temp_val = (float(len(hist_rate_tol))/n_rate_total)*100
met_string += '* rate tol fails / total rate (minus rate ref fails) = '+str(len(hist_rate_tol))+'/'+str(n_rate_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_rate_nan1))/n_rate_total)*100
met_string += '* rate nan1 fails / total rate (minus rate ref fails) = '+str(len(hist_rate_nan1))+'/'+str(n_rate_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_rate_nan2))/n_rate_total)*100
met_string += '* rate nan2 fails / total rate (minus rate ref fails) = '+str(len(hist_rate_nan2))+'/'+str(n_rate_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_rate_nans))/n_rate_total)*100
met_string += '* rate nan1 or nan2 fails / total rate (minus rate ref fails) = '+str(len(hist_rate_nans))+'/'+str(n_rate_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_rate_tot))/n_rate_total)*100
met_string += '* rate tolerance or nan1 or nan2 (any) fails / total rate (minus rate ref fails) = '+str(len(hist_rate_tot))+'/'+str(n_rate_total)+' = '+str(temp_val)+'%\n'
# ------------------------------
met_string += '======================\n'
temp_val = (float(len(hist_ts_tol))/n_ts_total)*100
met_string += '* ts tol fails / total ts (minus ts ref fails) = '+str(len(hist_ts_tol))+'/'+str(n_ts_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_ts_nan1))/n_ts_total)*100
met_string += '* ts nan1 fails / total ts (minus ts ref fails) = '+str(len(hist_ts_nan1))+'/'+str(n_ts_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_ts_nan2))/n_ts_total)*100
met_string += '* ts nan2 fails / total ts (minus ts ref fails) = '+str(len(hist_ts_nan2))+'/'+str(n_ts_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_ts_nans))/n_ts_total)*100
met_string += '* ts nan1 or nan2 fails / total ts (minus ts ref fails) = '+str(len(hist_ts_nans))+'/'+str(n_ts_total)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(hist_ts_tot))/n_ts_total)*100
met_string += '* ts tolerance or nan1 or nan2 (any) fails / total ts (minus ts ref fails) = '+str(len(hist_ts_tot))+'/'+str(n_ts_total)+' = '+str(temp_val)+'%\n'
# ------------------------------
n1 = hist_ts_tol+hist_rate_tol
n2 = hist_ts_nan1+hist_rate_nan1
n3 = hist_ts_nan2+hist_rate_nan2
n4 = hist_ts_nans+hist_rate_nans
n5 = hist_ts_tot+hist_rate_tot
# ------------------------------
d1 = n_ts_total+n_rate_total        # denominator. all the same

met_string += '======================\n'
temp_val = (float(len(n1))/d1)*100
met_string += '* all (rate or ts) tol fails / total rate & ts (minus rate & ts ref fails) = '+str(len(n1))+'/'+str(d1)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(n2))/d1)*100
met_string += '* all (rate or ts) nan1 fails / total rate & ts (minus rate & ts ref fails) = '+str(len(n2))+'/'+str(d1)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(n3))/d1)*100
met_string += '* all (rate or ts) nan2 fails / total rate & ts (minus rate & ts ref fails) = '+str(len(n3))+'/'+str(d1)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(n4))/d1)*100
met_string += '* all (rate or ts) nan1 or nan2 fails / total rate & ts (minus rate & ts ref fails) = '+str(len(n4))+'/'+str(d1)+' = '+str(temp_val)+'%\n'
temp_val = (float(len(n5))/d1)*100
met_string += '* all (rate or ts) tolerance or nan1 or nan2 (any) fails / total rate & ts (minus rate & ts ref fails) = '+str(len(n5))+'/'+str(d1)+' = '+str(temp_val)+'%\n'

out_fp.write(met_string)

# everything that goes in the logfile is done. close it
out_fp.close()

if do_hist:
    # plot histograms... store in /histogram folder of out_direct
    histo_d = join(out_direct, 'histograms')
    if not os.path.isdir(histo_d):
        os.mkdir(histo_d)

    plt.close()     # close overwritten ratemap figure...

    # rate
    fig2 = plt.figure()
    plt.hist(hist_rate_tol, 100, range=[0, 100])
    plt.title('Linear rate tolerance errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_rate_total)+' equiv. configurations (minus ref. mismatches)')
    fig2.savefig(join(histo_d, 'rate_tol.png'), bbox_inches='tight')
    plt.close()

    fig3 = plt.figure()
    plt.hist(hist_rate_nan1, 100, range=[0, 100])
    plt.title('Linear rate NaN1 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_rate_total)+' equiv. configurations (minus ref. mismatches)')
    fig3.savefig(join(histo_d, 'rate_nan1.png'), bbox_inches='tight')
    plt.close()

    fig4 = plt.figure()
    plt.hist(hist_rate_nan2, 100, range=[0, 100])
    plt.title('Linear rate NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_rate_total)+' equiv. configurations (minus ref. mismatches)')
    fig4.savefig(join(histo_d, 'rate_nan2.png'), bbox_inches='tight')
    plt.close()

    fig5 = plt.figure()
    plt.hist(hist_rate_nans, 100, range=[0, 100])
    plt.title('Linear rate NaN1 or NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_rate_total)+' equiv. configurations (minus ref. mismatches)')
    fig5.savefig(join(histo_d, 'rate_nans.png'), bbox_inches='tight')
    plt.close()

    fig6 = plt.figure()
    plt.hist(hist_rate_tot, 100, range=[0, 100])
    plt.title('Linear rate tolerance or NaN1 or NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_rate_total)+' equiv. configurations (minus ref. mismatches)')
    fig6.savefig(join(histo_d, 'rate_tot.png'), bbox_inches='tight')
    plt.close()

    # ts algorithm
    fig7 = plt.figure()
    plt.hist(hist_ts_tol, 100, range=[0, 100])
    plt.title('Time series tolerance errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_ts_total)+' equiv. configurations (minus ref. mismatches)')
    fig7.savefig(join(histo_d, 'ts_tol.png'), bbox_inches='tight')
    plt.close()

    fig8 = plt.figure()
    plt.hist(hist_ts_nan1, 100, range=[0, 100])
    plt.title('Time series NaN1 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_ts_total)+' equiv. configurations (minus ref. mismatches)')
    fig8.savefig(join(histo_d, 'ts_nan1.png'), bbox_inches='tight')
    plt.close()

    fig9 = plt.figure()
    plt.hist(hist_ts_nan2, 100, range=[0, 100])
    plt.title('Time series NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_ts_total)+' equiv. configurations (minus ref. mismatches)')
    fig9.savefig(join(histo_d, 'ts_nan2.png'), bbox_inches='tight')
    plt.close()

    fig10 = plt.figure()
    plt.hist(hist_ts_nans, 100, range=[0, 100])
    plt.title('Time series NaN1 or NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_ts_total)+' equiv. configurations (minus ref. mismatches)')
    fig10.savefig(join(histo_d, 'ts_nans.png'), bbox_inches='tight')
    plt.close()

    fig11 = plt.figure()
    plt.hist(hist_ts_tot, 100, range=[0, 100])
    plt.title('Time series tolerance or NaN1 or NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of '+str(n_ts_total)+' equiv. configurations (minus ref. mismatches)')
    fig11.savefig(join(histo_d, 'ts_tot.png'), bbox_inches='tight')
    plt.close()

    # both algorithms
    fig12 = plt.figure()
    plt.hist(hist_ts_tol+hist_rate_tol, 100, range=[0, 100])
    plt.title('All tolerance errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of all '+str(n_rate_total+n_ts_total)+' rate & time series tests (minus ref. mismatches)')
    fig12.savefig(join(histo_d, 'both_tol.png'), bbox_inches='tight')
    plt.close()

    fig13 = plt.figure()
    plt.hist(hist_ts_nan1+hist_rate_nan1, 100, range=[0, 100])
    plt.title('All NaN1 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of all '+str(n_rate_total+n_ts_total)+' rate & time series tests (minus ref. mismatches)')
    fig13.savefig(join(histo_d, 'both_nan1.png'), bbox_inches='tight')
    plt.close()

    fig14 = plt.figure()
    plt.hist(hist_ts_nan2+hist_rate_nan2, 100, range=[0, 100])
    plt.title('All NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of all '+str(n_rate_total+n_ts_total)+' rate & time series tests (minus ref. mismatches)')
    fig14.savefig(join(histo_d, 'both_nan2.png'), bbox_inches='tight')
    plt.close()

    fig15 = plt.figure()
    plt.hist(hist_ts_nans+hist_rate_nans, 100, range=[0, 100])
    plt.title('All NaN1 or NaN2 errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of all '+str(n_rate_total+n_ts_total)+' rate & time series tests (minus ref. mismatches)')
    fig15.savefig(join(histo_d, 'both_nans.png'), bbox_inches='tight')
    plt.close()

    fig16 = plt.figure()
    plt.hist(hist_ts_tot+hist_rate_tot, 100, range=[0, 100])
    plt.title('All tolerance or NaN1 or NaN2 (any) errors histogram')
    plt.xlabel('Percentage error')
    plt.ylabel('Quantity out of all '+str(n_rate_total+n_ts_total)+' rate & time series tests (minus ref. mismatches)')
    fig16.savefig(join(histo_d, 'both_tot.png'), bbox_inches='tight')
    plt.close()

    #plt.show()