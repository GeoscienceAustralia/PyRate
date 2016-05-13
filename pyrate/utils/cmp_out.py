'''
* problems with getopt means won't have traditional command line "args"
* everything has to be an opt and param... if you leave out the opt. it does the default...
* options
    * directory
    * ignore
    * just
    * output file
'''

import sys
import os
from os.path import join
import time
import numpy as np
import pickle
import scipy.io     # scipy is the package. need must always import the module (io)
from decimal import Decimal

# =================================================================
# ARRAY CHECKING FUNCTIONS...
# these could really be made into just 1 function
# =================================================================
def chk_out_mats(m1, m2):
    # checks pirate/pyrate (m1/m2) output arrays
    # will determine if 2d or 3d... can only deal with either of those two
    # returns an error string that you just plonk into the error file
    relative_tolerance=0.1      # todo: have this as function and program argument

    er_str = ''
    if m1.shape != m2.shape:
        return ' '*12+'* ARRAY SHAPES DO NOT MATCH. DOES NOT MAKE SENSE TO CHECK FURTHER.\n'
    '''
    print len(m1.shape)
    while True: pass
    '''
    fun_mode = len(m1.shape)              # function mode. 2d or 3d matrices. don't bother checking otherwise. assume can't happen
    if fun_mode == 2:
        n_r, n_c = m1.shape
        n_e = 1
    else:   # gauranteed to be 3. else is faster than testing for 3
        n_r, n_c, n_e = m1.shape        # number of rows, cols, epochs

    it_r, it_c, it_e = (0,0,0)      # iterators. not pythonic, but whatever
    while it_e < n_e:
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
                        # not within tolerance...
                        er_str += ' '*16+'* do not match to tolerance @ '+str(it)+'\n'
                        er_str += ' '*20+'* m1 = '+str(Decimal(m1[it].item()))+'\n'
                        er_str += ' '*20+'* m2 = '+str(Decimal(m2[it].item()))+'\n'
                    '''
                    # my original ideas about checking if floats are equal...
                    # check if match to a certain number of decimal places
                    if truncate(m1[it].item(),1) != truncate(m2[it].item(),1):
                        er_str += ' '*8+'* do not match to precision @ '+str(it)+'\n'
                        er_str += ' '*12+'* m1 = '+str(Decimal(m1[it].item()))+'\n'
                        er_str += ' '*12+'* m2 = '+str(Decimal(m2[it].item()))+'\n'
                    '''
                else:   # atleast one of py[it], m2[it] is NaN
                    if (not np.isnan(m1[it])) and (np.isnan(m2[it])):
                        er_str += ' '*16+'* m1 should be NaN (but is not)'+'\n'
                    if (np.isnan(m1[it])) and (not np.isnan(m2[it])):
                        er_str += ' '*16+'* m1 should not be NaN (but is)'+'\n'
                it_c += 1
            it_r += 1
        it_e += 1
    return er_str
# =================================================================
# =================================================================

import getopt
opts, args = getopt.getopt(sys.argv[1:], 'd:i:j:', ['od=','of='])
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
# checking for valid output options
if out_direct != False and out_fn != False:
    print 'can only have one or the other or none of --od and --of'
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
print 'writing output to '+out_fn
out_fp = open(name=out_fn, mode='w')
out_fp.write('cmp_out.py log file\n')
out_fp.write('===================\n')

# these used in below loop to keep track of what indent writing to in log file
in_tst = ''
in_dat = ''
write_tst = False
write_dat = False

for path, dirs, files in os.walk(root_direct):
    rest, cont_fld = os.path.split(path)    # containing folder
    if cont_fld.isdigit():                  # if the containing folder is a number folder ==> will have stuff to compare...
        # check if we should skip this folder...
        rest, dat_fld = os.path.split(rest)
        _, tst_fld = os.path.split(rest)
        if jst_list != []:
            want_list = list(set(jst_list)-set(ign_list))   # todo: could do this more efficient
            if not (tst_fld in want_list):
                continue
        else:
            if tst_fld in ign_list:
                continue

        # determine if need to write indents
        write_tst = False
        write_dat = False
        if in_tst != tst_fld:
            in_tst = tst_fld
            write_tst = True
        if in_dat != dat_fld or write_tst:
            in_dat = dat_fld
            write_dat = True

        # good to go... check for errors on this num_fld
        # load files/data
        # todo: put these in a try and give a better error message
        pick_f = open(join(path, 'python', 'out_py.pkl'), 'r')
        out_py = pickle.load(pick_f)
        out_mt = scipy.io.loadmat(join(path, 'matlab', 'out_mt.mat'))

        ers = []
        py = out_py['tsincr']   # data is double
        mt = out_mt['tsincr']   # data is single
        py = py[:py.shape[0]-1,:,:]     # remove last row... have to check this happens for all inputs
        ers.append(('tsincr', chk_out_mats(m1=py, m2=mt)))

        py = out_py['rate']   # data is double
        mt = out_mt['stackmap']   # data is single
        py = py[:py.shape[0]-1,:]     # remove last row... have to check this happens for all inputs
        ers.append(('rate/stackmap', chk_out_mats(m1=py, m2=mt)))

        # check if actually have an error to write
        if ['', ''] != [er[1] for er in ers]:
            # write indents if need to
            if write_tst:
                out_fp.write('* '+tst_fld+'\n')
            if write_dat:
                out_fp.write(' '*4+'* '+dat_fld+'\n')
            out_fp.write(' '*8+'* '+cont_fld+'\n')
            # write errors
            for er in ers:
                if er[1] != '':
                    out_fp.write(' '*12+'* '+er[0]+'\n')
                    out_fp.write(er[1])     # <-- pass an indent level parameter

#out_fp.write('test!!!!')
out_fp.close()

