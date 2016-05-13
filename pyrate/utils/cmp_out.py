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
def chk_3darrs(a1, a2):
    # returns an error string that you just plonk into the error file
    er_str = ''

    if a1.shape != a2.shape:
        return ' '*12+'* ARRAY SHAPES DO NOT MATCH. DOES NOT MAKE SENSE TO CHECK FURTHER.\n'

    # otherwise test individual elements
    #relative_tolerance=1e-10
    relative_tolerance=0.1
    n_r, n_c, n_e = a1.shape        # number of rows, cols, epochs (a2 and py should be equal now)
    it_r, it_c, it_e = (0,0,0)      # iterators. not practice, but whatever
    while it_e < n_e:
        it_r = 0
        while it_r < n_r:
            it_c = 0
            while it_c < n_c:
                it = (it_r, it_c, it_e)     # this can be used to index ndarray
                if (not np.isnan(a1[it])) and (not np.isnan(a2[it])):  # both floats
                    # Duncan's way of checking if floats are equal...
                    if abs(a1[it].item() - a2[it].item()) > relative_tolerance:
                        # not within tolerance...
                        er_str += ' '*12+'* do not match to tolerance @ '+str(it)+'\n'
                        er_str += ' '*16+'* a1 = '+str(Decimal(a1[it].item()))+'\n'
                        er_str += ' '*16+'* a2 = '+str(Decimal(a2[it].item()))+'\n'
                    '''
                    # my original ideas about checking if floats are equal...
                    # check if match to a certain number of decimal places
                    if truncate(a1[it].item(),1) != truncate(a2[it].item(),1):
                        er_str += ' '*8+'* do not match to precision @ '+str(it)+'\n'
                        er_str += ' '*12+'* a1 = '+str(Decimal(a1[it].item()))+'\n'
                        er_str += ' '*12+'* a2 = '+str(Decimal(a2[it].item()))+'\n'
                    '''
                else:   # atleast one of py[it], a2[it] is NaN
                    if (not np.isnan(a1[it])) and (np.isnan(a2[it])):
                        er_str += ' '*12+'* a1 should be NaN (but is not)'+'\n'
                    if (np.isnan(a1[it])) and (not np.isnan(a2[it])):
                        er_str += ' '*12+'* a1 should not be NaN (but is)'+'\n'
                it_c += 1
            it_r += 1
        it_e += 1
    return er_str
# =============================================================================
# =============================================================================
def chk_2darrs(a1, a2):
    # returns an error string that you just plonk into the error file
    er_str = ''

    if a1.shape != a2.shape:
        return ' '*12+'* ARRAY SHAPES DO NOT MATCH. DOES NOT MAKE SENSE TO CHECK FURTHER.\n'

    # otherwise test individual elements
    #relative_tolerance=1e-10
    relative_tolerance=0.1
    n_r, n_c = a1.shape        # number of rows, cols, epochs (a2 and py should be equal now)
    it_r, it_c = (0,0)      # iterators. not practice, but whatever
    while it_r < n_r:
        it_c = 0
        while it_c < n_c:
            it = (it_r, it_c)     # this can be used to index ndarray
            if (np.isnan(a1[it])) and (np.isnan(a2[it])):
                print 'GOT NAN MATCH'
            if (not np.isnan(a1[it])) and (not np.isnan(a2[it])):  # both floats
                # Duncan's way of checking if floats are equal...
                if abs(a1[it].item() - a2[it].item()) > relative_tolerance:
                    # not within tolerance...
                    er_str += ' '*12+'* do not match to tolerance @ '+str(it)+'\n'
                    er_str += ' '*16+'* PyRate = '+str(Decimal(a1[it].item()))+'\n'
                    er_str += ' '*16+'* Matlab = '+str(Decimal(a2[it].item()))+'\n'
                '''
                # my original ideas about checking if floats are equal...
                # check if match to a certain number of decimal places
                if truncate(a1[it].item(),1) != truncate(a2[it].item(),1):
                    er_str += ' '*8+'* do not match to precision @ '+str(it)+'\n'
                    er_str += ' '*12+'* a1 = '+str(Decimal(a1[it].item()))+'\n'
                    er_str += ' '*12+'* a2 = '+str(Decimal(a2[it].item()))+'\n'
                '''
            else:   # atleast one of py[it], a2[it] is NaN
                if (not np.isnan(a1[it])) and (np.isnan(a2[it])):
                    er_str += ' '*12+'* PyRate should be NaN (but is not) @ '+str(it)+'\n'
                if (np.isnan(a1[it])) and (not np.isnan(a2[it])):
                    er_str += ' '*12+'* PyRate should not be NaN (but is) @ '+str(it)+'\n'
                    er_str += ' '*16+'* PyRate = '+str(Decimal(a1[it].item()))+'\n'
                    er_str += ' '*16+'* Matlab = '+str(Decimal(a2[it].item()))+'\n'
            it_c += 1
        it_r += 1
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

class GFile(file):
    def __init__(self, file):
        file.__init__()

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
        if in_dat != dat_fld:
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
        ers.append(('tsincr', chk_3darrs(a1=py, a2=mt)))

        py = out_py['rate']   # data is double
        mt = out_mt['stackmap']   # data is single
        py = py[:py.shape[0]-1,:]     # remove last row... have to check this happens for all inputs
        ers.append(('rate/stackmap', chk_2darrs(a1=py, a2=mt)))

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
                    out_fp.write(' '*10+'* '+er[0]+'\n')
                    out_fp.write(er[1])     # <-- pass an indent level parameter

#out_fp.write('test!!!!')
out_fp.close()

