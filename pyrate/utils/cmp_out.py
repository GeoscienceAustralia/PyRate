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
#out_fp.write('test!!!!')
out_fp.close()