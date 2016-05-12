'''
USAGE:
------
python run_cfgs.py [-m, -b] <full directory path>
e.g.

what this script does
---------------------
* recurses the <full directory path> and searches for *.pcfg and *.mcfg files and runs them (with PyRate and Matlab, respectively)
    * it makes sense to run gen_cfgs.py first, but this will work on any directory structure containing *.pcfg and *.mcfg files
* works on both unix and windows
    * you will need the "Matlab Compiled Runtimes" to run the Matlab scripts
        * you are unable to install these on GA windows workstations (see below)
        * <instructions...>

GA MS windows problems
----------------------
* since you can't install the runtimes on GA Windows workstations and option is to run this script on unix (and it only
  makes sense to do this once because the output won't change (you can't and shouldn't need to modify the Matlab program)
* there are command line switches for

Requirements
------------
* PYRATEPATH needs to be set
* ensure matlab runtimes environment variables are set

'''

import sys
import os
import threading
import subprocess
import time

from os.path import join

from copy import deepcopy

N_PROCS = 10    # max number of procs want to allow. TODO: make this adaptive
procs = []

# use getop now, since we need switches
# make default just PyRate (-b for both, -m for just matlab)
import getopt
opts, args = getopt.getopt(sys.argv[1:], 'bm')
if len(args) != 1:
    print 'must supply just 1 argument: the directory structure you wish to run'
    sys.exit(0)
if not (opts == [] or len(opts) == 1):
    print 'supply no switches or either one of -m or -b'
    sys.exit(0)
if opts == []:
    p_only = True
    m_only = False
    both = False
else:
    p_only = False
    if opts[0][0] == '-m':
        m_only = True
        both = False
    else:
        m_only = False
        both = True

# make it so it walks from this and runs anything in the walk
root_d = sys.argv[1:][0]    # root directory. will walk from here
if not os.path.isdir(root_d):
    print 'argument not a valid directory'
    sys.exit(0)

if not ('PYRATEPATH' in os.environ):
    print 'you need to set your PYRATEPATH'
    sys.exit(0)

PYRATEPATH = os.environ['PYRATEPATH']

# configuration for running matlab...
# hard code in for now
mcc_exec_pth = '/nas/gemd/insar/pr_testing/mcc/v1/pirate'
mcc_env = deepcopy(os.environ)
mcc_env['LD_LIBRARY_PATH'] = "/nas/gemd/insar/pr_testing/MATLAB_Runtime/v85/runtime/glnxa64:/nas/gemd/insar/pr_testing/MATLAB_Runtime/v85/bin/glnxa64:/nas/gemd/insar/pr_testing/MATLAB_Runtime/v85/sys/os/glnxa64"

PR, ML = 0, 1

def spawner(spawn_type, prog):
    global spawner_fin, procs  # modify these so global
    for path, dirs, files in os.walk(root_d):
        cont_fld = os.path.split(path)[1]    # containing folder
        if cont_fld == 'python' or cont_fld == 'matlab':
            continue

        for file in files:
            #print file[-4:]
            if (file[-5:] == '.pcfg' and spawn_type == PR and (p_only or both)) or (file[-5:] == '.mcfg' and spawn_type == ML and (m_only or both)):     # TODO: make this check better
                # we want to run this config file. but only if there is space in the queue to do so
                while len(procs) >= N_PROCS:
                    time.sleep(1)   # TODO: wat set this to
                # since this is the only spawner, queue won't be refilled

                if spawn_type == PR:
                    procs.append(subprocess.Popen(['python', join(PYRATEPATH,'pyrate','scripts', prog), join(path,file)]))
                else:   # run compiled matlab
                    pass    # for now
    spawner_fin = True

def oversee():
    while True:
        if procs != []:
            for proc in procs:
                # poll the proc... just to see what we get... and if the object stays
                if proc.poll() == None:     # process hasn't returned yet
                    time.sleep(0.1)     # may not need to sleep here
                else:
                    print proc.poll()   # better to assign to variable
                    # inspect return code, but for now just remove from list
                    procs.remove(proc)
        elif spawner_fin:
            break

spawner_fin = False
thr_oversee = threading.Thread(target=oversee)
thr_spawner = threading.Thread(target=spawner, kwargs={'spawn_type': PR, 'prog': 'run_prepifg.py'})
thr_oversee.start()
thr_spawner.start()
thr_oversee.join()
thr_spawner.join()

# ML prepifg here

spawner_fin = False
thr_oversee = threading.Thread(target=oversee)
thr_spawner = threading.Thread(target=spawner, kwargs={'spawn_type': PR, 'prog': 'run_pyrate.py'})
thr_oversee.start()
thr_spawner.start()
thr_oversee.join()
thr_spawner.join()