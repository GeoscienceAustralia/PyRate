'''
rewritten cfg file generator script
-----------------------------------
* PR_TESTING_DATA environment variable to be set
    * data will be copied over to cfg folders (runs faster if each process has its own data)
* notation:
    * <>_d = directory
    * <>_fn = filename (that is not a directory (i.e. a regular file))
    * <>_<>_r = relative (i.e. not full path)
* todo: VALUE should really by VALUES
'''

# imports ------------------------
import sys
import os
import shutil
from numpy import arange
from subprocess import call
from copy import deepcopy
from os.path import join
# --------------------------------

# SELECT WHAT DATA YOU WANT TO TEST ===================================================================
# =====================================================================================================
# check environment variable is set
if not ('PR_TESTING_DATA' in os.environ):
    print 'you need to set your PR_TESTING_DATA environment variable'
    sys.exit(0)
dat_root_d = os.environ['PR_TESTING_DATA']     # check at start of script that this exists
dat_drs = ['syd_g']                            # manually select specific folders

# PARAM SWEEPS (YOU MODIFY THESE TO WHAT YOU WANT =====================================================
# =====================================================================================================
base = [
    ['networkx_or_matlab', [0]],
    ['luigi', [0]],

    ['atmfit', [0]],     # Matt says don't use this
                         # normally would have to specify atmfitmethod, atmfitlksx... afterward

    ['basepflag', [0]],  # needed for Matlab to work, Python will ignore it
    ['unwflag', [0]],    # not used in Python. leave as 0
    ['prjflag', [3]],    # Re-project data from Line of sight, 1 = vertical, 2 = horizontal, 3 = no conversion
                         # 3 or run into problems. try others later

    ['noDataAveragingThreshold', [0.5]],    # Sudipta says this is for network... doesn't look like it's even used

    ['noDataValue', [0.0]],     # depends on data. convert these values to NaNs
    ['nan_conversion', [1]],    # if you have a noDataValue convert it to NaN

    ['ifglks', [(1, 1)]],       # tuples for x and y looks

    # stacking/velocity parameters
    ['nsig', [2]],      # Matlab default 2
    ['pthr', [12]],     # dependant on ifgsize
    ['maxsig', [2]],    # Matlab default 2

    ['vcmtmethod', [1]],    # Python uses general method 1... keep it at this
    ['vcmsmethod', [1]],    # is calculated but not even used, only one that has looks, defaults to 10 looks
    ['vcmslksx', [10]],     # what to set for matlab
    ['vcmslksy', [10]],
]

#  BRANCHES -----------------------------------------------
# ---------------------------------------------------------
crop_b = [
    [
        ['ifgcropopt', [1]],        # all ifgs the same size. don't need to do any processing
    ],
    #[
    #    ['ifgcropopt', [1, 2]],
    #],
    #[
    #    ['ifgcropopt', [3]],
    #    ['ifgxfirst', [3]],
    #    ['ifgxlast', [3]],
    #    ['ifgyfirst', [3]],
    #    ['ifgylast', [3]],
    #],
]

orb_b = [
    [
        ['orbfit', [0]],            # on or off
    ],
    #[
    #    ['orbfit', [1]],
    #    ['orbfitmethod', [1, 2]],
    #    ['orbfitdegrees', [1]],    # default 1
    #    ['orbfitlks', [(5, 5)]],       # will probably depend on data
    #],
]

ref_b = [
    [
        ['refest', [2]],            # 1 is stupid
        ['refx', [0]],
        ['refy', [0]],
        ['refnx', [10]],            # x/10
        ['refny', [10]],            # y/10
        ['refchipsize', [5]],       # ?
        ['refminfrac', [0.8]],      # default for matlab
    ],
]

ts_b = [
    #[
    #    ['tscal', [0]],        # will never test this (will always want to create time series)
    #],
    [   # laplcian smoothing
        ['tscal', [1]],
        ['tsmethod', [1]],

        ['smorder', [1]],       # pyrate config.py says this and smfactor "NOT CURRENTLY USED"
        ['smfactor', [1]],

        ['ts_pthr', [12]],      # based on ifg size (12 is for sydney data)

        ['ts_interp', [0]],     # Matlab laplacian only. default 0. Matt says keep this at 0

        ['stepest', [0]],       # Matlab only, default 0
    ],
    [
        # SVD (note pyrate required smorder & smfactor parameters (to run), even though they won't be used ==> dummy param
        ['tscal', [1]],
        ['tsmethod', [2]],

        ['smorder', [0]],
        ['smfactor', [0]],

        ['ts_pthr', [12]],      # based on ifg size (12 is for sydney data)

        ['stepest', [0]],       # matlab only, default 0
    ],
]

ML, PY = 0, 1                   # to indicate Matlab or Python
NAME, VALUE, STATE = 0, 1, 2    # used in recursive function to generate cfg files

# combine branches to create work
work = []
for i1 in crop_b:
    for i2 in orb_b:
        for i3 in ref_b:
            for i4 in ts_b:
                work.append(deepcopy(base + i1 + i2 + i3 + i4))

# set state to 0
for b_comb in work:     # b_com = branch combination, which is what it is
    [n_v.append(0) for n_v in b_comb]

# figure out how many permutations (and hence how many 0's we need)
n_perms = 0
for b_comb in work:
    n_perms_per_b = 1
    for n_v_s in b_comb:
        n_perms_per_b *= len(n_v_s[1])
    n_perms += n_perms_per_b

print 'generating '+str(n_perms)+' configuration files...'
n_zs = len(str(n_perms))

# FUNCTIONS ===============================================
# =========================================================
def rec_loop(n_v_s, loop=0):
    if loop != 0:
        if n_v_s[loop][STATE] == len(n_v_s[loop][VALUE]):
            n_v_s[loop][STATE] = 0
    if loop < (len(n_v_s)-1):
        while n_v_s[loop][STATE] < len(n_v_s[loop][VALUE]):
            rec_loop(n_v_s, loop+1)
            n_v_s[loop][STATE] += 1
    if loop == (len(n_v_s)-1):
        while (n_v_s[loop][STATE] < len(n_v_s[loop][VALUE])):
            global counter

            z_str = str(counter).zfill(n_zs)

            # some helper references
            num_d = join(script_d, dat_dr, z_str)      # the number directory
            ml_d = join(num_d, 'matlab')
            py_d = join(num_d, 'python')

            # check if the following directories exist. if not create them
            ds = []
            ds.append(join(script_d, dat_dr))    # create folder with name of actual test data (used for lumping together tests on data)
            ds.append(num_d)
            ds.append(ml_d)
            ds.append(py_d)
            ds.append(join(ml_d, 'logs'))
            for d in ds:
                if not os.path.exists(d):
                    os.mkdir(d)

            # copy the data over for python and matlab (each cfg gets its own data files)
            ml_dat_d = join(ml_d, 'data')      # todo: require this above. make reference
            py_dat_d = join(py_d, 'data')
            shutil.copytree(dat_d, ml_dat_d)    # docs say don't need to create the directory beforehand
            shutil.copytree(dat_d, py_dat_d)

            # some helper references...
            ml_obs_d = join(ml_dat_d, 'obs')
            py_obs_d = join(py_dat_d, 'obs')
            ml_dem_d = join(ml_dat_d, 'dem')
            py_dem_d = join(py_dat_d, 'dem')

            cfg = []   # 2 x n list for storing n configuration strings for ml and py configuration files

            # set obsdir, demdir, ifgfilelist appropriately (header 1)
            h1 = ['', '']
            h1[ML] += '### FILE DIRECTORIES ###\n'
            h1[ML] += 'obsdir: '+ml_obs_d+os.path.sep+'\n'
            h1[ML] += 'demdir: '+ml_dem_d+os.path.sep+'\n'
            h1[ML] += 'outdir: '+ml_d+os.path.sep+'\n'
            # -----------------------------------------------------------------------------------------
            h1[PY] += '### FILE DIRECTORIES ###\n'
            h1[PY] += 'obsdir: '+py_obs_d+os.path.sep+'\n'
            h1[PY] += 'demdir: '+py_dem_d+os.path.sep+'\n'
            h1[PY] += 'outdir: '+py_d+os.path.sep+'\n'
            cfg.append(h1)

            h2 = ['', '']
            h2[ML] += '### OTHER ###\n'
            h2[ML] += 'ifgfilelist: '+ifgfilelist_fr+'\n'                   # Matlab expects relative path (to obs directory)
            h2[ML] += 'processor: 1\n'
            h2[ML] += 'demfile: '+demfile_fr+'\n'
            h2[ML] += 'parfile: '+parfile_fr+'\n'
            # -----------------------------------------------------------------------------------------
            h2[PY] = '### OTHER ###\n'
            h2[PY] += 'ifgfilelist: '+join(py_obs_d, ifgfilelist_fr)+'\n'   # Python wants absolute path
            h2[PY] += 'processor: 1\n'
            h2[PY] += 'demfile: '+join(py_dem_d, demfile_fr)+'\n'
            h2[PY] += 'demHeaderFile: '+join(py_dem_d, parfile_fr)+'\n'     # different name for Python
            cfg.append(h2)

            h3 = ['', '']
            h3[ML] += '### NEEDED FOR MATLAB TO WORK ###'
            h3[ML] += 'wavelength: '+wl+'\n'
            h3[ML] += 'gmtfaultfile: xxxx.gmt\n'
            h3[ML] += 'ifgname_prefix_len: 0\n'
            h3[ML] += 'ifgname_date_format: 8\n'
            h3[ML] += 'auxdir: logs'+os.path.sep+'\n'
            # below stuff isn't used but is needed (doesn't matter if directory is created or not)
            h3[ML] += 'ratedir: ratemap'+os.path.sep+'\n'
            h3[ML] += 'tsdir: timeseries'+os.path.sep+'\n'
            h3[ML] += 'profdir: prof'+os.path.sep+'\n'
            # -----------------------------------------------------------------------------------------
            # this header isn't relevant for Python. leave as ''
            cfg.append(h3)

            # stringify the parameters that are recursed in this function
            body = ['', '']
            for item in n_v_s:
                # check if name is ifglks or orbfitlks and set appropriately
                if item[NAME] == 'ifglks':
                    body[ML] += 'ifglksx: '+str(item[VALUE][item[STATE]][0])+'\n'
                    body[ML] += 'ifglksy: '+str(item[VALUE][item[STATE]][1])+'\n'
                elif item[NAME] == 'orbfitlks':
                    body[ML] += 'orbfitlksx: '+str(item[VALUE][item[STATE]][0])+'\n'
                    body[ML] += 'orbfitlksy: '+str(item[VALUE][item[STATE]][1])+'\n'
                # other data doesn't need processing
                else:
                    body[ML] += item[NAME]+': '+str(item[VALUE][item[STATE]])+'\n'
                    body[PY] += item[NAME]+': '+str(item[VALUE][item[STATE]])+'\n'
            cfg.append(body)

            # finished. write the two cfg files with diff extensions\
            mt_cfg_fp = open(join(num_d, z_str)+'.mcfg', 'w')
            py_cfg_fp = open(join(num_d, z_str)+'.pcfg', 'w')
            for block in cfg:
                mt_cfg_fp.write(block[ML])
                py_cfg_fp.write(block[PY])
            mt_cfg_fp.close()
            py_cfg_fp.close()

            counter += 1
            n_v_s[loop][STATE] += 1 # this must be at bottom of while loop

# [END OF] FUNCTIONS ======================================
# =========================================================
# create a folder of the same name as this script in the same directory as the script
test_d, script_fr = os.path.split(os.path.realpath(__file__))    # will give the directory script is in
script_d = join(test_d, script_fr[:-3])
if not os.path.exists(script_d):
        os.mkdir(script_d)

for dat_dr in dat_drs:
    # start by iterating through each data folder...
    # set references that can be used in rec_loop
    # name of ifgfilelist (has to be in obs_d (convention))
    dat_d = join(dat_root_d, dat_dr)
    obs_d = join(dat_d, 'obs')
    dem_d = join(dat_d, 'dem')

    for fr in os.listdir(obs_d):
        if fr[:4] == 'ifms':        # convention. ifg file list has to begin with 'ifms'. e.g. ifms_17
            ifgfilelist_fr = fr     # set here because i may use 'fr' again
            break
    else:
        print 'ifg file list could not be found. must be in "obs" directory and must begin with "ifms" (e.g. ifms_17)'
        sys.exit(0)

    # get the relative demfile name
    for fr in os.listdir(dem_d):
        if fr[-4:] == '.dem':       # convention
            demfile_fr = fr         # set here because i may use 'fr' again
            break
    else:
        print '<demfile error message>'
        sys.exit(0)

    # get the relative parfile name
    for fr in os.listdir(dem_d):
        if fr[-4:] == '.par':       # convention
            parfile_fr = fr         # set here because i may use 'fr' again
            break
    else:
        print '<parfile error message>'
        sys.exit(0)

    # get wavelength from special file
    # todo: read from par file
    wl_fp = open(join(dat_d, 'wavelength.txt'), 'r')
    wl = wl_fp.readline()
    wl_fp.close()

    counter = 0
    for b_comb in work:
        rec_loop(n_v_s=b_comb)
        # set state back to 0 for next folder
        for n_v_s in b_comb:
            n_v_s[STATE] = 0

