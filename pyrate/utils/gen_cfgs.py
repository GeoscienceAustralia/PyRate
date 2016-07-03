'''
v2
==
* TODO: maybe don't have PR_TESTING_DATA as an environment variable. just tell the user to modify it in the script?

v1
==
* more details in `PyRate/doc/garrick_report.docx`
* generate PyRate & Pirate configuration files based on selection of combinations of parameter 'sweeps'
* also gives each configuration its own copy of the data
* the data directory needs to be specified in the PR_TESTING_DATA environment variable

the testing procedure
---------------------
1. generate configuration files with with gen_cfgs.py
2. run those configurations with run_cfgs.py
3. compare the PyRate and Pirate outputs with cmp_out.py

how to use this script?
-----------------------
* copy the script into the directory where you want the generated configuration files stored
* select the data you want to use from line 52 or 53
* modify parameter sweeps (garrick_report.docx for explanation) from line 55 onwards
* perhaps rename the script reflecting the sweeps you selected
* running `$ python gen_cfgs.py` will produce the following directory structure:
     * ./ --|-- gen_cfgs.py
            |-- gen_cfgs --|-- data_direct_1 --|-- 0 --|-- PyRate and Pirate configuration files for configuration 0
                           |                   |       |-- data and output folders for configuration 0
                           |                   |-- 1 --|-- PyRate and Pirate configuration files for configuration 1
                           |                           |-- data and output folders for configuration 1
                           |-- data_direct_2 --|-- 0 --|-- PyRate and Pirate configuration files for configuration 0
                                               |       |-- data and output folders for configuration 0
                                               |-- 1 --|-- PyRate and Pirate configuration files for configuration 1
                                                       |-- data and output folders for configuration 1
'''

# IMPORTS =============================================================================================
# =====================================================================================================
import sys
import os
import shutil
from copy import deepcopy
from os.path import join

# SELECT WHAT DATA YOU WANT TO TEST ===================================================================
# =====================================================================================================
# check environment variable is set
if not ('PR_TESTING_DATA' in os.environ):
    print 'you need to set your PR_TESTING_DATA environment variable'
    sys.exit(0)
# TODO: check if is a valid directory
dat_root_d = os.environ['PR_TESTING_DATA']
dat_drs = ['syd_g']                            # manually select specific folders
#dat_drs = os.listdir(dat_root_d)               # select all folders in the

# PARAM SWEEPS (YOU MODIFY THESE TO WHAT YOU WANT) ====================================================
# =====================================================================================================
base = [
    ['networkx_or_matlab', [0]],
    ['luigi', [0]],

    ['atmfit', [0]],     # Matt says don't use this
                         # normally would have to specify atmfitmethod, atmfitlksx...

    ['basepflag', [0]],  # needed for Matlab to work. Python will ignore it
    ['unwflag', [0]],    # not used in Python. leave as 0
    ['prjflag', [3]],    # Re-project data from Line of sight, 1 = vertical,
                         # 2 = horizontal, 3 = no conversion
                         # 3 or run into problems. try others later

    ['noDataAveragingThreshold', [0.5]],    # Sudipta says this is for network...
                                            # doesn't look like it's even used

    ['noDataValue', [0.0]],     # depends on data. convert these values to NaNs
    ['nan_conversion', [1]],    # if you have a noDataValue convert it to NaN

    ['ifglks', [(1, 1)]],       # tuples for x and y looks

    # stacking/rate parameters
    ['nsig', [2]],      # Matlab default 2
    ['pthr', [12]],     # don't exceed number of epochs
    ['maxsig', [2]],    # Matlab default 2

    ['vcmtmethod', [1]],    # Python uses general method 1... keep it at this
    ['vcmsmethod', [1]],    # is calculated but not even used, only one that has looks
    ['vcmslksx', [10]],     # what to set for Matlab
    ['vcmslksy', [10]],
]

# BRANCHES ------------------------------------------------
# ---------------------------------------------------------
crop_b = [
    [
        ['ifgcropopt', [1]],        # all ifgs the same size (sydney test data)
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
    #    ['orbfitlks', [(5, 5)]],   # will probably depend on data
    #],
]

ref_b = [
    [
        ['refest', [2]],
        ['refx', [0]],
        ['refy', [0]],
        ['refn', [(7, 4)]],         # good values are res(x)/10, res(y)/10
                                    # above will de-tuple to refnx, refny
        ['refchipsize', [5]],       # can't have a:
                                    # * even chipsize (refpixel.py line 137),
                                    # * < 3,
                                    # * greater than width
        ['refminfrac', [0.8]],      # default for Matlab
    ],
]

ts_b = [
    #[
    #    ['tscal', [0]],        # will never test this (want to create time series)
    #],
    #[   # laplcian smoothing
    #    ['tscal', [1]],
    #    ['tsmethod', [1]],
    #
    #    ['smorder', [1]],       # PyRate says this and smfactor "NOT CURRENTLY USED"
    #                            # (but I think they are implemented)
    #                            # smorder: order of smoothing operator
    #                            # (1: first-order difference;
    #                            # 2: second-order difference)
    #    ['smfactor', [1]],      # smfactor: smoothing factor
    #                            # (0: calculate & plot L-curve;
    #                            # others: using the specific smoothing factor)
    #
    #    ['ts_pthr', [12]],      # don't exceed number of epochs
    #
    #    ['ts_interp', [0]],     # Matlab laplacian only. Matt says keep this at 0
    #
    #    ['stepest', [0]],       # Matlab only, default 0
    #],
    [   # SVD (note PyRate requires smorder & smfactor parameters (to run),
        # even though they won't be used
        ['tscal', [1]],
        ['tsmethod', [2]],

        ['smorder', [1]],       # aren't use for SVD, but Pyrate requires them (bug)
        ['smfactor', [1]],

        ['ts_pthr', [12]],      # don't exceed number of epochs

        ['stepest', [0]],       # Matlab only, default 0
    ],
]

ML, PY = 0, 1                   # to indicate Matlab or Python
NAME, VALUE, STATE = 0, 1, 2    # used in recursive function to generate config files

# combine branches to create work
work = []
for i1 in crop_b:
    for i2 in orb_b:
        for i3 in ref_b:
            for i4 in ts_b:
                work.append(deepcopy(base + i1 + i2 + i3 + i4))

# set state to 0
for b_comb in work:     # b_com = branch combination
    [n_v.append(0) for n_v in b_comb]

# figure out how many permutations (and hence how many 0's we need for file names)
n_perms = 0
for b_comb in work:
    n_perms_per_b = 1
    for n_v_s in b_comb:
        n_perms_per_b *= len(n_v_s[1])
    n_perms += n_perms_per_b

# inform user how many times data will be copied and ask if want to proceeed
print 'generating and copying data for **'+str(n_perms)+'** configuration file(s)...'
while True:
    ans = raw_input('Do you wish to continue (Y/n): ')
    if ans == 'Y':
        break
    if ans == 'n':
        sys.exit(0)

n_zs = len(str(n_perms))

# FUNCTIONS ===========================================================================================
# =====================================================================================================
def rec_loop(n_v_s, loop=0):
    '''
    * big recursive function to generate the combinations of configuration parameters
    * see later in script where the function is called to get an idea of the data structure this function will modify and work on
    :param n_v_s: stands for name, value, state
    :param loop: keeping track of what 'name'/parameter the loop is up to
    :return: Nothing
    '''
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

            # copy the data over for Python and Matlab (each cfg gets its own data files)
            ml_dat_d = join(ml_d, 'data')      # TODO: require this above. make reference
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
            h1[ML] += '\n'
            # -----------------------------------------------------------------------------------------
            h1[PY] += '### FILE DIRECTORIES ###\n'
            h1[PY] += 'obsdir: '+py_obs_d+os.path.sep+'\n'
            h1[PY] += 'demdir: '+py_dem_d+os.path.sep+'\n'
            h1[PY] += 'outdir: '+py_d+os.path.sep+'\n'
            h1[PY] += '\n'
            cfg.append(h1)

            h2 = ['', '']
            h2[ML] += '### OTHER ###\n'
            h2[ML] += 'ifgfilelist: '+ifgfilelist_fr+'\n'                   # Matlab expects relative path (to obs directory)
            h2[ML] += 'processor: 1\n'
            h2[ML] += 'demfile: '+demfile_fr+'\n'
            h2[ML] += 'parfile: '+parfile_fr+'\n'
            h2[ML] += '\n'
            # -----------------------------------------------------------------------------------------
            h2[PY] = '### OTHER ###\n'
            h2[PY] += 'ifgfilelist: '+join(py_obs_d, ifgfilelist_fr)+'\n'   # Python wants absolute path
            h2[PY] += 'processor: 1\n'
            h2[PY] += 'demfile: '+join(py_dem_d, demfile_fr)+'\n'
            h2[PY] += 'demHeaderFile: '+join(py_dem_d, parfile_fr)+'\n'     # different name for Python
            h2[PY] += '\n'
            cfg.append(h2)

            h3 = ['', '']
            h3[ML] += '### NEEDED FOR MATLAB TO WORK ###\n'
            h3[ML] += 'wavelength: '+wl+'\n'
            h3[ML] += 'gmtfaultfile: xxxx.gmt\n'
            h3[ML] += 'ifgname_prefix_len: 0\n'
            h3[ML] += 'ifgname_date_format: 8\n'
            h3[ML] += 'auxdir: logs'+os.path.sep+'\n'
            # below stuff isn't used but is needed (doesn't matter if directory is created or not)
            h3[ML] += 'ratedir: ratemap'+os.path.sep+'\n'
            h3[ML] += 'tsdir: timeseries'+os.path.sep+'\n'
            h3[ML] += 'profdir: prof'+os.path.sep+'\n'
            h3[ML] += '\n'
            # -----------------------------------------------------------------------------------------
            h3[PY] += '### NEEDED FOR PYTHON TO WORK ###\n'
            h3[PY] += 'pickle: 1\n'       # cmp_out.py needs pickled output     # TODO: change cmp_cfgs.py to read in output files as opposed to pickle so don't need this line
            h3[PY] += '\n'
            cfg.append(h3)

            # stringify the parameters that are recursed in this function
            body = ['', '']
            body[ML] += '### BODY ###\n'
            body[PY] += '### BODY ###\n'
            for item in n_v_s:
                # check if name is ifglks or orbfitlks any other tuples and set appropriately
                if item[NAME] == 'ifglks':
                    for vrsn in [ML, PY]:
                        body[vrsn] += 'ifglksx: '+str(item[VALUE][item[STATE]][0])+'\n'
                        body[vrsn] += 'ifglksy: '+str(item[VALUE][item[STATE]][1])+'\n'
                elif item[NAME] == 'orbfitlks':
                    for vrsn in [ML, PY]:
                        body[vrsn] += 'orbfitlksx: '+str(item[VALUE][item[STATE]][0])+'\n'
                        body[vrsn] += 'orbfitlksy: '+str(item[VALUE][item[STATE]][1])+'\n'
                elif item[NAME] == 'refn':
                    for vrsn in [ML, PY]:
                        body[vrsn] += 'refnx: '+str(item[VALUE][item[STATE]][0])+'\n'
                        body[vrsn] += 'refny: '+str(item[VALUE][item[STATE]][1])+'\n'
                # TODO: consider doing similar with vcmslks
                # other data doesn't need processing
                else:
                    body[ML] += item[NAME]+': '+str(item[VALUE][item[STATE]])+'\n'
                    body[PY] += item[NAME]+': '+str(item[VALUE][item[STATE]])+'\n'
            cfg.append(body)

            # finished. write the two cfg files with diff extensions
            mt_cfg_fp = open(join(num_d, z_str)+'.mcfg', 'w')
            py_cfg_fp = open(join(num_d, z_str)+'.pcfg', 'w')
            for block in cfg:
                mt_cfg_fp.write(block[ML])
                py_cfg_fp.write(block[PY])
            mt_cfg_fp.close()
            py_cfg_fp.close()

            counter += 1
            n_v_s[loop][STATE] += 1     # this must be at bottom of while loop

# =====================================================================================================
# =====================================================================================================
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
            ifgfilelist_fr = fr     # set here because may use 'fr' again
            break
    else:
        print 'ifg file list could not be found. must be in "obs" directory and must begin with "ifms" (e.g. ifms_17)'
        sys.exit(0)

    # get the relative demfile name
    for fr in os.listdir(dem_d):
        if fr[-4:] == '.dem':       # convention
            demfile_fr = fr
            break
    else:
        print '<demfile error message>'         # TODO:
        sys.exit(0)

    # get the relative parfile name
    for fr in os.listdir(dem_d):
        if fr[-4:] == '.par':       # convention
            parfile_fr = fr         # set here because i may use 'fr' again
            break
    else:
        print '<parfile error message>'     # TODO:
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

