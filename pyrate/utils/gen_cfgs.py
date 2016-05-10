import sys
import os
if not ('PR_TESTING_DATA' in os.environ):
    print 'need to set your PR_TESTING_DATA environment variable'
    sys.exit(0)
'''
# may not want this. used pwd because will have the file and then the directory
if not ('PR_TESTING_OUT_DIR' in os.environ):
    print 'need to set your PR_TESTING_DATA environment variable'
    sys.exit(0)
'''
# know where the data and out directory is. can proceed

from numpy import arange
from subprocess import call
from copy import deepcopy

from os.path import join        # better than above method

PR, ML = 0, 1                   # to indicate PyRate or Matlab
NAME, VALUE, STATE = 0, 1, 2    # used in recursive function to generate cfg files

# PARAM SWEEPS ============================================
# =========================================================
base = [
    ['networkx_or_matlab', [0]],
    ['luigi', [0]],

    ['atmfit', [0]],     # Matt says don't use this
                         # normally would have to specify atmfitmethod, atmfitlksx,y afterward

    ['basepflag', [0]],  # PYTHON = NOT CURRENTLY USED, matlab is being used
    ['unwflag', [0]],    # closure check and mask unwrapping errors flag added on 17/09/2012 (MATLAB ONLY. 1 DEFAULT)
                         # not implemented in python. turn it off.
    ['prjflag', [3]],    # Re-project data from Line of sight, 1 = vertical, 2 = horizontal, 3 = no conversion
                         # 3 or run into problems. try others later

    #['noDataAveragingThreshold', arange(0, 1.0, 0.50).tolist()],
    ['noDataAveragingThreshold', [0.5]],    # Sudipta says this is for network... doesn't look like it's even used

    ['noDataValue', [0.0]],     # confusing... seems like no value...
    ['nan_conversion', [1]],    # if you have a noDataValue convert it to NaN?

    ['ifglksx', [1]],       # don't rescale...
    ['ifglksy', [1]],

    # stacking/velocity parameters
    ['nsig', [2]],      # matlab default 2
    ['pthr', [12]],     # dependant on ifgsize
    ['maxsig', [2]],    # <-- matlab default 2

    # Matt says vcm isn't being used, but would be base sweep if was
    # needed for windows version
    ['vcmtmethod', [1]],    # pyrate uses general method1... keep it @ this
    ['vcmsmethod', [1]],    # is calculated by not even used, only one that has looks, defaults to 10 looks
    ['vcmslksx', [10]],      # what to set for matlab
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
    #    ['orbfitmethod', [1,2]],
    #    ['orbfitdegrees', [1]],    # default 1
    #    ['orbfitlksx', [5]],       # will probably depend on data
    #    ['orbfitlksy', [5]],
    #],
]

# Todo, how to pick manual reference point
# does some data have a optimized, predetermined reference point
ref_b = [
    [
        ['refest', [2]],    # 1 is stupid
        ['refx', [0]],
        ['refy', [0]],
        ['refnx', [10]],        # x/10
        ['refny', [10]],        # y/10
        ['refchipsize', [5]],   # ?
        ['refminfrac', [0.8]],      # default for matlab
    ],
]

ts_b = [
    #[
    #    ['tscal', [0]],        # will never test this
    #],
    [   # laplcian smoothing
        ['tscal', [1]],
        ['tsmethod', [1]],

        ['smorder', [1]],     # pyrate config.py says this and below "NOT CURRENTLY USED"
        ['smfactor', [1]],      # MATLAB = smfactor: smoothing factor(0: calculate & plot L-curve; others: smoothing factor (10^(smf))
                                # PYRATE = Laplacian smoothing factor (0: calculate & plot L-curve; others: using the specific smoothing factor 10**smfactor) NOT CURRENTLY USED
                                # 0 smfactor screws up matlab if tspar.smf==0
                                #                               tspar.smf_min = getpar('smf_min:',parmat,'n');
        ['ts_pthr', [12]],      # based on ifg size

        ['ts_interp', [0]],     # MATLAB LAPLACIAN only. default 0. Matt says keep this @ 0

        ['stepest', [0]],       # matlab only, default 0
    ],
    [
        # SVD (note pyrate required smorder & smfactor parameters (to run), even though they won't be used ==> dummy param
        ['tscal', [1]],
        ['tsmethod', [2]],

        ['smorder', [0]],
        ['smfactor', [0]],

        ['ts_pthr', [12]],      # based on ifg size

        ['stepest', [0]],       # matlab only, default 0
    ],
]

# combine branches to create work
work = []
for i1 in crop_b:
    for i2 in orb_b:
        for i3 in ref_b:
            for i4 in ts_b:
                work.append(deepcopy(base + i1 + i2 + i3 + i4))

# set state to 0 while we here...
for b_comb in work:     # b_com = branch combination, which is what it is
    [n_v.append(0) for n_v in b_comb]

# figure out how many permutations (and hence how many 0's we need)
n_perms = 0
for b_comb in work:
    n_perms_per_b = 1
    for n_v_s in b_comb:
        n_perms_per_b *= len(n_v_s[1])
    n_perms += n_perms_per_b

print n_perms
n_zs = len(str(n_perms))

# FUNCTIONS ===============================================
# =========================================================
def rec_loop(cnfg_headers, n_v_s, loop=0):
    """

    :param cnfg_headers: list containing header strings for unix then windows (i.e. just 2 elements)
    :param n_v_s: [param name string, [values that name takes], what value we are up to in the recursive loop]
    :param loop:
    :return: nothing
    """
    if loop != 0:
        if n_v_s[loop][STATE] == len(n_v_s[loop][VALUE]):
            n_v_s[loop][STATE] = 0
    if loop < (len(n_v_s)-1):
        while n_v_s[loop][STATE] < len(n_v_s[loop][VALUE]):
            rec_loop(cnfg_headers, n_v_s, loop+1)
            n_v_s[loop][STATE] += 1
    if loop == (len(n_v_s)-1):
        while (n_v_s[loop][STATE] < len(n_v_s[loop][VALUE])):
            global counter, out_d, n_zs     # out_d and n_zs aren't modified, so technically can leave out here. these are set below, not above.

            z_str = str(counter).zfill(n_zs)

            this_d = join(out_d, z_str)

            # check if the following directories exist. if not create them
            directs = []
            directs.append(this_d)
            directs.append(join(this_d, 'python'))
            directs.append(join(this_d, 'matlab'))
            directs.append(join(this_d, 'matlab', 'logs'))
            for direct in directs:
                if not os.path.exists(direct):
                    os.mkdir(direct)

            cnfg_str = ['', '']
            for item in n_v_s:
                cnfg_str[PR] += item[NAME] + ': ' + str(item[VALUE][item[STATE]]) + '\n'
                cnfg_str[ML] += item[NAME] + ': ' + str(item[VALUE][item[STATE]]) + '\n'

            header_str = ['', '']
            for i, opsys in enumerate(cnfg_headers):
                for n_v in opsys:
                    header_str[i] += n_v[NAME] + ': ' + n_v[VALUE] + '\n'

            outdir_str = [None, None]
            outdir_str[PR] = 'outdir: '+this_d+os.path.sep+'python'+os.path.sep+'\n'
            outdir_str[ML] = 'outdir: '+this_d+os.path.sep+'matlab'+os.path.sep+'\n'

            cnfg_strs = [None, None]

            cnfg_strs[PR] = header_str[PR]+outdir_str[PR]+cnfg_str[PR]

            cnfg_strs[ML] =  '############# needed to work #############\n'
            cnfg_strs[ML] += 'gmtfaultfile: xxxx.gmt\n'
            #cnfg_strs[WIN] += 'atmfit: 0\n'
            cnfg_strs[ML] += 'ifgname_prefix_len: 0\n'
            cnfg_strs[ML] += 'ifgname_date_format: 8\n'
            cnfg_strs[ML] += 'auxdir: logs'+os.path.sep+'\n'
            # below stuff isn't used but is needed (doesn't matter if directory is created or not)
            cnfg_strs[ML] += 'ratedir: ratemap'+os.path.sep+'\n'
            cnfg_strs[ML] += 'tsdir: timeseries'+os.path.sep+'\n'
            cnfg_strs[ML] += 'profdir: prof'+os.path.sep+'\n'
            cnfg_strs[ML] += '##########################################\n'
            # add the important parts onto the "need to work" shit above
            cnfg_strs[ML] += header_str[ML]+outdir_str[ML]+cnfg_str[ML]

            # finished. write the two config files with diff extensions
            cnfg_file_PR = open(this_d + os.path.sep + z_str + '.pcfg', 'w')
            cnfg_file_PR.write(cnfg_strs[PR])
            cnfg_file_ML = open(this_d + os.path.sep + z_str + '.mcfg', 'w')
            cnfg_file_ML.write(cnfg_strs[ML])

            counter += 1
            n_v_s[loop][STATE] += 1 # <-- this must be @ the bottom of while

# SETUP FOLDERS
# -------------

dat_d = os.environ['PR_TESTING_DATA']   # directory where ifg data stored. see file docstring for required structure
#dat_flds = os.listdir(dat_d)           # <-- want to use everything in here
dat_flds = ['syd_g']                    # manually select specific folders

# create a folder of the same name as this script in the same directory as the script
test_d, script_fn = os.path.split(os.path.realpath(__file__))    # will give the directory script is in
script_d = join(test_d, script_fn[:-3])
if not os.path.exists(script_d):
        os.mkdir(script_d)

for fld in dat_flds:
    # cnfg file info to do with directories and stuff that won't get sweeped in rec_loop function
    # list of [name, value]
    # these work automatically for the Python version. a object will be created for the Matlab cnfg file
    header_n_v = []

    fld_abs = join(dat_d, fld)      #cp_app(dat_d, fld)
    obsdir =  join(fld_abs, 'obs')
    header_n_v.append(['obsdir', obsdir+os.path.sep])   # seperator on end to indicated directory (also might be required for matlab)

    for file_name in os.listdir(obsdir):
        if file_name[:4] == 'ifms':
            ifgfilelist = join(obsdir, file_name)       #cp_app(obsdir, file_name)
            break
    header_n_v.append(['ifgfilelist', ifgfilelist])     # this has to be changed to relative path for Matlab

    dem_dir = join(fld_abs, 'dem')                      #cp_app(fld_abs, 'dem')
    for file_name in os.listdir(dem_dir):
        if file_name[-4:] == '.dem':
            demfile = join(dem_dir, file_name)          #cp_app(dem_dir, file_name)
            break
    header_n_v.append(['demfile', demfile])             # relative for Matlab

    # figure out if roipac or gamma pre-processing by file extensions
    # TEST FILES WILL ONLY BE GAMMA FOR NOW
    for file_name in os.listdir(obsdir):
        if file_name[-4:] == '.rsc':
            processor = '0'   # is ROIPAC
            break
        if file_name[-4:] == '.par':
            processor = '1'   # is gamma
            break
    header_n_v.append(['processor', processor])

    if processor == '1':
        demHeaderFile = deepcopy(demfile)
        demHeaderFile += '.par'
        resourceHeader = ''
    else:
        demHeaderFile = ''
        resourceHeader = deepcopy(demfile)
        resourceHeader += '.rsc'
    header_n_v.append(['demHeaderFile', demHeaderFile])
    header_n_v.append(['resourceHeader', resourceHeader])

    # make folder for this specific test data...
    out_d = join(script_d, fld)
    if not os.path.exists(out_d):
        os.mkdir(out_d)
        
    counter = 0     # what cfg file out of all cfg files
    
    headers_n_v = [None, None]
    headers_n_v[PR] = deepcopy(header_n_v)     # pyrate version shouldn't have to change
    temp = []
    for n_v in header_n_v:
        if n_v[NAME] == 'ifgfilelist':
            n_v[VALUE] = os.path.split(n_v[VALUE])[1]
        if n_v[NAME] == 'demfile':
            n_v[VALUE] = os.path.split(n_v[VALUE])[1]
        if n_v[NAME] == 'demHeaderFile':
            n_v[NAME] = 'parfile'
            n_v[VALUE] = os.path.split(n_v[VALUE])[1]
        if n_v[NAME] == 'resourceHeader':
            continue    # we don't want this
        if n_v[NAME] == 'processor':
            continue
        temp.append(n_v)
    temp.append(['demdir', dem_dir+os.path.sep])       # stupid matlab requires dem_dir
    headers_n_v[ML] = temp

    # TODO: might need to add slashes onto end of directories...

    # add wavelength for WIN and MCC
    # ensure there is a file that specifies wavelength for each data folder... (needed for matlab)
    wl_f = open(join(fld_abs, 'wavelength.txt'), 'r')
    wl = wl_f.readline()
    wl_f.close()
    headers_n_v[ML].append(['wavelength', wl])

    for b_comb in work:
        rec_loop(cnfg_headers=headers_n_v, n_v_s=b_comb)
        # set state back to 0 for next folder
        for n_v_s in b_comb:
            n_v_s[STATE] = 0
