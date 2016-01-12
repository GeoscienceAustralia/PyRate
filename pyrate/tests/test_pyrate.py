"""
Tests end-to-end runs and functions in PyRate workflow module.

Created on 17/09/2012

.. codeauthor:: Ben Davies
"""

import os, glob, shutil, logging, unittest
from os.path import join
from pyrate import shared, config, prepifg
from pyrate.scripts import run_pyrate
from pyrate import config as cf

# testing constants
BASE_DIR = '/tmp/pyrate/workflow'
BASE_OUT_DIR = join(BASE_DIR, 'out')
BASE_DEM_DIR = join(BASE_DIR, 'dem')
BASE_CFG_FILE = join(BASE_DIR, 'pyrate_workflow_test.conf')
BASE_DEM_FILE = join(BASE_DEM_DIR, 'sydney_trimmed.tif')

TEST_CORE = join(os.environ['PYRATEPATH'], 'tests', 'sydney_test')

CURRENT_DIR = os.getcwd()


# FIXME: make test module to run PyRate with complex options (will be slow)
# TODO: test external ifglist as arg

WORKFLOW_CONF = '''#------------------------------------
# input/output parameters
obsdir:       tif/
ifgfilelist:  tif/ifms_17
#atmdir:
#atmfilelist:
simdir:
demfile:      dem/sydney_trimmed.tif
#eqfile:
#eqdate:
basepflag:    0
# pscorflag - remove postseismic model
pscorflag:    0
outdir:       out/

#------------------------------------
# input/output parameters
# Directory for the unwrapped geocoded interferograms.
obsdir:       tif/

# File containing the list of unwrapped geocoded interferograms.
ifgfilelist:  tif/ifms_17

# InSAR processor ROIPAC = 0, GAMMA = 1
processor:    0

# The DEM used by the interferometric processing software
demfile:      dem/sydney_trimmed.tif

# the projection of the interferograms. This
# is only required if the InSAR processor was
# ROIPAC (specified by 'processor = 1').
#projection

# The resource header used for reading ROIPAC interferograms.
resourceHeader: dem/sydney_trimmed.dem.rsc

# The name of the DEM header file from GAMMA.
# This is only required if the InSAR processor
# was GAMMA (specified by 'processor = 1')
demHeaderFile: dem/sydney_trimmed.dem.par

# Where to write PyRate outputs
outdir:       out/

# No data averaging threshold for prepifg
noDataAveragingThreshold: 0.5

# The no data value in the interferograms
noDataValue:  0.0

# Run NetworkX version of MST (1) or matlab version (0)
networkx_or_matlab:  0

# Nan conversion flag
nan_conversion: 1

#------------------------------------
#interferogram crop option
#1: minimum 2: maximum 3: customize 4: all ifms already same size
#lksx,lksy: looks in east and north direction respectively
#xfirst,yfirst: x,y of top-left corner
#xlast,ylast: x,y of bottom-right corner
ifgcropopt:   4
ifglksx:      1
ifglksy:      1
#ifgxfirst:    150.92
#ifgxlast:     150.94
#ifgyfirst:    -34.18
#ifgylast:     -34.22

#------------------------------------
# reference point options
# refx/y: coordinates of reference points
# refnx/y: number tie points in x/y direction
# refchipsize: chip size of the reference point window
# refminfrac: minimum fraction of coherence pixels
refx:          -1
refy:          -1
refnx:         5
refny:         5
refchipsize:   5
refminfrac:    0.8
refest:        1

#------------------------------------
# orbital errors fitting parameters
# orbfitmethod = 1: interferogram by interferogram; orbfitmethod = 2: epoch by epoch
# degrees: polynomial degrees for the orbital error fitting (1: planar; 2: quadratic; 3: part-cubic)
# orbfitlooksx/y: multi-look processing for orbital correction
# orbrefest: remove reference phase
# orbmaskflag: mask some patches for orbital correction
orbfit:        1
orbfitmethod:  1
orbfitdegrees: 1
orbfitlksx:    0
orbfitlksy:    0
#orbrefest:     1
#orbmaskflag:   0

#------------------------------------
# atmospheric delay errors fitting parameters
# atmfitmethod = 1: interferogram by interferogram; atmfitmethod = 2, epoch by epoch
# atmfitlooksx/y: multi-look processing for atm correction
# atmrefest: remove reference phase
##atmfit:        0 NOT IMPLEMENTED
##atmfitmethod:  2 NOT IMPLEMENTED
##atmfitlksx:    1 NOT IMPLEMENTED
##atmfitlksy:    1 NOT IMPLEMENTED
##atmrefest:     1 NOT IMPLEMENTED
##atmmaskflag:   0 NOT IMPLEMENTED

#------------------------------------
# Spatio-temporal APS filter parameters
# apsest: (1: APS estimation; 0: not)
##apsest:         0 NOT IMPLEMENTED
# spatially correlated noise low-pass filter parameters
# slpfmethod: filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 for both butterworth and gaussian filters in the same unit with pixel size
# slpforder: order n for butterworth filter (default 1)
##slpfmethod:     1 NOT IMPLEMENTED
##slpfcutoff:     0 NOT IMPLEMENTED
##slpforder:      1 NOT IMPLEMENTED
# temporal high-pass filter parameters
# thpfcutoff: cutoff t0 for gaussian filter in year;
# thpfpthr: valid pixel threshold;
##thpfmethod:   1 NOT IMPLEMENTED
##thpfcutoff:   0.25 NOT IMPLEMENTED
##thpfpthr:     1 NOT IMPLEMENTED

#------------------------------------
#vcm estimation parameters
#vcmtmethod = 1: general method; 2: network estimation for vcm_t;
#vcmsmethod = 1: general method; 2: sparse matrix for the first line; 3: sparse matrix for vcm_s;
#vcmslksx/y: looks also used for slip rate estimation
#vcmtmethod:    1
##vcmsmethod:    2 NOT IMPLEMENTED
##vcmslksx:      5 NOT IMPLEMENTED
##vcmslksy:      5 NOT IMPLEMENTED

#------------------------------------
# time series parameters
# tscal: time series calculation (1: ts calculation; 0: no)
# tsmethod: time series method (1: Laplacian smoothing, 2: SVD)
# smorder: order of smoothing operator(1: first-order difference; 2: second-order difference)
# smfactor: smoothing factor(0: calculate & plot L-curve; others: using the specific smoothing factor)
# smf_min/max/int: region of smoothing factors for L-curve calculation, the exact region will be calculated by 10.^(smf)
# lcurve_lksx/lksy: looks number for L-curve calculation
# ts_pthr: pixel number threshold for time series inversion
# ts_interp: interpolate across epoch gaps 0: no (default); 1: yes
tscal:         1
ts_pthr:       10
ts_interp:     1
smfactor:     -0.25
smorder:       2
##tsmethod:      1 ONLY SVD METHOD IMPLEMENTED
##smorder:       2 ONLY SVD METHOD IMPLEMENTED
##smfactor:     -0.25 ONLY SVD METHOD IMPLEMENTED
##smf_min:      -3 ONLY SVD METHOD IMPLEMENTED
##smf_max:      -1 ONLY SVD METHOD IMPLEMENTED
##smf_int:     0.25 ONLY SVD METHOD IMPLEMENTED
##lcurv_lksx:    4 ONLY SVD METHOD IMPLEMENTED
##lcurv_lksy:    4 ONLY SVD METHOD IMPLEMENTED

#------------------------------------
# stacking parameters
# pthr: minimum number of coherent ifgs for each pixel
# nsig: n-sigma used as residuals threshold for ILSM stacking
# maxsig: maximum residual used as threshold for the stacked interferogram
nsig:          3
pthr:          20
maxsig:        2
'''


def test_transform_params():
    params = {config.IFG_LKSX: 3, config.IFG_LKSY: 2, config.IFG_CROP_OPT: 1}
    assert run_pyrate.transform_params(params) == (3, 2, 1)


def test_warp_required():
    nocrop = prepifg.ALREADY_SAME_SIZE
    assert run_pyrate.warp_required(xlooks=2, ylooks=1, crop=nocrop)
    assert run_pyrate.warp_required(xlooks=1, ylooks=2, crop=nocrop)
    assert run_pyrate.warp_required(xlooks=1, ylooks=1, crop=nocrop)
    assert not run_pyrate.warp_required(xlooks=1, ylooks=1, crop=None)

    for c in prepifg.CROP_OPTIONS[:-1]:
        assert run_pyrate.warp_required(xlooks=1, ylooks=1, crop=c)


def test_original_ifg_paths():
    ifgdir = join(TEST_CORE, 'tif')
    ifglist_path = join(ifgdir, 'ifms_17')
    paths = run_pyrate.original_ifg_paths(ifglist_path)
    assert paths[0] == join(ifgdir, 'geo_060619-061002.tif'), str(paths[0])
    assert paths[-1] == join(ifgdir, 'geo_070709-070813.tif')


# def test_working_ifg_paths():
#     # working ifg paths is not used anymore
#     # base file list changes to the correct mlooked file even though
#     # CROPOUT is the same as the original tif
#     src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
#
#     assert run_pyrate.working_ifg_paths(src_paths, 1, 1, 4) == src_paths
#     assert run_pyrate.working_ifg_paths(src_paths, 0, 0, 4) == src_paths


def test_working_ifg_paths_mlooked_exists_failure():
    src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
    try:
        _ = run_pyrate.working_ifg_paths(src_paths, 2, 2, 4)
        raise Exception("Didn't fail with missing files")
    except IOError:
        return


def test_dest_ifg_paths():
    # given source ifgs to process, get paths of ifgs in out dir
    src_paths = ['tif/ifg0.tif', 'tif/ifg1.tif']
    dest_paths = run_pyrate.dest_ifg_paths(src_paths, outdir='out')
    assert dest_paths == ['out/ifg0.tif', 'out/ifg1.tif']


# FIXME: change to read output ifgs
def get_ifgs(_open=True):
    paths = glob.glob(join(BASE_OUT_DIR, 'geo_*-*.tif'))
    ifgs = [shared.Ifg(p) for p in paths]
    assert len(ifgs) == 17, 'Got %s' % ifgs

    if _open:
        for i in ifgs:
            i.open(readonly=False)
    return ifgs


class PyRateTests(unittest.TestCase):
    # Initialise & run workflow from class setup, ignoring multilooking as it is
    # a separate step. Unit tests verify different steps have completed

    @classmethod
    def setUpClass(cls):
        # start each full test run cleanly
        shutil.rmtree(BASE_DIR, ignore_errors=True)

        try:
            # copy source data (treat as prepifg already run)
            os.makedirs(BASE_OUT_DIR)

            for path in glob.glob(join(TEST_CORE, 'tif/*')):
                dest = join(BASE_OUT_DIR, os.path.basename(path))
                shutil.copy(path, dest)
                os.chmod(dest, 0660)

            os.makedirs(BASE_DEM_DIR)
            orig_dem = join(TEST_CORE, 'dem', 'sydney_trimmed.tif')
            os.symlink(orig_dem, BASE_DEM_FILE)
            os.chdir(BASE_DIR)

            # manually start logging as main() is being bypassed
            run_pyrate.init_logging(logging.DEBUG)

            params = config._parse_conf_file(WORKFLOW_CONF)
            params[cf.SIM_DIR] = os.environ['PYRATEPATH']
            paths = glob.glob(join(BASE_OUT_DIR, 'geo_*-*.tif'))
            run_pyrate.process_ifgs(paths, params)
            # TODO: add matlab mst path tests
            os.chdir(CURRENT_DIR)
        except:
            # revert working dir & avoid paths busting other tests
            os.chdir(CURRENT_DIR)
            raise

    def setUp(self):
        # performance: to save constantly opening ifgs
        if not hasattr(self, 'ifgs'):
            self.ifgs = get_ifgs()

    def get_logfile_path(self):
        logpaths = glob.glob(join(BASE_DIR, '*.log'))

        if len(logpaths) != 1:
            msg = 'Log not generated. Use --nologcapture if running nosetests'
            self.fail(msg)
        return logpaths[0]

    def key_check(self, ifg, key, value):
        'Helper to check for metadata flags'
        md = ifg.dataset.GetMetadata()
        self.assertTrue(key in md, 'Missing %s in %s' % (key, ifg.data_path))
        self.assertTrue(md[key], value)

    def test_basic_outputs(self):
        self.assertTrue(os.path.exists(BASE_DIR))
        self.assertTrue(os.path.exists(BASE_OUT_DIR))

        for i in self.ifgs:
            self.assertFalse(i.is_read_only)

        log_path = self.get_logfile_path()
        st = os.stat(log_path)
        self.assertTrue(st.st_size > 0)

    def test_wavelength_conversion(self):
        # ensure phase has been converted from metres to millimetres
        key = 'PHASE_UNITS'
        value = 'MILLIMETRES'

        for i in self.ifgs:
            self.key_check(i, key, value)

    def test_orbital_correction(self):
        key = 'ORBITAL_ERROR'
        value = 'REMOVED'

        for i in self.ifgs:
            self.key_check(i, key, value)

    #def test_mst_matrix(self):
        # TODO: should matrix be written to disk to be read in?
        #raise NotImplementedError

    def test_refpixel_found(self):
        log_path = self.get_logfile_path()
        for line in file(log_path):
            if 'Reference pixel coordinate:' in line:
                return

        self.fail('No reference pixel found')


if __name__ == "__main__":
    unittest.main()
