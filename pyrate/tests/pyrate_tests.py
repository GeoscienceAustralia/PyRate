"""
Tests end-to-end runs and functions in PyRate workflow module.

Created on 17/09/2012
@author: Ben Davies, NCI
"""

import os
from os.path import join

import glob
import shutil
import logging
import unittest

from pyrate import pyrate, shared, config, prepifg


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


def test_transform_params():
	params = {config.IFG_LKSX: 3, config.IFG_LKSY: 2, config.IFG_CROP_OPT: 1}
	assert pyrate.transform_params(params) == (3, 2, 1)


def test_warp_required():
	nocrop = prepifg.ALREADY_SAME_SIZE
	assert pyrate.warp_required(xlooks=2, ylooks=1, crop=nocrop)
	assert pyrate.warp_required(xlooks=1, ylooks=2, crop=nocrop)
	assert not pyrate.warp_required(xlooks=1, ylooks=1, crop=nocrop)
	assert not pyrate.warp_required(xlooks=1, ylooks=1, crop=None)

	for c in prepifg.CROP_OPTIONS[:-1]:
		assert pyrate.warp_required(xlooks=1, ylooks=1, crop=c)


def test_src_ifg_paths():
	ifgdir = join(TEST_CORE, 'tif')
	ifglist_path = join(ifgdir, 'ifms_17')
	paths = pyrate.src_ifg_paths(ifglist_path)
	assert paths[0] == join(ifgdir, 'geo_060619-061002.tif'), str(paths[0])
	assert paths[-1] == join(ifgdir, 'geo_070709-070813.tif')


def test_working_ifg_paths():
	# base file list shouldn't change if mlooking not required
	src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
	assert pyrate.working_ifg_paths(src_paths, 1, 1, 4) == src_paths
	assert pyrate.working_ifg_paths(src_paths, 0, 0, 4) == src_paths


def test_working_ifg_paths_mlooked_exists_failure():
	src_paths = ['temp/ifg0.tif', 'temp/ifg1.tif']
	try:
		_ = pyrate.working_ifg_paths(src_paths, 2, 2, 4)
		raise Exception("Didn't fail with missing files")
	except IOError:
		return


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
			pyrate.init_logging(logging.DEBUG)

			params = config._parse_conf_file(WORKFLOW_CONF)
			paths = glob.glob(join(BASE_OUT_DIR, 'geo_*-*.tif'))
			pyrate.process_ifgs(paths, params)
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
		self.assertEqual(len(logpaths), 1)
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


WORKFLOW_CONF = '''#------------------------------------
# input/output parameters
obsdir:       tif/
ifgfilelist:  tif/ifms_17
#atmdir:
#atmfilelist:
#simdir:
demfile:      dem/sydney_trimmed.tif
#eqfile:
#eqdate:
basepflag:    0
# pscorflag - remove postseismic model
pscorflag:    0
outdir:       out/

#------------------------------------
#initial model parameters
#slip rate for the initial forward model in mm/yr
#add more models following nmodels if nmodels>0
#nmodels:        1
#inisliprate_m1: 1
#fmodelfile_m1:

#------------------------------------
#interferogram crop option
#1: minimum 2: maximum 3: customize 4: all ifms already same size
#lksx,lksy: looks in east and north direction respectively
#xfirst,yfirst: x,y of top-left corner
#xlast,ylast: x,y of bottom-right corner
#ifgconv:      1
ifgcropopt:   4
ifglksx:      1
ifglksy:      1
#ifgxfirst:
#ifgxlast:
#ifgyfirst:
#ifgylast:

#------------------------------------
# simulation parameters
# realdata: enter 0 for synthetic/monte-carlo data, 1 for a single set of real data.
# nsets: number of sets of data (multiple sets required for Monte Carlo error estimates
realdata:     1
nsets:        1

#------------------------------------
# stacking parameters
# pthr: minimum number of coherent ifgs for each pixel
# nsig: n-sigma used as residuals threshold for ILSM stacking
# maxsig: maximum residual used as threshold for the stacked interferogram
nsig:          3
pthr:          20
maxsig:        2
#------------------------------------
# orbital errors fitting parameters
# orbfitmethod = 1: interferogram by interferogram; orbfitmethod = 2: epoch by epoch
# degrees: polynomial degrees for the orbital error fitting (1: planar; 2: quadratic)
# orbfitlooksx/y: multi-look processing for orbital correction
# orbrefest: remove reference phase
# orbmaskflag: mask some patches for orbital correction
orbfit:        0
orbfitmethod:  1
orbfitdegrees: 1
orbfitlksx:    0
orbfitlksy:    0
#orbrefest:
#orbmaskflag:

#------------------------------------
# atmospheric delay errors fitting parameters
# atmfitmethod = 1: interferogram by interferogram; atmfitmethod = 2, epoch by epoch
# atmfitlooksx/y: multi-look processing for atm correction
# atmrefest: remove reference phase
#atmfit:        0
#atmfitmethod:  2
#atmfitlksx:    1
#atmfitlksy:    1
#atmrefest:     1
#atmmaskflag:   0

#------------------------------------
#vcm estimation parameters
#vcmtmethod = 1: general method; 2: network estimation for vcm_t;
#vcmsmethod = 1: general method; 2: sparse matrix for the first line; 3: sparse matrix for vcm_s;
#vcmslksx/y: looks also used for slip rate estimation
#vcmtmethod:    1
#vcmsmethod:    2
#vcmslksx:      5
#vcmslksy:      5

#------------------------------------
# reference point options
# refx/y: coordinates of reference points
# refnx/y: number tie points in x/y direction
# refchipsize: chip size of the reference point window
# refminfrac: minimum fraction of coherence pixels
refx:          -1
refy:          -1
#refx:
#refy:
refnx:         5
refny:         5
refchipsize:   5
refminfrac:    0.8

#------------------------------------
#apsest: (1: APS estimation; 0: not)
apsest:         0

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
tsmethod:      1
smorder:       2
smfactor:     -0.25
smf_min:      -3
smf_max:      -1
smf_int:     0.25
lcurv_lksx:    4
lcurv_lksy:    4
ts_pthr:       10
ts_interp:     1

#------------------------------------
# spatially correlated noise low-pass filter parameters
# slpfmethod: filter method (1: butterworth; 2: gaussian)
# slpfcutoff: cutoff d0 for both butterworth and gaussian filters in the same unit with pixel size
# slpforder: order n for butterworth filter (default 1)
slpfmethod:     1
slpfcutoff:     0
slpforder:      1
# temporal high-pass filter parameters
# thpfcutoff: cutoff t0 for gaussian filter in year;
# thpfpthr: valid pixel threshold;
thpfmethod:   1
thpfcutoff:   0.25
thpfpthr:     1

#profile parameters
#provide either pixel numbers or lat / long coordinates. the other will be calculated
#proffmt: 1) lat/long supplied (default) 2) pixel location supplied
#prof=[16,417;316,48];  %unit (pixel number, left profile)
#swath=15;              %unit (km)
#step=1;                %unit (km)
proffmt:      1
#profx0:     165
#profy0:    2920
#profx1:     915
#profy1:       1
profx0:   89.60
profy0:   26.84
profx1:   93.12
profy1:   40.55
profswath:   100
profstep:    10
gmtfaultfile:
#gmtfaultfile:

#------------------------------------
# profile parameters
# swath/step: unit (km)
#make_prof:                1
#profswath:                5
#profstep:                 1
#gmtfaultfile:

#------------------------------------
# iterative algorithm parameters
# tol: tolerance for the convergence
# maxiter: maximam iterative number
tol:          0.2
maxiter:      10
'''
