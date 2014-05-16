'''
Collection of generic testing utils and mock objs for PyRate
Author: Ben Davies
'''

import glob
from os.path import join

from shared import Ifg
from numpy import ndarray, float32, isnan, sum as nsum


SYD_TEST_DIR = "../tests/sydney_test"
SYD_TEST_OBS = join(SYD_TEST_DIR, 'obs')
SYD_TEST_DEM = join(SYD_TEST_DIR, 'dem')

PREP_TEST_DIR = '../tests/prepifg'
PREP_TEST_OBS = join(PREP_TEST_DIR, 'obs')

SINGLE_TEST_DIR = '../tests/single'


# small dummy ifg list to limit overall # of ifgs
IFMS5 = """geo_060828-061211.unw
geo_061106-061211.unw
geo_061106-070115.unw
geo_061106-070326.unw
geo_070326-070917.unw
"""

# TODO: get rid of first returned arg?
def sydney_test_setup():
	'''Returns Ifg objs for the files in the sydney test dir'''
	datafiles = glob.glob(join(SYD_TEST_OBS, "*.unw") )
	ifgs = [Ifg(i) for i in datafiles]
	for i in ifgs:
		i.open()

	return SYD_TEST_OBS, ifgs


def sydney5_ifgs():
	'''Convenience func to return a subset of 5 linked Ifgs from the testdata'''
	return [Ifg(join(SYD_TEST_OBS, p)) for p in IFMS5.split()]


def sydney5_mock_ifgs(xs=3, ys=4):
	'''Returns smaller mocked version of sydney Ifgs for testing'''
	ifgs = sydney5_ifgs()
	for i in ifgs: i.open()
	mocks = [MockIfg(i, xs, ys) for i in ifgs]
	for i,m in zip(ifgs, mocks):
		m.phase_data = i.phase_data[:ys,:xs]

	return mocks



class MockIfg(object):
	'''Mock Ifg for detailed testing'''

	def __init__(self, ifg, xsize=None, ysize=None):
		'''Creates mock ifg based on another interferogram. Size args specify the
		dimensions of the phase band (so the mock ifg can be resized differently to
		the source interferogram for smaller test datasets).
		'''
		self.MASTER = ifg.MASTER
		self.SLAVE = ifg.SLAVE
		self.DATE12 = ifg.DATE12

		self.FILE_LENGTH = ysize
		self.WIDTH = xsize
		self.X_SIZE = ifg.X_SIZE
		self.Y_SIZE = ifg.Y_SIZE
		self.X_STEP = ifg.X_STEP
		self.Y_STEP = ifg.Y_STEP
		self.num_cells = ysize * xsize
		self.phase_data = ndarray((self.FILE_LENGTH, self.WIDTH), dtype=float32)
		self.nan_fraction = ifg.nan_fraction # use existing overall nan fraction

	def open(self):
		pass # can't open anything!

	@property
	def nan_count(self):
		return nsum(isnan(self.phase_data))

	@property
	def shape(self):
		return (self.FILE_LENGTH, self.WIDTH)
