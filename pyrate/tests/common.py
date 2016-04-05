"""
Collection of generic testing utils and mock objs for PyRate.

.. codeauthor:: Ben Davies, Sudipta Basak
"""

import errno
import os, glob
import shutil
from os.path import join
from numpy import isnan, sum as nsum
import tempfile
from pyrate.shared import Ifg

TEMPDIR = tempfile.gettempdir()
BASE_TEST = join(os.environ['PYRATEPATH'], "tests")
SYD_TEST_DIR = join(BASE_TEST, "sydney_test")
SYD_TEST_OBS = join(SYD_TEST_DIR, 'obs')
SYD_TEST_OUT = join(SYD_TEST_DIR, 'out')
SYD_TEST_TIF = join(SYD_TEST_DIR, 'tif')

SYD_TEST_DEM_DIR = join(SYD_TEST_DIR, 'dem')
SYD_TEST_MATLAB_MST_DIR = join(SYD_TEST_DIR, 'matlab_mst')
SYD_TEST_MATLAB_PREPIFG_DIR = join(SYD_TEST_DIR, 'matlab_preifg_output')
SYD_TEST_MATLAB_ORBITAL_DIR = join(SYD_TEST_DIR,
                                   'matlab_orbital_error_correction')
SYD_TEST_DEM_UNW = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.dem')
SYD_TEST_DEM_HDR = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.dem.rsc')
SYD_TEST_DEM_TIF = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.tif')

PREP_TEST_DIR = join(BASE_TEST, 'prepifg')
PREP_TEST_OBS = join(PREP_TEST_DIR, 'obs')
PREP_TEST_TIF = join(PREP_TEST_DIR, 'tif')

HEADERS_TEST_DIR = join(BASE_TEST, 'headers')
INCID_TEST_DIR = join(BASE_TEST, 'incidence')

GAMMA_TEST_DIR = join(BASE_TEST, "gamma")



# small dummy ifg list to limit overall # of ifgs
IFMS5 = """geo_060828-061211.tif
geo_061106-061211.tif
geo_061106-070115.tif
geo_061106-070326.tif
geo_070326-070917.tif
"""

UNWS5 = """geo_060828-061211.unw
geo_061106-061211.unw
geo_061106-070115.unw
geo_061106-070326.unw
geo_070326-070917.unw
"""

IFMS16 = ['geo_060619-061002.tif',
        'geo_060828-061211.tif',
        'geo_061002-070219.tif',
        'geo_061002-070430.tif',
        'geo_061106-061211.tif',
        'geo_061106-070115.tif',
        'geo_061106-070326.tif',
        'geo_061211-070709.tif',
        'geo_061211-070813.tif',
        'geo_070115-070326.tif',
        'geo_070115-070917.tif',
        'geo_070219-070430.tif',
        'geo_070219-070604.tif',
        'geo_070326-070917.tif',
        'geo_070430-070604.tif',
        'geo_070604-070709.tif']


def sydney_data_setup(datafiles=None):
    """Returns Ifg objs for the files in the sydney test dir
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if datafiles:
        for i, d in enumerate(datafiles):
            datafiles[i] = os.path.join(SYD_TEST_TIF, d)
    else:
        datafiles = glob.glob(join(SYD_TEST_TIF, "*.tif"))
    datafiles.sort()
    ifgs = [Ifg(i) for i in datafiles]
    
    for i in ifgs: 
        i.open()

    return ifgs


def sydney_data_setup_ifg_file_list(datafiles=None):
    """Returns the file list of all the .tif files after prepifg conversion
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if datafiles:
        for i, d in enumerate(datafiles):
            datafiles[i] = os.path.join(SYD_TEST_TIF, d)
    else:
        datafiles = glob.glob(join(SYD_TEST_TIF, "*.tif"))
    return datafiles


def sydney_data_setup_unw_file_list():
    """Returns unw file list before prepifg operation
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    return glob.glob(join(SYD_TEST_OBS, "*.unw"))


def sydney5_ifgs():
    """Convenience func to return a subset of 5 linked Ifgs from the testdata"""
    return [Ifg(join(SYD_TEST_TIF, p)) for p in IFMS5.split()]


def sydney5_mock_ifgs(xs=3, ys=4):
    '''Returns smaller mocked version of sydney Ifgs for testing'''
    ifgs = sydney5_ifgs()
    for i in ifgs:
        i.open()

    return [MockIfg(i, xs, ys) for i in ifgs]


class MockIfg(object):
    """Mock Ifg for detailed testing"""

    def __init__(self, ifg, xsize=None, ysize=None):
        """
        Creates mock Ifg based on a given interferogram. Size args specify the
        dimensions of the phase band (so the mock ifg can be resized differently
        to the source interferogram for smaller test datasets).
        """
        self.master = ifg.master
        self.slave = ifg.slave

        self.nrows = ysize
        self.ncols = xsize
        self.x_size = ifg.x_size
        self.y_size = ifg.y_size
        self.x_step = ifg.x_step
        self.y_step = ifg.y_step
        self.num_cells = self.ncols * self.nrows
        self.phase_data = ifg.phase_data[:ysize, :xsize]
        self.nan_fraction = ifg.nan_fraction # use existing overall nan fraction

    def __repr__(self, *args, **kwargs):
        return 'MockIfg: %s -> %s' % (self.master, self.slave)

    def open(self):
        # TODO: could move some of the init code here to mimic Ifgs
        pass # can't actually open anything!

    @property
    def nan_count(self):
        return nsum(isnan(self.phase_data))

    @property
    def shape(self):
        return (self.nrows, self.ncols)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def move_files(source_dir, dest_dir, file_type='*.tif'):
    for filename in glob.glob(os.path.join(source_dir, file_type)):
        shutil.move(filename, dest_dir)