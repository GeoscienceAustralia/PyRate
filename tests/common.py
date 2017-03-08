#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains generic utilities and mock objects for use in the
PyRate test suite.
"""

import glob
import os
import shutil
import tempfile
from os.path import join

import numpy as np
from numpy import isnan, sum as nsum
from osgeo import gdal

from pyrate import config as cf, mst, timeseries, matlab_mst, algorithm, \
    ifgconstants as ifc, linrate
from pyrate.shared import Ifg, pre_prepare_ifgs, get_projection_info, \
    write_output_geotiff

TEMPDIR = tempfile.gettempdir()
BASE_TEST = join(os.environ['PYRATEPATH'], "tests", "test_data")
SYD_TEST_DIR = join(BASE_TEST, "sydney_test")
SYD_TEST_OBS = join(SYD_TEST_DIR, 'obs')  # roipac processed unws
SYD_TEST_OUT = join(SYD_TEST_DIR, 'out')
SYD_TEST_TIF = join(SYD_TEST_DIR, 'tif')
SYD_TEST_GAMMA = join(SYD_TEST_DIR, 'gamma_sydney_test')  # gamma processed unws

SYD_TEST_DEM_DIR = join(SYD_TEST_DIR, 'dem')
SYD_TEST_MATLAB_MST_DIR = join(SYD_TEST_DIR, 'matlab_mst')
SYD_TEST_MATLAB_PREPIFG_DIR = join(SYD_TEST_DIR, 'matlab_preifg_output')
SYD_TEST_MATLAB_ORBITAL_DIR = join(SYD_TEST_DIR,
                                   'matlab_orbital_error_correction')
SYD_TEST_DEM_ROIPAC = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.dem')
SYD_TEST_DEM_GAMMA = join(SYD_TEST_GAMMA, '20060619_utm.dem')
SYD_TEST_INCIDENCE = join(SYD_TEST_GAMMA, '20060619_utm.inc')
SYD_TEST_ELEVATION = join(SYD_TEST_GAMMA, '20060619_utm.lv_theta')
SYD_TEST_DEM_HDR_GAMMA = join(SYD_TEST_GAMMA, '20060619_utm_dem.par')
SYD_TEST_DEM_HDR = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.dem.rsc')
SYD_TEST_DEM_TIF = join(SYD_TEST_DEM_DIR, 'sydney_trimmed.tif')
SYDNEY_TEST_CONF = join(SYD_TEST_DIR, 'pyrate_system_test.conf')

PREP_TEST_DIR = join(BASE_TEST, 'prepifg')
PREP_TEST_OBS = join(PREP_TEST_DIR, 'obs')
PREP_TEST_TIF = join(PREP_TEST_DIR, 'tif')

HEADERS_TEST_DIR = join(BASE_TEST, 'headers')
INCID_TEST_DIR = join(BASE_TEST, 'incidence')

GAMMA_TEST_DIR = join(BASE_TEST, "gamma")

# small dummy ifg list to limit overall # of ifgs
IFMS5 = """geo_060828-061211_unw.tif
geo_061106-061211_unw.tif
geo_061106-070115_unw.tif
geo_061106-070326_unw.tif
geo_070326-070917_unw.tif
"""

UNWS5 = """geo_060828-061211.unw
geo_061106-061211.unw
geo_061106-070115.unw
geo_061106-070326.unw
geo_070326-070917.unw
"""

IFMS16 = ['geo_060619-061002_unw.tif',
        'geo_060828-061211_unw.tif',
        'geo_061002-070219_unw.tif',
        'geo_061002-070430_unw.tif',
        'geo_061106-061211_unw.tif',
        'geo_061106-070115_unw.tif',
        'geo_061106-070326_unw.tif',
        'geo_061211-070709_unw.tif',
        'geo_061211-070813_unw.tif',
        'geo_070115-070326_unw.tif',
        'geo_070115-070917_unw.tif',
        'geo_070219-070430_unw.tif',
        'geo_070219-070604_unw.tif',
        'geo_070326-070917_unw.tif',
        'geo_070430-070604_unw.tif',
        'geo_070604-070709_unw.tif']


def sydney_data_setup(datafiles=None, is_dir=False):
    """Returns Ifg objs for the files in the sydney test dir
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if is_dir:
        datafiles = glob.glob(join(datafiles, "*.tif"))
    else:
        if datafiles:
            for i, d in enumerate(datafiles):
                datafiles[i] = os.path.join(SYD_TEST_TIF, d)
        else:
            datafiles = glob.glob(join(SYD_TEST_TIF, "*.tif"))
    datafiles.sort()
    ifgs = [Ifg(i) for i in datafiles]
    
    for i in ifgs: 
        i.open()
        i.nodata_value = 0

    return ifgs


def sydney_ifg_file_list(datafiles=None):
    """Returns the file list of all the .tif files after prepifg conversion
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if datafiles:
        for i, d in enumerate(datafiles):
            datafiles[i] = os.path.join(SYD_TEST_TIF, d)
    else:
        datafiles = glob.glob(join(SYD_TEST_TIF, "*.tif"))
    datafiles.sort()
    return datafiles


def sydney_data_roipac_unws():
    """Returns unw file list before prepifg operation
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    return glob.glob(join(SYD_TEST_OBS, "*.unw"))


def sydney_data_setup_gamma_unws():
    """Returns unw file list before prepifg operation
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    return glob.glob(join(SYD_TEST_GAMMA, "*.unw"))


def sydney5_ifgs():
    """Convenience func to return a subset of 5 linked Ifgs from the testdata"""
    BASE_DIR = tempfile.mkdtemp()
    data_paths = [os.path.join(SYD_TEST_TIF, p) for p in IFMS5.split()]
    new_data_paths = [os.path.join(BASE_DIR, os.path.basename(d))
                      for d in data_paths]
    for d in data_paths:
        shutil.copy(d, os.path.join(BASE_DIR, os.path.basename(d)))

    return [Ifg(p) for p in new_data_paths]


def sydney5_mock_ifgs(xs=3, ys=4):
    '''Returns smaller mocked version of sydney Ifgs for testing'''
    ifgs = sydney5_ifgs()
    for i in ifgs:
        i.open()
        i.nodata_value = 0

    return [MockIfg(i, xs, ys) for i in ifgs]


class MockIfg(object):
    """Mock Ifg for detailed testing"""

    def __init__(self, ifg, xsize=None, ysize=None):
        """
        Creates mock Ifg based on a given interferogram. Size args specify the
        dimensions of the phase band (so the mock ifg can be resized differently
        to the source interferogram for smaller test datasets).
        """
        self.dataset = ifg.dataset
        self.master = ifg.master
        self.slave = ifg.slave
        self.data_path = ifg.data_path
        self.nrows = ysize
        self.ncols = xsize
        self.x_size = ifg.x_size
        self.y_size = ifg.y_size
        self.x_step = ifg.x_step
        self.y_step = ifg.y_step
        self.num_cells = self.ncols * self.nrows
        self.phase_data = ifg.phase_data[:ysize, :xsize]
        self.nan_fraction = ifg.nan_fraction # use existing overall nan fraction
        self.is_open = False

    def __repr__(self, *args, **kwargs):
        return 'MockIfg: %s -> %s' % (self.master, self.slave)

    def open(self):
        # TODO: could move some of the init code here to mimic Ifgs
        pass  # can't actually open anything!

    @property
    def nan_count(self):
        return nsum(isnan(self.phase_data))

    @property
    def shape(self):
        return self.nrows, self.ncols

    def write_modified_phase(self):  #dummy
        pass

    def close(self):  # dummy
        pass


def reconstruct_linrate(shape, tiles, output_dir, out_type):
    rate = np.zeros(shape=shape, dtype=np.float32)
    for t in tiles:
        rate_file = os.path.join(output_dir, out_type +
                                 '_{}.npy'.format(t.index))
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y:t.bottom_right_y,
             t.top_left_x:t.bottom_right_x] = rate_tile
    return rate


def reconstruct_mst(shape, tiles, output_dir):
    mst_file_0 = os.path.join(output_dir, 'mst_mat_{}.npy'.format(0))
    shape0 = np.load(mst_file_0).shape[0]

    mst = np.empty(shape=((shape0,) + shape), dtype=np.float32)
    for i, t in enumerate(tiles):
        mst_file_n = os.path.join(output_dir, 'mst_mat_{}.npy'.format(i))
        mst[:, t.top_left_y:t.bottom_right_y,
                t.top_left_x: t.bottom_right_x] = np.load(mst_file_n)
    return mst


def move_files(source_dir, dest_dir, file_type='*.tif'):
    for filename in glob.glob(os.path.join(source_dir, file_type)):
        shutil.move(filename, dest_dir)


def assert_ifg_phase_equal(ifg_path1, ifg_path2):
    ds1 = gdal.Open(ifg_path1)
    ds2 = gdal.Open(ifg_path2)
    np.testing.assert_array_almost_equal(ds1.ReadAsArray(), ds2.ReadAsArray())
    ds1 = None
    ds2 = None


def prepare_ifgs_without_phase(ifg_paths, params):
    ifgs = [Ifg(p) for p in ifg_paths]
    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        nan_conversion = params[cf.NAN_CONVERSION]
        if nan_conversion:  # nan conversion happens here in networkx mst
            # if not ifg.nan_converted:
            i.nodata_value = params[cf.NO_DATA_VALUE]
            i.convert_to_nans()
    return ifgs


def mst_calculation(ifg_paths_or_instance, params):
    if isinstance(ifg_paths_or_instance, list):
        ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)
        mst_grid = mst.mst_parallel(ifgs, params)
    else:
        nan_conversion = params[cf.NAN_CONVERSION]
        assert isinstance(ifg_paths_or_instance, matlab_mst.IfgListPyRate)
        ifgs = ifg_paths_or_instance.ifgs
        for i in ifgs:
            if not i.mm_converted:
                i.nodata_value = params[cf.NO_DATA_VALUE]
                i.convert_to_mm()
        ifg_instance_updated, epoch_list = \
            matlab_mst.get_nml(ifg_paths_or_instance,
                               nodata_value=params[cf.NO_DATA_VALUE],
                               nan_conversion=nan_conversion)
        mst_grid = matlab_mst.matlab_mst_bool(ifg_instance_updated)

    # write mst output to a file
    mst_mat_binary_file = join(params[cf.OUT_DIR], 'mst_mat')
    np.save(file=mst_mat_binary_file, arr=mst_grid)

    for i in ifgs:
        i.close()
    return mst_grid


def compute_time_series(ifgs, mst_grid, params, vcmt):
    # Calculate time series
    tsincr, tscum, tsvel = calculate_time_series(
        ifgs, params, vcmt=vcmt, mst=mst_grid)

    # tsvel_file = join(params[cf.OUT_DIR], 'tsvel.npy')
    tsincr_file = join(params[cf.OUT_DIR], 'tsincr.npy')
    tscum_file = join(params[cf.OUT_DIR], 'tscum.npy')
    np.save(file=tsincr_file, arr=tsincr)
    np.save(file=tscum_file, arr=tscum)
    # np.save(file=tsvel_file, arr=tsvel)

    # TODO: write tests for these functions
    write_timeseries_geotiff(ifgs, params, tsincr, pr_type='tsincr')
    write_timeseries_geotiff(ifgs, params, tscum, pr_type='tscuml')
    # write_timeseries_geotiff(ifgs, params, tsvel, pr_type='tsvel')
    return tsincr, tscum, tsvel


def calculate_time_series(ifgs, params, vcmt, mst):
    res = timeseries.time_series(ifgs, params, vcmt, mst)
    for r in res:
        if len(r.shape) != 3:
            raise timeseries.TimeSeriesError

    tsincr, tscum, tsvel = res
    return tsincr, tscum, tsvel


# This is not used anywhere now, but may be useful
def time_series_interpolation(ifg_instance_updated, params):

    edges = matlab_mst.get_sub_structure(ifg_instance_updated,
                                         np.zeros(len(ifg_instance_updated.id),
                                                  dtype=bool))

    _, _, ntrees = matlab_mst.matlab_mst_kruskal(edges, ntrees=True)
    # if ntrees=1, no interpolation; otherwise interpolate
    if ntrees > 1:
        params[cf.TIME_SERIES_INTERP] = 1
    else:
        params[cf.TIME_SERIES_INTERP] = 0

    return params


def write_timeseries_geotiff(ifgs, params, tsincr, pr_type):
    # setup metadata for writing into result files
    gt, md, wkt = get_projection_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]

    for i in range(tsincr.shape[2]):
        md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
        md['SEQUENCE_POSITION'] = i+1  # sequence position

        data = tsincr[:, :, i]
        dest = join(params[cf.OUT_DIR], pr_type + "_" +
                    str(epochlist.dates[i + 1]) + ".tif")
        md[ifc.DATA_TYPE] = pr_type
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)


def calculate_linear_rate(ifgs, params, vcmt, mst_mat=None):
    # log.info('Calculating linear rate')
    res = linrate.linear_rate(ifgs, params, vcmt, mst_mat)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    rate, error, samples = res
    write_linrate_tifs(ifgs, params, res)
    # log.info('Linear rate calculated')
    return rate, error, samples


def write_linrate_tifs(ifgs, params, res):
    # log.info('Writing linrate results')
    rate, error, samples = res
    gt, md, wkt = get_projection_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]
    dest = join(params[cf.OUT_DIR], "linrate.tif")
    md[ifc.EPOCH_DATE] = epochlist.dates
    md[ifc.DATA_TYPE] = ifc.LINRATE
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "linerror.tif")
    md[ifc.DATA_TYPE] = ifc.LINERROR
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "linsamples.tif")
    md[ifc.DATA_TYPE] = ifc.LINSAMP
    write_output_geotiff(md, gt, wkt, samples, dest, np.nan)
    write_linrate_numpy_files(error, rate, samples, params)


def write_linrate_numpy_files(error, rate, samples, params):
    rate_file = join(params[cf.OUT_DIR], 'rate.npy')
    error_file = join(params[cf.OUT_DIR], 'error.npy')
    samples_file = join(params[cf.OUT_DIR], 'samples.npy')
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)
