#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
import stat
import tempfile
import logging
from os.path import join

import numpy as np
from numpy import isnan, sum as nsum
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from osgeo import gdalconst
from osgeo import gdal_array

import sys

PYRATE_MODULE_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pyrate")
TEST_MODULE_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(PYRATE_MODULE_PATH)
sys.path.append(TEST_MODULE_PATH)

from core import algorithm, ifgconstants as ifc, config as cf, timeseries, mst, stack
from core.shared import Ifg, nan_and_mm_convert, get_geotiff_header_info, write_output_geotiff
from core.logger import pyratelogger as log

PYRATEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPDIR = tempfile.gettempdir()
BASE_TEST = join(PYRATEPATH, "tests", "test_data")
SML_TEST_DIR = join(BASE_TEST, "small_test")
SML_TEST_OBS = join(SML_TEST_DIR, "roipac_obs")  # roipac processed unws
SML_TEST_OUT = join(SML_TEST_DIR, "out")
SML_TEST_TIF = join(SML_TEST_DIR, "tif")
SML_TEST_GAMMA = join(SML_TEST_DIR, "gamma_obs")  # gamma processed unws
SML_TEST_CONF = join(SML_TEST_DIR, "conf")
SML_TEST_GAMMA_HEADER_LIST = join(SML_TEST_GAMMA, "headers")

SML_TEST_DEM_DIR = join(SML_TEST_DIR, "dem")
SML_TEST_LEGACY_PREPIFG_DIR = join(SML_TEST_DIR, "prepifg_output")
SML_TEST_LEGACY_ORBITAL_DIR = join(SML_TEST_DIR, "orbital_error_correction")
SML_TEST_DEM_ROIPAC = join(SML_TEST_DEM_DIR, "roipac_test_trimmed.dem")
SML_TEST_DEM_GAMMA = join(SML_TEST_GAMMA, "20060619_utm.dem")
SML_TEST_INCIDENCE = join(SML_TEST_GAMMA, "20060619_utm.inc")
SML_TEST_ELEVATION = join(SML_TEST_GAMMA, "20060619_utm.lv_theta")
SML_TEST_DEM_HDR_GAMMA = join(SML_TEST_GAMMA, "20060619_utm_dem.par")
SML_TEST_DEM_HDR = join(SML_TEST_DEM_DIR, "roipac_test_trimmed.dem.rsc")
SML_TEST_DEM_TIF = join(SML_TEST_DEM_DIR, "roipac_test_trimmed.tif")

SML_TEST_COH_DIR = join(SML_TEST_DIR, "coherence")
SML_TEST_COH_LIST = join(SML_TEST_COH_DIR, "coherence_17")

TEST_CONF_ROIPAC = join(SML_TEST_CONF, "pyrate_roipac_test.conf")
TEST_CONF_GAMMA = join(SML_TEST_CONF, "pyrate_gamma_test.conf")

PREP_TEST_DIR = join(BASE_TEST, "prepifg")
PREP_TEST_OBS = join(PREP_TEST_DIR, "obs")
PREP_TEST_TIF = join(PREP_TEST_DIR, "tif")

HEADERS_TEST_DIR = join(BASE_TEST, "headers")
INCID_TEST_DIR = join(BASE_TEST, "incidence")

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

IFMS16 = [
    "geo_060619-061002_unw.tif",
    "geo_060828-061211_unw.tif",
    "geo_061002-070219_unw.tif",
    "geo_061002-070430_unw.tif",
    "geo_061106-061211_unw.tif",
    "geo_061106-070115_unw.tif",
    "geo_061106-070326_unw.tif",
    "geo_061211-070709_unw.tif",
    "geo_061211-070813_unw.tif",
    "geo_070115-070326_unw.tif",
    "geo_070115-070917_unw.tif",
    "geo_070219-070430_unw.tif",
    "geo_070219-070604_unw.tif",
    "geo_070326-070917_unw.tif",
    "geo_070430-070604_unw.tif",
    "geo_070604-070709_unw.tif",
]


def remove_tifs(path):
    """

    Args:
      path: 

    Returns:

    """
    tifs = glob.glob(os.path.join(path, "*.tif"))
    for tif in tifs:
        os.remove(tif)


def small_data_setup(datafiles=None, is_dir=False):
    """Returns Ifg objs for the files in the small test dir input phase data is
    in radians; these ifgs are in radians - not converted to mm

    Args:
      datafiles: Default value = None)
      is_dir: Default value = False)

    Returns:

    """
    if is_dir:
        datafiles = glob.glob(join(datafiles, "*.tif"))
    else:
        if datafiles:
            for i, d in enumerate(datafiles):
                datafiles[i] = os.path.join(SML_TEST_TIF, d)
        else:
            datafiles = glob.glob(join(SML_TEST_TIF, "*.tif"))
    datafiles.sort()
    ifgs = [Ifg(i) for i in datafiles]

    for i in ifgs:
        i.open()
        i.nodata_value = 0

    return ifgs


def small_ifg_file_list(datafiles=None):
    """Returns the file list of all the .tif files after prepifg conversion
    input phase data is in radians; these ifgs are in radians - not converted to
    mm

    Args:
      datafiles: Default value = None)

    Returns:

    """
    if datafiles:
        for i, d in enumerate(datafiles):
            datafiles[i] = os.path.join(SML_TEST_TIF, d)
    else:
        datafiles = glob.glob(join(SML_TEST_TIF, "*.tif"))
    datafiles.sort()
    return datafiles


def small_data_roipac_unws():
    """Returns unw file list before prepifg operation input phase data is in
    radians; these ifgs are in radians - not converted to mm

    Args:

    Returns:

    """
    return glob.glob(join(SML_TEST_OBS, "*.unw"))


def small_data_setup_gamma_unws():
    """Returns unw file list before prepifg operation input phase data is in
    radians; these ifgs are in radians - not converted to mm

    Args:

    Returns:

    """
    return glob.glob(join(SML_TEST_GAMMA, "*.unw"))


def small5_ifgs():
    """Convenience func to return a subset of 5 linked Ifgs from the testdata"""
    BASE_DIR = tempfile.mkdtemp()
    data_paths = [os.path.join(SML_TEST_TIF, p) for p in IFMS5.split()]
    new_data_paths = [os.path.join(BASE_DIR, os.path.basename(d)) for d in data_paths]
    for d in data_paths:
        shutil.copy(d, os.path.join(BASE_DIR, os.path.basename(d)))

    return [Ifg(p) for p in new_data_paths]


def small5_mock_ifgs(xs=3, ys=4):
    """Returns smaller mocked version of small Ifgs for testing

    Args:
      xs: Default value = 3)
      ys: Default value = 4)

    Returns:

    """
    ifgs = small5_ifgs()
    for i in ifgs:
        i.open()
        i.nodata_value = 0

    return [MockIfg(i, xs, ys) for i in ifgs]


class MockIfg(object):
    """Mock Ifg for detailed testing"""

    def __init__(self, ifg, xsize=None, ysize=None):
        """Creates mock Ifg based on a given interferogram. Size args specify
        the dimensions of the phase band (so the mock ifg can be resized
        differently to the source interferogram for smaller test datasets).

        Args:
            ifg:
            xsize:
            ysize:
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
        self.nan_fraction = ifg.nan_fraction  # use existing overall nan fraction
        self.is_open = False

    def __repr__(self, *args, **kwargs):
        """
        Args:
            *args:
            **kwargs:
        """
        return "MockIfg: %s -> %s" % (self.master, self.slave)

    def open(self):
        """ """
        # TODO: could move some of the init code here to mimic Ifgs
        pass  # can't actually open anything!

    @property
    def nan_count(self):
        """ """
        return nsum(isnan(self.phase_data))

    @property
    def shape(self):
        """ """
        return self.nrows, self.ncols

    def write_modified_phase(self):  # dummy
        """ """
        pass

    def close(self):  # dummy
        """ """
        pass


def reconstruct_stack_rate(shape, tiles, output_dir, out_type):
    """

    Args:
      shape: param tiles:
      output_dir: param out_type:
      tiles: 
      out_type: 

    Returns:

    """
    rate = np.zeros(shape=shape, dtype=np.float32)
    for t in tiles:
        rate_file = os.path.join(output_dir, out_type + "_{}.npy".format(t.index))
        rate_tile = np.load(file=rate_file)
        rate[t.top_left_y : t.bottom_right_y, t.top_left_x : t.bottom_right_x] = rate_tile
    return rate


def reconstruct_mst(shape, tiles, output_dir):
    """

    Args:
      shape: param tiles:
      output_dir: 
      tiles: 

    Returns:

    """
    mst_file_0 = os.path.join(output_dir, "mst_mat_{}.npy".format(0))
    shape0 = np.load(mst_file_0).shape[0]

    mst = np.empty(shape=((shape0,) + shape), dtype=np.float32)
    for i, t in enumerate(tiles):
        mst_file_n = os.path.join(output_dir, "mst_mat_{}.npy".format(i))
        mst[:, t.top_left_y : t.bottom_right_y, t.top_left_x : t.bottom_right_x] = np.load(mst_file_n)
    return mst


def move_files(source_dir, dest_dir, file_type="*.tif"):
    """

    Args:
      source_dir: param dest_dir:
      file_type: Default value = "*.tif")
      dest_dir: 

    Returns:

    """
    for filename in glob.glob(os.path.join(source_dir, file_type)):
        shutil.move(filename, dest_dir)


def assert_ifg_phase_equal(ifg_path1, ifg_path2):
    """

    Args:
      ifg_path1: param ifg_path2:
      ifg_path2: 

    Returns:

    """
    ds1 = gdal.Open(ifg_path1)
    ds2 = gdal.Open(ifg_path2)
    np.testing.assert_array_almost_equal(ds1.ReadAsArray(), ds2.ReadAsArray())
    ds1 = None
    ds2 = None


def prepare_ifgs_without_phase(ifg_paths, params):
    """Args:
        ifg_paths:

    Args:
      ifg_paths: param params:
      params: 

    Returns:

    """
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
    """Args:
        ifg_paths_or_instance:

    Args:
      ifg_paths_or_instance: param params:
      params: 

    Returns:

    """
    if isinstance(ifg_paths_or_instance, list):
        ifgs = pre_prepare_ifgs(ifg_paths_or_instance, params)
        mst_grid = mst.mst_parallel(ifgs, params)
        # write mst output to a file
        mst_mat_binary_file = join(params[cf.OUT_DIR], "mst_mat")
        np.save(file=mst_mat_binary_file, arr=mst_grid)

        for i in ifgs:
            i.close()
        return mst_grid
    return None


def get_nml(ifg_list_instance, nodata_value, nan_conversion=False):
    """

    Args:
      ifg_list_instance(xxx(eg str): xxxx
      nodata_value(float): No data value in image
      nan_conversion(bool, optional): Convert NaNs (Default value = False)

    Returns:
      list: ifg_list_instance: replaces in place

    """
    _epoch_list, n = algorithm.get_epochs(ifg_list_instance.ifgs)
    ifg_list_instance.reshape_n(n)
    if nan_conversion:
        ifg_list_instance.update_nan_frac(nodata_value)
        # turn on for nan conversion
        ifg_list_instance.convert_nans(nan_conversion=nan_conversion)
    ifg_list_instance.make_data_stack()
    return ifg_list_instance, _epoch_list


def compute_time_series(ifgs, mst_grid, params, vcmt):
    """

    Args:
      ifgs: param mst_grid:
      params: param vcmt:
      mst_grid: 
      vcmt: 

    Returns:

    """
    # Calculate time series
    """
    Args:
        ifgs:
        mst_grid:
        params:
        vcmt:
    """
    tsincr, tscum, tsvel = calculate_time_series(ifgs, params, vcmt=vcmt, mst=mst_grid)

    # tsvel_file = join(params[cf.OUT_DIR], 'tsvel.npy')
    tsincr_file = join(params[cf.OUT_DIR], "tsincr.npy")
    tscum_file = join(params[cf.OUT_DIR], "tscum.npy")
    np.save(file=tsincr_file, arr=tsincr)
    np.save(file=tscum_file, arr=tscum)
    # np.save(file=tsvel_file, arr=tsvel)

    # TODO: write tests for these functions
    write_timeseries_geotiff(ifgs, params, tsincr, pr_type="tsincr")
    write_timeseries_geotiff(ifgs, params, tscum, pr_type="tscuml")
    # write_timeseries_geotiff(ifgs, params, tsvel, pr_type='tsvel')
    return tsincr, tscum, tsvel


def calculate_time_series(ifgs, params, vcmt, mst):
    """Args:
        ifgs:

    Args:
      vcmt: param mst:
      ifgs: param params:
      params: 
      mst: 

    Returns:

    """
    res = timeseries.time_series(ifgs, params, vcmt, mst)
    for r in res:
        if len(r.shape) != 3:
            raise timeseries.TimeSeriesError

    tsincr, tscum, tsvel = res
    return tsincr, tscum, tsvel


def write_timeseries_geotiff(ifgs, params, tsincr, pr_type):
    """

    Args:
      ifgs: param params:
      tsincr: param pr_type:
      params: 
      pr_type: 

    Returns:

    """
    # setup metadata for writing into result files
    """
    Args:
        ifgs:
        params:
        tsincr:
        pr_type:
    """
    gt, md, wkt = get_geotiff_header_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]

    for i in range(tsincr.shape[2]):
        md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
        md["SEQUENCE_POSITION"] = i + 1  # sequence position

        data = tsincr[:, :, i]
        dest = join(params[cf.OUT_DIR], pr_type + "_" + str(epochlist.dates[i + 1]) + ".tif")
        md[ifc.DATA_TYPE] = pr_type
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)


def calculate_stacked_rate(ifgs, params, vcmt, mst_mat=None):
    """

    Args:
      ifgs: param params:
      vcmt: param mst_mat:  (Default value = None)
      params: 
      mst_mat:  (Default value = None)

    Returns:

    """
    # log.info('Calculating stacked rate')
    """
    Args:
        ifgs:
        params:
        vcmt:
        mst_mat:
    """
    res = stack.stack_rate(ifgs, params, vcmt, mst_mat)
    for r in res:
        if r is None:
            raise ValueError("TODO: bad value")

    rate, error, samples = res
    write_stack_tifs(ifgs, params, res)
    # log.info('Stacked rate calculated')
    return rate, error, samples


def write_stack_tifs(ifgs, params, res):
    """

    Args:
      ifgs: param params:
      res: 
      params: 

    Returns:

    """
    # log.info('Writing stacking results')
    """
    Args:
        ifgs:
        params:
        res:
    """
    rate, error, samples = res
    gt, md, wkt = get_geotiff_header_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]
    dest = join(params[cf.OUT_DIR], "stack_rate.tif")
    md[ifc.EPOCH_DATE] = epochlist.dates
    md[ifc.DATA_TYPE] = ifc.STACKRATE
    write_output_geotiff(md, gt, wkt, rate, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "stack_error.tif")
    md[ifc.DATA_TYPE] = ifc.STACKERROR
    write_output_geotiff(md, gt, wkt, error, dest, np.nan)
    dest = join(params[cf.OUT_DIR], "stack_samples.tif")
    md[ifc.DATA_TYPE] = ifc.STACKSAMP
    write_output_geotiff(md, gt, wkt, samples, dest, np.nan)
    write_stack_numpy_files(error, rate, samples, params)


def write_stack_numpy_files(error, rate, samples, params):
    """Args:
        error:
        rate:
        samples:

    Args:
      error: param rate:
      samples: param params:
      rate: 
      params: 

    Returns:

    """
    rate_file = join(params[cf.OUT_DIR], "rate.npy")
    error_file = join(params[cf.OUT_DIR], "error.npy")
    samples_file = join(params[cf.OUT_DIR], "samples.npy")
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)


def copytree(src, dst, symlinks=False, ignore=None):
    """Copy entire contents of src directory into dst directory. See:
    http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth?lq=1

    Args:
      src(str): source directory path
      dst(str): destination directory path (created if does not exist)
      symlinks(bool, optional): Whether to copy symlink or not (Default value = False)
      ignore(bool, optional): Default value = None)

    Returns:

    """
    if not os.path.exists(dst):  # pragma: no cover
        os.makedirs(dst)
    shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if symlinks and os.path.islink(s):  # pragma: no cover
            if os.path.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except AttributeError:
                pass  # lchmod not available
        elif os.path.isdir(s):  # pragma: no cover
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def pre_prepare_ifgs(ifg_paths, params):
    """Open ifg for reading
    
    Args:
        ifg_paths:

    Args:
      ifg_paths: param params:
      params: 

    Returns:

    """
    ifgs = [Ifg(p) for p in ifg_paths]
    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        nan_and_mm_convert(i, params)
    log.info("Opened ifg for reading")
    return ifgs
