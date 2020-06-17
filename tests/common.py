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
from os.path import join
from subprocess import check_output
from pathlib import Path

import numpy as np
from numpy import isnan, sum as nsum
from osgeo import gdal

from pyrate.core import algorithm, ifgconstants as ifc, config as cf, timeseries, mst, stack
from pyrate.core.shared import (Ifg, nan_and_mm_convert, get_geotiff_header_info,
                                write_output_geotiff, dem_or_ifg)
from pyrate.constants import PYRATEPATH
from pyrate.configuration import Configuration

TRAVIS = True if 'TRAVIS' in os.environ else False
PYTHON3P6 = True if ('TRAVIS_PYTHON_VERSION' in os.environ and os.environ['TRAVIS_PYTHON_VERSION'] == '3.6') else False
PYTHON3P7 = True if ('TRAVIS_PYTHON_VERSION' in os.environ and os.environ['TRAVIS_PYTHON_VERSION'] == '3.7') else False
PYTHON3P8 = True if ('TRAVIS_PYTHON_VERSION' in os.environ and os.environ['TRAVIS_PYTHON_VERSION'] == '3.8') else False
GDAL_VERSION = check_output(["gdal-config", "--version"]).decode(encoding="utf-8").split('\n')[0]


TEMPDIR = tempfile.gettempdir()
TESTDIR = join(PYRATEPATH, 'tests')
BASE_TEST = join(PYRATEPATH, "tests", "test_data")
SML_TEST_DIR = join(BASE_TEST, "small_test")
SML_TEST_OBS = join(SML_TEST_DIR, 'roipac_obs')  # roipac processed unws
SML_TEST_OUT = join(SML_TEST_DIR, 'out')
SML_TEST_TIF = join(SML_TEST_DIR, 'tif')
SML_TEST_GAMMA = join(SML_TEST_DIR, 'gamma_obs')  # gamma processed unws
SML_TEST_ROIPAC = join(SML_TEST_DIR, 'roipac_obs')  # gamma processed unws
SML_TEST_CONF = join(SML_TEST_DIR, 'conf')
SML_TEST_GAMMA_HEADER_LIST = join(SML_TEST_GAMMA, 'headers')
SML_TEST_ROIPAC_HEADER_LIST = join(SML_TEST_ROIPAC, 'headers')

SML_TEST_DEM_DIR = join(SML_TEST_DIR, 'dem')
SML_TEST_LEGACY_PREPIFG_DIR = join(SML_TEST_DIR, 'prepifg_output')
SML_TEST_LEGACY_ORBITAL_DIR = join(SML_TEST_DIR, 'orbital_error_correction')
SML_TEST_DEM_ROIPAC = join(SML_TEST_OBS, 'roipac_test_trimmed.dem')
SML_TEST_DEM_GAMMA = join(SML_TEST_GAMMA, '20060619_utm.dem')
SML_TEST_INCIDENCE = join(SML_TEST_GAMMA, '20060619_utm.inc')
SML_TEST_ELEVATION = join(SML_TEST_GAMMA, '20060619_utm.lv_theta')
SML_TEST_DEM_HDR_GAMMA = join(SML_TEST_GAMMA, '20060619_utm_dem.par')
SML_TEST_DEM_HDR = join(SML_TEST_OBS, 'roipac_test_trimmed.dem.rsc')
SML_TEST_DEM_TIF = join(SML_TEST_DEM_DIR, 'roipac_test_trimmed.tif')

SML_TEST_COH_DIR = join(SML_TEST_DIR, 'coherence')
SML_TEST_COH_LIST = join(SML_TEST_COH_DIR, 'coherence_17')

TEST_CONF_ROIPAC = join(SML_TEST_CONF, 'pyrate_roipac_test.conf')
TEST_CONF_GAMMA = join(SML_TEST_CONF, 'pyrate_gamma_test.conf')

ROIPAC_SYSTEM_CONF = PYRATEPATH.joinpath("tests", "test_data", "system", "roipac", "input_parameters.conf")
GAMMA_SYSTEM_CONF = PYRATEPATH.joinpath("tests", "test_data", "system", "gamma", "input_parameters.conf")
GEOTIF_SYSTEM_CONF = PYRATEPATH.joinpath("tests", "test_data", "system", "geotiff", "input_parameters.conf")


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
    tifs = glob.glob(os.path.join(path, '*.tif'))
    for tif in tifs:
        os.remove(tif)


def small_data_setup(datafiles=None, is_dir=False):
    """Returns Ifg objs for the files in the small test dir
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if is_dir:
        datafiles = glob.glob(join(datafiles, "*.tif"))
    else:
        if datafiles:
            for i, d in enumerate(datafiles):
                datafiles[i] = os.path.join(SML_TEST_TIF, d)
        else:
            datafiles = glob.glob(join(SML_TEST_TIF, "*.tif"))
    datafiles.sort()
    ifgs = [dem_or_ifg(i) for i in datafiles]
    
    for i in ifgs: 
        i.open()
        i.nodata_value = 0

    return ifgs


def assert_tifs_equal(tif1, tif2):
    mds = gdal.Open(tif1)
    sds = gdal.Open(tif2)

    md_mds = mds.GetMetadata()
    md_sds = sds.GetMetadata()
    # meta data equal
    assert md_mds == md_sds

    d1 = mds.ReadAsArray()
    d2 = sds.ReadAsArray()
    # phase equal
    np.testing.assert_array_almost_equal(d1,  d2, decimal=3)

    mds = None  # close datasets
    sds = None


def copy_small_ifg_file_list():
    temp_dir = tempfile.mkdtemp()
    move_files(SML_TEST_TIF, temp_dir, file_type='*.tif', copy=True)
    datafiles = glob.glob(join(temp_dir, "*.tif"))
    for d in datafiles:
        Path(d).chmod(0o664)  # assign write permission as conv2tif output is readonly
    return temp_dir, datafiles


def copy_and_setup_small_data():
    temp_dir, datafiles = copy_small_ifg_file_list()
    datafiles.sort()
    ifgs = [dem_or_ifg(i) for i in datafiles]

    for i in ifgs:
        i.open()
        i.nodata_value = 0
    return temp_dir, ifgs


def small_ifg_file_list(datafiles=None):
    """Returns the file list of all the .tif files after prepifg conversion
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    if datafiles:
        for i, d in enumerate(datafiles):
            datafiles[i] = os.path.join(SML_TEST_TIF, d)
    else:
        datafiles = glob.glob(join(SML_TEST_TIF, "*.tif"))
    datafiles.sort()
    return datafiles


def small_data_roipac_unws():
    """Returns unw file list before prepifg operation
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    return glob.glob(join(SML_TEST_OBS, "*.unw"))


def small_data_setup_gamma_unws():
    """Returns unw file list before prepifg operation
    input phase data is in radians; these ifgs are in radians - not converted to mm"""
    return glob.glob(join(SML_TEST_GAMMA, "*.unw"))


def small5_ifgs():
    """Convenience func to return a subset of 5 linked Ifgs from the testdata"""
    BASE_DIR = tempfile.mkdtemp()
    data_paths = [os.path.join(SML_TEST_TIF, p) for p in IFMS5.split()]
    new_data_paths = [os.path.join(BASE_DIR, os.path.basename(d))
                      for d in data_paths]
    for d in data_paths:
        shutil.copy(d, os.path.join(BASE_DIR, os.path.basename(d)))

    return [Ifg(p) for p in new_data_paths]


def small5_mock_ifgs(xs=3, ys=4):
    '''Returns smaller mocked version of small Ifgs for testing'''
    ifgs = small5_ifgs()
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


def reconstruct_stack_rate(shape, tiles, output_dir, out_type):
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


def move_files(source_dir, dest_dir, file_type='*.tif', copy=False):
    for filename in glob.glob(os.path.join(source_dir, file_type)):
        if copy:
            shutil.copy(filename, dest_dir)
        else:
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
        # write mst output to a file
        mst_mat_binary_file = join(params[cf.OUT_DIR], 'mst_mat')
        np.save(file=mst_mat_binary_file, arr=mst_grid)

        for i in ifgs:
            i.close()
        return mst_grid
    return None


def get_nml(ifg_list_instance, nodata_value, nan_conversion=False):
    """
    :param xxx(eg str, tuple, int, float...) ifg_list_instance: xxxx
    :param float nodata_value: No data value in image
    :param bool nan_conversion: Convert NaNs
    
    :return: ifg_list_instance: replaces in place
    :rtype: list
    :return: _epoch_list: list of epochs
    :rtype: list
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


def write_timeseries_geotiff(ifgs, params, tsincr, pr_type):
    # setup metadata for writing into result files
    gt, md, wkt = get_geotiff_header_info(ifgs[0].data_path)
    epochlist = algorithm.get_epochs(ifgs)[0]

    for i in range(tsincr.shape[2]):
        md[ifc.EPOCH_DATE] = epochlist.dates[i + 1]
        md['SEQUENCE_POSITION'] = i+1  # sequence position

        data = tsincr[:, :, i]
        dest = join(params[cf.OUT_DIR], pr_type + "_" +
                    str(epochlist.dates[i + 1]) + ".tif")
        md[ifc.DATA_TYPE] = pr_type
        write_output_geotiff(md, gt, wkt, data, dest, np.nan)


def calculate_stack_rate(ifgs, params, vcmt, mst_mat=None):
    # log.info('Calculating stacked rate')
    res = stack.stack_rate_array(ifgs, params, vcmt, mst_mat)
    for r in res:
        if r is None:
            raise ValueError('TODO: bad value')

    r, e, samples = res
    rate, error = stack.mask_rate(r, e, params['maxsig'])
    write_stackrate_tifs(ifgs, params, res)
    # log.info('Stacked rate calculated')
    return rate, error, samples


def write_stackrate_tifs(ifgs, params, res):
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
    write_stackrate_numpy_files(error, rate, samples, params)


def write_stackrate_numpy_files(error, rate, samples, params):
    rate_file = join(params[cf.OUT_DIR], 'rate.npy')
    error_file = join(params[cf.OUT_DIR], 'error.npy')
    samples_file = join(params[cf.OUT_DIR], 'samples.npy')
    np.save(file=rate_file, arr=rate)
    np.save(file=error_file, arr=error)
    np.save(file=samples_file, arr=samples)


def copytree(src, dst, symlinks=False, ignore=None):
    # pylint: disable=line-too-long
    """
    Copy entire contents of src directory into dst directory.
    See: http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth?lq=1

    :param str src: source directory path
    :param str dst: destination directory path (created if does not exist)
    :param bool symlinks: Whether to copy symlink or not
    :param bool ignore:
    """
    # pylint: disable=invalid-name
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
    """
    nan and mm convert ifgs
    """
    ifgs = [Ifg(p) for p in ifg_paths]
    for i in ifgs:
        if not i.is_open:
            i.open(readonly=False)
        nan_and_mm_convert(i, params)
    return ifgs


def assert_two_dirs_equal(dir1, dir2, ext, num_files=None):

    dir1_files = list(Path(dir1).glob(ext))
    dir2_files = list(Path(dir2).glob(ext))  # MultiProcess files
    dir1_files.sort()
    dir2_files.sort()
    # 17 unwrapped geotifs
    # 17 cropped multilooked tifs + 1 dem
    if num_files is not None:
        assert len(dir1_files) == num_files
        assert len(dir2_files) == num_files
    else:
        assert len(dir1_files) == len(dir2_files)
    if dir1_files[0].suffix == '.tif':
        for m_f, s_f in zip(dir1_files, dir2_files):
            assert m_f.name == s_f.name
            assert_tifs_equal(m_f.as_posix(), s_f.as_posix())

    elif dir1_files[0].suffix == '.npy':
        for m_f, s_f in zip(dir1_files, dir2_files):
            assert m_f.name == s_f.name
            np.testing.assert_array_almost_equal(np.load(m_f), np.load(s_f))
    elif dir1_files[0].suffix in {'.kml', '.png'}:
        return
    else:
        raise


def assert_same_files_produced(dir1, dir2, dir3, ext, num_files=None):
    assert_two_dirs_equal(dir1, dir2, ext, num_files)
    assert_two_dirs_equal(dir1, dir3, ext, num_files)


def manipulate_test_conf(conf_file, temp_obs_dir):
    params = cf.get_config_params(conf_file)
    copytree(params[cf.OBS_DIR], temp_obs_dir)
    # manipulate params
    params[cf.OBS_DIR] = temp_obs_dir.as_posix()
    outdir = temp_obs_dir.joinpath('out')
    outdir.mkdir(exist_ok=True)
    params[cf.OUT_DIR] = outdir.as_posix()
    params[cf.DEM_FILE] = temp_obs_dir.joinpath(Path(params[cf.DEM_FILE]).name).as_posix()
    params[cf.DEM_HEADER_FILE] = temp_obs_dir.joinpath(Path(params[cf.DEM_HEADER_FILE]).name).as_posix()
    params[cf.HDR_FILE_LIST] = temp_obs_dir.joinpath(Path(params[cf.HDR_FILE_LIST]).name).as_posix()
    params[cf.SLC_DIR] = temp_obs_dir.as_posix()
    params[cf.IFG_FILE_LIST] = temp_obs_dir.joinpath(Path(params[cf.IFG_FILE_LIST]).name).as_posix()
    params[cf.COH_FILE_DIR] = temp_obs_dir.as_posix()
    params[cf.APS_INCIDENCE_MAP] = temp_obs_dir.joinpath(Path(params[cf.APS_INCIDENCE_MAP]).name).as_posix()
    params[cf.TMPDIR] = temp_obs_dir.joinpath(Path(params[cf.TMPDIR]).name).as_posix()
    return params
