import os
import copy
from collections import namedtuple
import scipy.io as sio
import numpy as np
import pytest

from pyrate.algorithm import get_epochs
from pyrate.aps import temporal_low_pass_filter as tlpfilter
from pyrate.aps import _slp_filter, spatial_low_pass_filter
from pyrate.aps import spatio_temporal_filter
from pyrate import config as cf
from pyrate.compat import pickle, PY3
from tests.common import SML_TEST_DIR, TEST_CONF_GAMMA, small_data_setup

# tsincr matrix from matlab using svd timeseries method
tsincr_svd = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                      'tsincr_svd.mat'))
# prepread ifgs pickle file
preread_pk = 'preread_ifgs.pk' if PY3 else 'preread_ifgs_py2.pk'
ifgs_pk = os.path.join(SML_TEST_DIR, 'matlab_aps', preread_pk)
ifgs_pk = pickle.load(open(ifgs_pk, 'rb'))

# tlp filter output from matlab
ts_hp = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps', 'ts_hp.mat'))
ts_hp_m2 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                    'ts_hp_m2.mat'))
ts_hp_m3 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                    'ts_hp_m3.mat'))
epochlist = get_epochs(ifgs_pk)[0]
params = cf.get_config_params(os.path.join(TEST_CONF_GAMMA))


@pytest.fixture(params=[1, 2, 3])
def tlpfilter_method(request):
    return request.param


@pytest.fixture(params=[1, 2])
def slpfilter_method(request):
    return request.param


def test_tlpfilter_matlab(tlpfilter_method):
    params[cf.TLPF_METHOD] = tlpfilter_method
    ts = ts_hp if tlpfilter_method == 1  \
        else ts_hp_m2 if tlpfilter_method == 2 \
        else ts_hp_m3
    tsincr = tsincr_svd['tsincr']
    tsfilt_incr = tlpfilter(tsincr, epochlist, params)
    tsfilt_incr_matlab = ts['ts_hp']
    np.testing.assert_almost_equal(tsfilt_incr_matlab,
                                   tsfilt_incr, decimal=4)

tsincr = tsincr_svd['tsincr']

params[cf.TLPF_METHOD] = 3
ts_hp_before_slpfilter = tsincr - tlpfilter(tsincr, epochlist, params)

# convert nan's into zeros
ts_hp_before_slpfilter[np.isnan(ts_hp_before_slpfilter)] = 0

ts_aps_m1 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                     'ts_aps.mat'))['ts_aps']
ts_aps_m2 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                     'ts_aps_m2.mat'))['ts_aps']
ts_aps_m_auto = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                         'ts_aps_auto_cutoff.mat'))['ts_aps']
xpsize = 76.834133036409  # copied from matlab since pyrate's don't match
ypsize = 92.426722191659


def test_slpfilter_matlab(slpfilter_method):
    # TODO: write tests for alpha = 0 special case when alpha (1/cutoff) is
    # computed

    params[cf.SLPF_METHOD] = slpfilter_method
    ifgs = small_data_setup()
    ts_aps_m = ts_aps_m1 if slpfilter_method == 1 else ts_aps_m2
    ts_aps = np.zeros_like(ts_aps_m)
    rows, cols = ifgs[0].shape
    for i in range(ts_aps.shape[2]):
        ts_aps[:, :, i] = _slp_filter(cutoff=params[cf.SLPF_CUTOFF],
                                      rows=rows, cols=cols,
                                      x_size=xpsize, y_size=ypsize,
                                      params=params,
                                      phase=ts_hp_before_slpfilter[:, :, i])
    for i in ifgs:
        i.close()
    np.testing.assert_array_almost_equal(ts_aps, ts_aps_m, decimal=4)


def test_slpfilter_accumulated(slpfilter_method):
    ts_aps_before = copy.copy(ts_hp_before_slpfilter)
    params[cf.SLPF_METHOD] = slpfilter_method
    ifgs = small_data_setup()
    ts_aps_m = ts_aps_m1 if slpfilter_method == 1 else ts_aps_m2
    Ifg = namedtuple('Ifg', 'x_centre, y_centre, x_size, y_size, shape')
    ifg = Ifg(x_size=xpsize, y_size=ypsize, shape=ifgs[0].shape,
              x_centre=ifgs[0].x_centre, y_centre=ifgs[0].y_centre)
    ts_aps = spatial_low_pass_filter(ts_aps_before,
                                     ifg, params=params)

    for i in ifgs:
        i.close()
    np.testing.assert_array_almost_equal(ts_aps, ts_aps_m, decimal=4)


def test_slpfilter_auto_cutoff(slpfilter_method=2):
    ts_aps_before = copy.copy(ts_hp_before_slpfilter)
    params[cf.SLPF_METHOD] = slpfilter_method
    params[cf.SLPF_CUTOFF] = 0
    ifgs = small_data_setup()
    Ifg = namedtuple('Ifg', 'x_centre, y_centre, x_size, y_size, shape')
    ifg = Ifg(x_size=xpsize, y_size=ypsize, shape=ifgs[0].shape,
              x_centre=ifgs[0].x_centre, y_centre=ifgs[0].y_centre)

    ts_aps = spatial_low_pass_filter(ts_aps_before,
                                     ifg, params=params)

    for i in ifgs:
        i.close()
    np.testing.assert_array_almost_equal(ts_aps, ts_aps_m_auto, decimal=4)


def test_spatio_temporal_filter():
    import tempfile
    from pyrate import shared
    from os.path import basename, join
    from collections import OrderedDict
    from pyrate.scripts import run_prepifg
    from osgeo import gdal
    from pyrate import ifgconstants as ifc
    ifg_out = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                          'ifg_spatio_temp_out.mat'))['ifg']

    tsincr = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                      'tsincr_svd.mat'))['tsincr']
    params = cf.get_config_params(os.path.join(TEST_CONF_GAMMA))
    params[cf.OUT_DIR] = tempfile.mkdtemp()
    params[cf.TMPDIR] = join(params[cf.OUT_DIR], cf.TMPDIR)
    shared.mkdir_p(params[cf.TMPDIR])
    params[cf.SLPF_METHOD] = 2
    params[cf.SLPF_CUTOFF] = 0
    params[cf.SLPF_ORDER] = 1
    params[cf.SLPF_NANFILL] = 0
    params[cf.TLPF_METHOD] = 3
    params[cf.TLPF_CUTOFF] = 0.25
    params[cf.TLPF_PTHR] = 5
    ifgs = small_data_setup()
    _ = [ifgs_pk.pop(k) for k in ['gt', 'epochlist', 'md', 'wkt']]
    preread_ifgs = {join(params[cf.OUT_DIR], basename(k)): v
                    for k, v in ifgs_pk.items()}
    preread_ifgs = OrderedDict(sorted(preread_ifgs.items()))

    run_prepifg.main(params)
    for k, v in preread_ifgs.items():
        v.path = join(params[cf.OUT_DIR], basename(k))
        preread_ifgs[k] = v
    ifg = ifgs[0]
    ifg.x_size = xpsize
    ifg.y_size = ypsize
    spatio_temporal_filter(tsincr, ifg, params, preread_ifgs)
    for i in ifgs:
        i.close()

    for (p, v), i in zip(preread_ifgs.items(), range(ifg_out.shape[2])):
        ds = gdal.Open(p)
        metadata = ds.GetMetadata()
        assert ifc.PYRATE_APS_ERROR in metadata
        assert metadata[ifc.PYRATE_APS_ERROR] == ifc.APS_REMOVED
        arr = ds.GetRasterBand(1).ReadAsArray()
        np.testing.assert_array_almost_equal(arr, ifg_out[:, :, i], decimal=3)


@pytest.fixture(params=['cubic', 'linear', 'nearest'])
def interp_method(request):
    return request.param


def test_interpolate_nans(interp_method):
    from pyrate.aps import _interpolate_nans
    from copy import copy
    a = np.arange(25*3).reshape((5, 5, 3)).astype(float)
    a[np.random.randint(2, size=(5, 5, 3)).astype(bool)] = np.nan
    a_copy = copy(a)

    _interpolate_nans(a, method=interp_method)
    assert np.sum(np.isnan(a)) == 0
    np.testing.assert_array_almost_equal(a[~np.isnan(a_copy)],
                                         a_copy[~np.isnan(a_copy)],
                                         decimal=4)
