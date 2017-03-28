import os
import scipy.io as sio
import numpy as np
import pytest

from pyrate.algorithm import get_epochs
from pyrate.aps.temporal import tlpfilter
from pyrate.aps.spatial import _slp_filter
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
ts_aps_m1 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                     'ts_aps.mat'))['ts_aps']
ts_aps_m2 = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                     'ts_aps_m2.mat'))['ts_aps']


def test_slpfilter_matlab(slpfilter_method):
    # TODO: write tests for alpha = 0 special case when alpha (1/cutoff) is
    # computed

    params[cf.TLPF_METHOD] = slpfilter_method
    ifgs = small_data_setup()
    xpsize = 76.834133036409  # copied from matlab since pyrate's don't match
    ypsize = 92.426722191659
    ts_aps = np.zeros_like(ts_aps_m1)
    rows, cols = ifgs[0].shape
    for i in range(ts_aps.shape[2]):
        ts_aps[:, :, i] = _slp_filter(cutoff=params[cf.SLPF_CUTOFF],
                                      rows=rows, cols=cols,
                                      x_size=xpsize, y_size=ypsize,
                                      params=params,
                                      phase=ts_hp_before_slpfilter[:, :, i])

    np.testing.assert_array_almost_equal(ts_aps, ts_aps_m1, decimal=4)
