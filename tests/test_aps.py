import os
import scipy.io as sio
import numpy as np
from pyrate.algorithm import get_epochs
from pyrate.aps import tlpfilter
from pyrate import config as cf
from pyrate.compat import pickle, PY3
from tests.common import SML_TEST_DIR, TEST_CONF_GAMMA

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



def test_tlpfilter_matlab():
    epochlist = get_epochs(ifgs_pk)[0]
    params = cf.get_config_params(os.path.join(TEST_CONF_GAMMA))
    tsincr = tsincr_svd['tsincr']
    tsfilt_incr = tlpfilter(tsincr, epochlist, params)
    tsfilt_incr_matlab = ts_hp['ts_hp']
    np.testing.assert_almost_equal(tsfilt_incr_matlab,
                                   tsfilt_incr, decimal=4)


def test_tlpfilter_triangular_method_matlab():
    epochlist = get_epochs(ifgs_pk)[0]
    params = cf.get_config_params(os.path.join(TEST_CONF_GAMMA))
    params[cf.TLPF_METHOD] = 2
    tsincr = tsincr_svd['tsincr']
    tsfilt_incr = tlpfilter(tsincr, epochlist, params)
    tsfilt_incr_matlab = ts_hp_m2['ts_hp']
    np.testing.assert_almost_equal(tsfilt_incr_matlab,
                                   tsfilt_incr, decimal=4)
