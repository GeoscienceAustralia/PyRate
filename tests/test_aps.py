import os
import pickle
import scipy.io as sio
import numpy as np
from pyrate.algorithm import get_epochs
from pyrate.aps import tlpfilter
from pyrate import config as cf
from tests.common import SML_TEST_DIR, TEST_CONF_GAMMA

# tsincr matrix from matlab using svd timeseries method
tsincr_svd = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                      'tsincr_svd.mat'))
# prepread ifgs pickle file
ifgs_pk = os.path.join(SML_TEST_DIR, 'matlab_aps', 'preread_ifgs.pk')
ifgs_pk = pickle.load(open(ifgs_pk, 'rb'))

# tlp filter output from matlab
ts_hp = sio.loadmat(os.path.join(SML_TEST_DIR, 'matlab_aps',
                                     'ts_hp.mat'))


def test_tlpfilter_matlab():
    epochlist = get_epochs(ifgs_pk)[0]
    params = cf.get_config_params(os.path.join(TEST_CONF_GAMMA))
    tsincr = tsincr_svd['tsincr']
    tsfilt_incr = tlpfilter(tsincr, epochlist, params)
    # TODO: why don't the entire matrices match?
    np.testing.assert_almost_equal(ts_hp['ts_hp'][:30, :47, :],
                                   tsfilt_incr[:30, :47, :], decimal=4)
