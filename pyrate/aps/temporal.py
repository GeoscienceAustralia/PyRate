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
Temporal low pass filter.
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import logging
import numpy as np
from numpy import isnan
from pyrate import config as cf

log = logging.getLogger(__name__)


# TODO: use tiles here and distribute amongst processes
def tlpfilter(tsincr, epochlist, params):
    """
    Temporal low pass filter

    Parameters
    ----------
    tsincr : ndarray
        time series array of size (ifg.shape, nvels)
    epochlist : list
        list of shared.EpochList class instances
    params : dict
        parameters dict

    Returns
    -------
    tsfilt_incr : ndarray
        temporal filtered time series of shape (ifg.shape, nvels)
    """
    log.info('Applying temporal low pass filter')
    nanmat = ~isnan(tsincr)
    tsfilt_incr = np.empty_like(tsincr, dtype=np.float32) * np.nan
    intv = np.diff(epochlist.spans)  # time interval for the neighboring epoch
    span = epochlist.spans[: tsincr.shape[2]] + intv/2  # accumulated time
    rows, cols = tsincr.shape[:2]
    cutoff = params[cf.TLPF_CUTOFF]
    method = params[cf.TLPF_METHOD]
    threshold = params[cf.TLPF_PTHR]
    if method == 1:  # gaussian filter
        func = gauss
    elif method == 2:  # triangular filter
        func = _triangle
    else:
        func = mean_filter

    _tlpfilter(cols, cutoff, nanmat, rows, span, threshold, tsfilt_incr,
               tsincr, func)
    log.info('Finished applying temporal low pass filter')
    return tsfilt_incr


gauss = lambda m, yr, cutoff: np.exp(-(yr / cutoff) ** 2 / 2)


def _triangle(m, yr, cutoff):
    # pylint: disable=unused-argument
    wgt = cutoff - abs(yr)
    wgt[wgt < 0] = 0
    return wgt


mean_filter = lambda m, yr, cutoff: np.ones(m)


def _tlpfilter(cols, cutoff, nanmat, rows, span, threshold, tsfilt_incr,
               tsincr, func):
    for i in range(rows):
        for j in range(cols):
            sel = np.nonzero(nanmat[i, j, :])[0]  # don't select if nan
            m = len(sel)
            if m >= threshold:
                for k in range(m):
                    yr = span[sel] - span[sel[k]]
                    wgt = func(m, yr, cutoff)
                    wgt /= np.sum(wgt)
                    tsfilt_incr[i, j, sel[k]] = np.sum(tsincr[i, j, sel] * wgt)
