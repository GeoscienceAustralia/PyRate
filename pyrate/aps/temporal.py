import logging
import numpy as np
from numpy import isnan
from pyrate import config as cf

log = logging.getLogger(__name__)


# TODO: use tiles here and distribute amongst processes
def tlpfilter(tsincr, epochlist, params):
    log.info('Applying temporal low pass filter')

    """temporal low pass filter"""
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
        func = triangle
    else:
        func = mean_filter

    _tlpfilter(cols, cutoff, nanmat, rows, span, threshold, tsfilt_incr,
               tsincr, func)
    log.info('Finished applying temporal low pass filter')
    return tsfilt_incr


gauss = lambda m, yr, cutoff: np.exp(-(yr / cutoff) ** 2 / 2)


def triangle(m, yr, cutoff):
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
