import logging
import numpy as np
from scipy.fftpack import fft2, ifft2, fftshift, ifftshift
from pyrate import config as cf
from pyrate.vcm import cvd_from_phase

log = logging.getLogger(__name__)


def spatial_low_pass_filter(ts_hp, ifg, params):
    log.info('Applying spatial low pass filter')
    for i in range(ts_hp.shape[2]):
        ts_hp[:, :, i] = slpfilter(ts_hp[:, :, i], ifg, params)

    log.info('Finished applying spatial low pass filter')
    return ts_hp


def slpfilter(phase, ifg, params):
    if np.all(np.isnan(phase)):  # return for nan matrix
        return phase
    cutoff = params[cf.SLPF_CUTOFF]
    phase[np.isnan(phase)] = 0  # need it here for cvd calc
    if cutoff == 0:
        maxvar, alpha = cvd_from_phase(phase, ifg, calc_alpha=True)
        cutoff = 1.0/alpha
    rows, cols = ifg.shape
    return _slp_filter(phase, cutoff, rows, cols,
                       ifg.x_size, ifg.y_size, params)


def _slp_filter(phase, cutoff, rows, cols, x_size, y_size, params):
    cx = np.floor(cols/2)
    cy = np.floor(rows/2)
    # fft for the input image
    phase[np.isnan(phase)] = 0
    imf = fftshift(fft2(phase))
    # calculate distance
    distfact = 1.0e3  # to convert into meters
    [xx, yy] = np.meshgrid(range(cols), range(rows))
    xx = (xx - cx) * x_size
    yy = (yy - cy) * y_size
    dist = np.sqrt(xx ** 2 + yy ** 2)/distfact

    if params[cf.SLPF_METHOD] == 1:  # butterworth low pass filter
        H = 1. / (1 + ((dist / cutoff) ** (2 * params[cf.SLPF_ORDER])))
    else:  # Gaussian low pass filter
        H = np.exp(-(dist ** 2) / (2 * cutoff ** 2))
    outf = imf * H
    out = np.real(ifft2(ifftshift(outf)))
    out[np.isnan(phase)] = np.nan
    return out
