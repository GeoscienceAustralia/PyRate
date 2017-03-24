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
    if np.sum(np.isnan(phase)) == 0:  # if there are no nans in the data
        return ifg
    cutoff = params[cf.SLPF_CUTOFF]
    if cutoff == 0:
        maxvar, alpha = cvd_from_phase(phase, ifg, calc_alpha=True)
        cutoff = 1/alpha

    # find the center of the imag
    rows, cols = ifg.shape
    cx = ifg.x_centre
    cy = ifg.y_centre

    # fft for the input image
    phase[np.isnan(phase)] = 0
    imf = fftshift(fft2(phase))

    # calculate distance
    distfact = 1000.0  # to convert into km
    [xx, yy] = np.meshgrid(range(cols), range(rows))
    xx = (xx - cx) * ifg.x_size
    yy = (yy-cy) * ifg.y_size
    dist = np.divide(np.sqrt(((xx - ifg.x_centre) * ifg.x_size) ** 2 +
                     ((yy - ifg.y_centre) * ifg.y_size) ** 2),
                     distfact)

    if params[cf.SLPF_METHOD] == 1:  # butterworth low pass filter
        H = 1./(1+((dist / cutoff) ** (2*params[cf.SLPF_ORDER])))
    else:  # Gaussian lowpass filter
        H = np.exp(-(dist ** 2)/(2 * cutoff**2))

    outf = imf * H
    out = np.real(ifft2(ifftshift(outf)))
    out[np.isnan(phase)] = np.nan

    return out
