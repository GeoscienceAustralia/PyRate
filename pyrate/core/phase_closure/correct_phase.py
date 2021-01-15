import numpy as np
from typing import Any
from nptyping import NDArray, Float32, Bool
import pyrate.constants as C
from pyrate.core.logger import pyratelogger as log


def recover_pixels(closure: NDArray[(Any, Any), Float32], params: dict) -> NDArray[(Any, ), Bool]:
    """
    :param closure: sum `closure` of each loop.
    :param params: params dict
    :return: indices_breaching_threshold, potentially corrected if opted
    """
    large_dev_thr = params[C.LARGE_DEV_THR] * np.pi

    if params[C.CORR_PHASE]:
        log.debug("Correcting phase closure!")
        nans = np.isnan(closure)
        closure_plus_2pi = np.abs(closure)  # initialize
        pixels_under_thr = closure_plus_2pi < large_dev_thr  # nan pixels are not selected
        max_closure = np.nanmax(np.abs(closure))
        if max_closure > 2 * np.pi - large_dev_thr:
            log.debug(f"Maximum closure value of {max_closure} detected! Will attempt phase correction!")
            multiples_of_2pi = int(max_closure/2/np.pi) + 1
            # max_closure is already computed, use it
            if max_closure > 2 * np.pi - large_dev_thr:  # otherwise there is no point doing the check
                for m in range(multiples_of_2pi):
                    closure_plus_2pi -= 2 * np.pi
                    recovered_pixels = np.abs(closure_plus_2pi) < large_dev_thr
                    num_pixels_recovered = np.sum(recovered_pixels)
                    log.debug(f"Recovered {num_pixels_recovered} pixels using a further 2*pi correction!")
                    pixels_under_thr += recovered_pixels

            # ~ pixels_under_thr contains nans, we know nothing about these pixels, so discard
            indices_breaching_threshold = np.logical_and(~pixels_under_thr, ~nans)
        else:
            log.debug(f"Maximum closure value of {max_closure} detected! Correction will not be attempted!")
            indices_breaching_threshold = np.absolute(closure) > large_dev_thr

    else:
        # this will deal with nans in `closure`, i.e., nans are not selected in indices_breaching_threshold
        indices_breaching_threshold = np.absolute(closure) > large_dev_thr

    return indices_breaching_threshold
