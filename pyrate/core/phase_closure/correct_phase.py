#   This Python module is part of the PyRate software package.
#
#   Copyright 2021 Geoscience Australia
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

import numpy as np
from typing import Any
from nptyping import NDArray, Float32, Bool
import pyrate.constants as C
from pyrate.core.logger import pyratelogger as log


def find_breached_pixels(closure: NDArray[(Any, Any), Float32], params: dict) -> NDArray[(Any,), Bool]:
    """
    :param closure: sum `closure` of each loop.
    :param params: params dict
    :return: indices_breaching_threshold, potentially corrected if opted
    """
    large_dev_thr = params[C.LARGE_DEV_THR] * np.pi

    if params[C.CORR_PHASE]:
        log.debug("Correcting phase closure!")
        nans = np.isnan(closure)
        abs_closure = np.abs(closure)  # initialize
        pixels_under_thr = abs_closure < large_dev_thr  # nan pixels are not selected
        max_closure = np.nanmax(abs_closure)
        # max_closure is already computed, use it
        if max_closure > 2 * np.pi - large_dev_thr:  # otherwise there is no point doing the check
            recovered_pixels = __recover_pixels(abs_closure, large_dev_thr, max_closure)
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


def __recover_pixels(abs_closure: NDArray[(Any, Any), Float32],
                     large_dev_thr,
                     max_closure) -> NDArray[(Any, Any), Bool]:
    log.debug(f"Maximum closure value of {max_closure} detected! Will attempt phase correction!")
    multiples_of_2pi = int(max_closure / 2 / np.pi) + 1

    recovered_pixels = np.zeros_like(abs_closure, dtype=np.bool)
    for m in range(multiples_of_2pi):
        abs_closure -= 2 * np.pi
        this_loop_recovered = np.abs(abs_closure) < large_dev_thr
        recovered_pixels += this_loop_recovered
        log.debug(f"Recovered {np.sum(this_loop_recovered)} pixels using a further 2*pi correction!")
    log.info(f"A total of {np.sum(recovered_pixels)} were recovered using phase closure correction")
    return recovered_pixels
