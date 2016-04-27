from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

cdef greater(float x, float y):
    return x > y

@cython.boundscheck(False)
def resample(np.ndarray[DTYPE_t, ndim=2] data, int xscale, int yscale, double thresh):
    """
    Resamples/averages 'data' to return an array from the averaging of blocks
    of several tiles in 'data'. NB: Assumes incoherent cells are NaNs.

    :param data: source array to resample to different size
    :param xscale: number of cells to average along X axis
    :param yscale: number of Y axis cells to average
    :param thresh: minimum allowable proportion of NaN cells (range from 0.0-1.0),
        eg. 0.25 = 1/4 or more as NaNs results in a NaN value for the output cell.
    """
    cdef unsigned int x, y
    cdef unsigned int ysize = data.shape[0]
    cdef unsigned int xsize = data.shape[1]
    cdef int tile_cell_count, xres, yres
    cdef DTYPE_t nan_fraction
    cdef np.ndarray[DTYPE_t, ndim=2] dest, tile

    xres = int(xsize / xscale)
    yres = int(ysize / yscale)
    dest = np.zeros((yres, xres), dtype=np.float32) * np.nan
    tile_cell_count = xscale * yscale

    # calc mean without nans (fractional threshold ignores tiles with excess NaNs)
    # for y, x in itertools.product(xrange(yres), xrange(xres)):
    for y in range(yres):
        for x in range(xres):
            tile = data[y * yscale: (y+1) * yscale, x * xscale: (x+1) * xscale]
            nan_fraction = np.sum(np.isnan(tile)) / tile_cell_count

            if greater(thresh, nan_fraction) or (nan_fraction == 0 and thresh == 0):
                dest[y, x] = np.nanmean(tile)

    return dest