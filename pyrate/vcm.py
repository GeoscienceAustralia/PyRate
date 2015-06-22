'''
Variance/Covariance matrix functionality for PyRate, ported from Hua Wang's
MATLAB code for Pirate.

.. codeauthor:: Ben Davies, Matt Garthwaite
'''

from copy import copy
from numpy import array, where, isnan, real, imag, sum, sqrt, meshgrid
from numpy import zeros, vstack, ceil, mean, exp, reshape
from numpy.linalg import norm
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.optimize import fmin

from pyrate.algorithm import master_slave_ids


def pendiffexp(alphamod, cvdav):
    '''
    Fits an exponential model to data.

    :param float alphamod: Exponential decay exponent.
    :param array cvdav: Function magnitude at 0 radius (2 col array of radius, variance)
    '''

    # maxvar usually at zero lag
    mx = cvdav[1,0]
    return norm(cvdav[1,:] - (mx * exp(-alphamod * cvdav[0,:])))


def unique_points(points):
    '''
    Returns unique points from a list of coordinates.

    :param points: Sequence of (y,x) or (x,y) tuples.
    '''

    return vstack([array(u) for u in set(points) ] )


def cvd(ifg):
    '''
    Calculate average covariance versus distance (autocorrelation) and its best fitting exponential function

    :param ifg: An interferogram.
    :type ifg: :py:class:`pyrate.shared.Ifg`.
    '''
    # distance division factor of 1000 converts to km and is needed to match Matlab code output
    distfact = 1000

    # calculate 2D auto-correlation of image using the spectral method (Wiener-Khinchin theorem)
    phase = where(isnan(ifg.phase_data), 0, ifg.phase_data)
    fft_phase = fft2(phase)
    pspec = real(fft_phase)**2 + imag(fft_phase)**2
    autocorr_grid = ifft2(pspec)
    nzc = sum(sum(phase != 0))
    autocorr_grid = fftshift(real(autocorr_grid)) / nzc
        
    # pixel distances from pixel at zero lag (image centre).
    xx,yy = meshgrid(range(ifg.ncols), range(ifg.nrows))
    r = sqrt(((xx-ifg.x_centre) * ifg.x_size)**2 + ((yy-ifg.y_centre) * ifg.y_size)**2) / distfact
  
    r = reshape(r, ifg.num_cells, 1)
    acg = reshape(autocorr_grid, ifg.num_cells, 1)    
  
    # Symmetry in image; keep only unique points
    tmp = unique_points(zip(acg, r))
    acg, r = tmp[:,0], tmp[:,1]
    
    # Alternative method to remove duplicate cells from Matlab Pirate
    #r = r[:ceil(len(r)/2)+nlines] # Reason for '+nlines' term unknown

    # TODO: write a better way of getting rid of duplicate {r, acg} pairs - sets?
    # eg. array([x for x in set([(1,1), (2,2), (1,1)])])
    # the above shortens r by some number of cells
    #print 'r',r

    acgorig = copy(acg)
    rorig = copy(r)
    
    # bin width for collecting data
    w = max(ifg.x_size, ifg.y_size) * 2 / distfact
    
    # pick the smallest axis to determine circle search radius
    #print 'ifg.X_CENTRE, ifg.Y_CENTRE=', ifg.x_centre, ifg.y_centre
    #print 'ifg.X_SIZE, ifg.Y_SIZE', ifg.x_size, ifg.y_size
    if (ifg.x_centre * ifg.x_size) < (ifg.y_centre * ifg.y_size):
        maxdist = ifg.x_centre * ifg.x_size / distfact
        #print "maxdist from X axis is", maxdist
    else:
        maxdist = ifg.y_centre * ifg.y_size/ distfact
        #print "maxdist from Y axis is", maxdist
    
    # filter out data where the of lag distance is greater than maxdist
    #r = array([e for e in rorig if e <= maxdist]) # MG: prefers to use all the data
    #acg = array([e for e in rorig if e <= maxdist])    
    r = rorig[rorig<maxdist]
    acg = acgorig[rorig<maxdist]
       
    # classify values of r according to bin number
    rbin = ceil(r / w).astype(int)
    maxbin = max(rbin) # consistent with Matlab code
    
    cvdav = zeros( (2, maxbin) )

    for b in range(maxbin):
        cvdav[0,b] = b * w # distance instead of bin number
        cvdav[1,b] = mean(acg[rbin == b]) # mean variance for that bin

    # calculate best fit function maxvar*exp(-alpha*r)
    alphaguess = 2 / (maxbin * w)
    alpha = fmin(pendiffexp, x0=alphaguess, args=(cvdav,), disp=0 )
    print "1st guess, alpha", alphaguess, alpha

    assert len(alpha) == 1
    # maximum variance usually at the zero lag
    maxvar = max(acg[:len(r)])
    return maxvar, alpha[0]


def get_vcmt(ifgs, maxvar):
    '''
    Returns the temporal variance/covariance matrix.
    '''

    # c=0.5 for common master or slave; c=-0.5 if master of one matches slave of another
    nifgs = len(ifgs)
    vcm_pat = zeros((nifgs, nifgs))

    dates = [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]
    ids = master_slave_ids(dates)

    for i, ifg in enumerate(ifgs):
        mas1, slv1 = ids[ifg.master], ids[ifg.slave]

        for j, ifg2 in enumerate(ifgs):
            mas2, slv2 = ids[ifg2.master], ids[ifg2.slave]
            if mas1 == mas2 or slv1 == slv2:
                vcm_pat[i,j] = 0.5

            if mas1 == slv2 or slv1 == mas2:
                vcm_pat[i,j] = -0.5

            if mas1 == mas2 and slv1 == slv2:
                vcm_pat[i,j] = 1.0 # handle testing ifg against itself

    # make covariance matrix in time domain
    std = sqrt(maxvar).reshape((nifgs,1))
    vcm_t = std * std.transpose()
    return vcm_t * vcm_pat
