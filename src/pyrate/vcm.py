# Variance/Covariance matrix

# Notes:
# Hua has hardwired cvdcalc to only use method=1, we only need to worry about
# doing the calc spectrally.


from copy import copy
from numpy import array, where, isnan, real, imag, sum, sqrt, meshgrid, reshape
from numpy import zeros, vstack, ceil, mean, exp, ones, diag, sqrt
from numpy.linalg import norm
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.optimize import fmin

from algorithm import master_slave_ids


def pendiffexp(alphamod, cvdav):
	'''TODO exp covariance model
	TODO: returns a float of magnitude of the residuals
	alphamod - TODO
	cvdav - 2 col array of radius, variance (can be plotted as X,Y line)'''
	mx = cvdav[:,1].max()
	return norm(cvdav[:,1] - (mx * exp(-alphamod * cvdav[:,0])))


def unique_points(points):
	'''Returns unique points from a list of coordinates.
	points - sequence of (y,x) or (x,y) tuples
	'''
	return vstack([array(u) for u in set(points) ] )


def cvd(ifg):
	phase = where(isnan(ifg.phase_data), 0, ifg.phase_data)

	fft_phase = fft2(phase)
	pspec = real(fft_phase)**2 + imag(fft_phase)**2
	autocorr_grid = ifft2(pspec)

	nzc = sum(sum(phase != 0))
	autocorr_grid = fftshift(real(autocorr_grid)) / nzc
	acg = reshape(autocorr_grid, ifg.num_cells, 1)

	xx,yy = meshgrid(range(ifg.WIDTH), range(ifg.FILE_LENGTH))
	r = sqrt(((xx-ifg.X_CENTRE) * ifg.X_SIZE)**2 + ((yy-ifg.Y_CENTRE) * ifg.Y_SIZE)**2)  # radius from centre???
	r = reshape(r, ifg.num_cells) #, 1)

	# ACG is on the Y axis, R is on the X axis
	tmp = unique_points(zip(acg, r))
	acg, r = tmp[:,0], tmp[:,1]

	#print 'acg', acg
	#print "/" * 40
	#print 'r', r
	#print "/" * 40

	#	TODO: the line below is a a method to remove duplicate cells
	#r = r[:ceil(len(r)/2)] # BUG comment? only need 1st half of image (symmetry)
	# NB: old one had '+nlines', reason unknown, MG says ignore it.

	#acg = acg[:len(r)]

	# TODO: write a better way of getting rid of duplicate {r, acg} pairs - sets?
	# eg. array([x for x in set([(1,1), (2,2), (1,1)])])
	# the above shortens r by some number of cells
	#print 'r',r

	acgorig = copy(acg)
	maxvar = max(acg[:len(r)])

	# pick the smallest axis to focus on a circle around the centre point????
	#print 'ifg.X_CENTRE, ifg.Y_CENTRE=', ifg.X_CENTRE, ifg.Y_CENTRE
	#print 'ifg.X_SIZE, ifg.Y_SIZE', ifg.X_SIZE, ifg.Y_SIZE
	if (ifg.X_CENTRE * ifg.X_SIZE) < (ifg.Y_CENTRE * ifg.Y_SIZE):
		maxdist = ifg.X_CENTRE * ifg.X_SIZE
		w = ifg.X_SIZE * 2  # bin width
		#print "maxdist from X axis is", maxdist
	else:
		maxdist = ifg.Y_CENTRE * ifg.Y_SIZE
		w = ifg.Y_SIZE * 2  # bin width
		#print "maxdist from Y axis is", maxdist

	#rsub = array([e for e in r if e < maxdist]) # MG: prefers to use all the data
	rsub = r
	acg=acgorig;
	#acg = [e for e in r if e < maxdist]

	rbin = ceil(rsub / w).astype(int)  # classify values of r according to bin number

	maxbin = max(rbin) + 1
	cvdav = zeros( (2, maxbin) )

	for b in range(maxbin):
		cvdav[0,b] = b * w # distance instead of bin #
		cvdav[1,b] = mean(acg[rbin == b]) # mean variance for that bin

	# calculate best fit function maxvar*exp(-alpha*r)
	#alpha = fmin(pendiffexp, [2 / (maxbin * w)], [], cvdav)
	alphaguess = 2 / (maxbin * w)
	alpha = fmin(pendiffexp, x0=alphaguess, args=(cvdav,) )
	print "1st guess, alpha", alphaguess, alpha

	assert len(alpha) == 1
	return maxvar, alpha[0]


def get_vcmt(ifgs, maxvar):
	'''Returns the temporal variance/covariance matrix.'''

	# c=0.5 for common master or slave; c=-0.5 if master of one matches slave of another
	nifgs = len(ifgs)
	vcm_pat = zeros((nifgs, nifgs))

	dates = [ifg.MASTER for ifg in ifgs] + [ifg.SLAVE for ifg in ifgs]
	ids = master_slave_ids(dates)

	for i, ifg in enumerate(ifgs): # TODO: do we need -1?
		mas1, slv1 = ids[ifg.MASTER], ids[ifg.SLAVE]

		for j, ifg2 in enumerate(ifgs): # TODO: do we need +1?
			mas2, slv2 = ids[ifg2.MASTER], ids[ifg2.SLAVE]
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
