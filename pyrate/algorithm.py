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
This Python module contains a collection of algorithms used in PyRate
"""
import logging
from numpy import sin, cos, unique, histogram, diag, dot
from scipy.linalg import qr, solve, lstsq
from pyrate.shared import EpochList, IfgException, PrereadIfg
from pyrate.ifgconstants import DAYS_PER_YEAR

log = logging.getLogger(__name__)


def is_square(arr):
    """
    Parameters
    ----------
    arr: ndarray
        numpy array
    Return
    ------
    whether arr is square or not
    """
    shape = arr.shape
    if len(shape) == 2 and (shape[0] == shape[1]):
        return True
    return False


def least_squares_covariance(A, b, v):
    """
    Least squares solution in the presence of known covariance.

    This function is known as lscov() in MATLAB.
    A: design matrix
    b: observations (vector of phase values)
    v: covariances (weights) in vector form
    """
    # pylint: disable=too-many-locals
    # X = LSCOV(A,b,V) returns the vector X that minimizes
    # (A*X-b)'*inv(V)*(A*X-b) for the case in which length(b) > length(X).
    # This is the over-determined least squares problem with covariance V.
    # The solution is found without needing to invert V which is a square
    # symmetric matrix with dimensions equal to length(b).
    #
    # The classical linear algebra solution to this problem is:
    #
    #     x = inv(A'*inv(V)*A)*A'*inv(V)*b
    #
    # Reference:
    # G. Strang, "Introduction to Applied Mathematics",
    # Wellesley-Cambridge, p. 398, 1986.
    # L. S hure 3-31-89

    # convert vectors to 2D singleton array
    if len(A.shape) != 2:
        raise ValueError('')

    m, n = A.shape
    if m <= n:
        raise ValueError('Problem must be over-determined')

    V = diag(1.0 / v.squeeze())
    q, r = qr(A)  # Orthogonal-triangular Decomposition
    efg = dot(q.T, dot(V, q))  # TODO: round it??
    g = efg[n:, n:]  # modified to 0 indexing
    cd = dot(q.T, b)  # q.T * b
    f = efg[:n, n:]  # TODO: check +1/indexing
    c = cd[:n]  # modified to 0 indexing
    d = cd[n:]  # modified to 0 indexing
    r = r[:n, :n]  # modified to 0 indexing

    func = solve if is_square(g) else lstsq
    tmp = func(g, d)
    func = solve if is_square(r) else lstsq
    return func(r, (c-f * tmp))


def los_conversion(phase_data, unit_vec):
    """
    Converts phase from line-of-sight (LOS) to horizontal/vertical components.
    phase_data - phase band data array (eg. ifg.phase_data)
    unit_vec - 3 component sequence, eg. [EW, NS, vertical]
    """

    # NB: currently not tested as implementation is too simple
    return phase_data * unit_vec


def unit_vector(incidence, azimuth):
    """
    Returns unit vector tuple (east_west, north_south, vertical).
    """

    vertical = cos(incidence)
    north_south = sin(incidence) * cos(azimuth)
    east_west = sin(incidence) * sin(azimuth)
    return east_west, north_south, vertical


def ifg_date_lookup(ifgs, date_pair):
    """
    Returns an Ifg which has a master/slave dates given in 'date_pair'.
    ifgs - list of Ifg objects to search in
    date_pair - a (datetime.date, datetime.date) tuple
    """
    if len(date_pair) != 2:
        msg = "Need (datetime.date, datetime.date) master/slave pair"
        raise IfgException(msg)

    # check master/slave dates are in order
    try:
        # TODO: Clarify: Is the comparison here for a different date?
        # Then it should be written in a more pythonic way
        # The if below is always true as long as the dates are different
        # and not in order
        if date_pair[0] > date_pair[1]:
            date_pair = date_pair[1], date_pair[0]
    except:
        raise ValueError("Bad date_pair arg to ifg_date_lookup()")

    for i in ifgs:
        if date_pair == (i.master, i.slave):
            return i

    raise ValueError("Cannot find Ifg with "
                     "master/slave of %s" % str(date_pair))


def ifg_date_index_lookup(ifgs, date_pair):
    """
    Returns an Ifg index which has a master/slave dates given in 'date_pair'.
    ifgs - list of Ifg objects to search in
    date_pair - a (datetime.date, datetime.date) tuple
    """

    if len(date_pair) != 2:
        msg = "Need (datetime.date, datetime.date) master/slave pair"
        raise IfgException(msg)

    # check master/slave dates are in order
    try:
        if date_pair[0] > date_pair[1]:
            date_pair = date_pair[1], date_pair[0]
    except:
        raise ValueError("Bad date_pair arg to ifg_date_lookup()")

    for i, _ in enumerate(ifgs):
        if date_pair == (ifgs[i].master, ifgs[i].slave):
            return i

    raise ValueError("Cannot find Ifg with "
                     "master/slave of %s" % str(date_pair))


def get_epochs(ifgs):
    """
    Returns an EpochList derived from all given interferograms.
    """
    log.info('Finding epochs')
    if isinstance(ifgs, dict):
        ifgs = [v for v in ifgs.values() if isinstance(v, PrereadIfg)]
    combined = get_all_epochs(ifgs)
    dates, n = unique(combined, False, True)
    repeat, _ = histogram(n, bins=len(set(n)))

    # absolute span for each date from the zero/start point
    span = [(dates[i] - dates[0]).days / DAYS_PER_YEAR
            for i in range(len(dates))]
    return EpochList(dates, repeat, span), n


def get_all_epochs(ifgs):
    """
    Returns sequence of all master and slave dates in given ifgs.
    """

    return [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]


def master_slave_ids(dates):
    """
    Returns dict of 'date:unique ID' for each date in 'dates'. IDs are ordered
    from oldest to newest, starting at 0.
    Replaces ifglist.mas|slvnum in Pirate.
    """

    dset = sorted(set(dates))
    return dict([(date_, i) for i, date_ in enumerate(dset)])
