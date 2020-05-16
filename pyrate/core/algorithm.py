#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
This Python module contains a collection of generic algorithms used in PyRate
"""
from typing import Union, Iterable, Dict, Tuple
from numpy import sin, cos, unique, histogram, diag, dot
from scipy.linalg import qr, solve, lstsq
from pyrate.core.shared import EpochList, IfgException, PrereadIfg
from pyrate.core.ifgconstants import DAYS_PER_YEAR


def is_square(arr):
    """
    Determines whether an array is square or not.

    :param ndarray arr: numpy array

    :return: condition
    :rtype: bool
    """

    shape = arr.shape
    if len(shape) == 2 and (shape[0] == shape[1]):
        return True
    return False


def least_squares_covariance(A, b, v):
    """
    Least squares solution in the presence of known covariance.

    :param ndarray A: Design matrix
    :param ndarray b: Observations (vector of phase values)
    :param ndarray v: Covariances (vector of weights)

    :return: solution
    :rtype: ndarray
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

    :param ndarray phase_data: Phase band data array (eg. ifg.phase_data)
    :param tuple unit_vec: 3 component unit vector e.g. (EW, NS, vertical)

    :return: converted_phase
    :rtype: ndarray
    """

    # NB: currently not tested as implementation is too simple
    return phase_data * unit_vec


def unit_vector(incidence, azimuth):
    """
    Returns unit vector tuple (east_west, north_south, vertical).

    :param float incidence: incidence angle w.r.t. nadir
    :param float azimuth: azimuth of looking vector

    :return: Unit vector (EW, NS, vertical).
    :rtype: tuple
    """

    vertical = cos(incidence)
    north_south = sin(incidence) * cos(azimuth)
    east_west = sin(incidence) * sin(azimuth)
    return east_west, north_south, vertical


def ifg_date_lookup(ifgs, date_pair):
    """
    Returns the Interferogram which has the master and slave dates given
    in 'date_pair'.

    :param list ifgs: List of interferogram objects to search
    :param tuple date_pair: A (datetime.date, datetime.date)

    :return: interferogram list
    :rtype: list
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
    Returns the Interferogram index which has the master and slave dates
    given in 'date_pair'.

    :param list ifgs: List of interferogram objects to search
    :param tuple date_pair: A (datetime.date, datetime.date)

    :return: interferogram index
    :rtype: int
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

    raise ValueError("Cannot find Ifg with master/slave of %s" % str(date_pair))


def get_epochs(ifgs: Union[Iterable, Dict]) -> Tuple[EpochList, int]:
    """
    Returns an EpochList derived from the given interferograms.

    :param ifgs: List of interferogram objects

    :return: EpochList
    :rtype: list
    """

    if isinstance(ifgs, dict):
        ifgs = [v for v in ifgs.values() if isinstance(v, PrereadIfg)]
    combined = get_all_epochs(ifgs)
    dates, n = unique(combined, False, True)
    repeat, _ = histogram(n, bins=len(set(n)))

    # absolute span for each date from the zero/start point
    span = [(dates[i] - dates[0]).days / DAYS_PER_YEAR for i in range(len(dates))]
    return EpochList(dates, repeat, span), n


def get_all_epochs(ifgs):
    """
    Returns a sequence of all master and slave dates in given interferograms.

    :param list ifgs: List of interferogram objects

    :return: master and slave dates
    :rtype: list
    """

    return [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]


def master_slave_ids(dates):
    """
    Returns a dictionary of 'date:unique ID' for each date in 'dates'.
    IDs are ordered from oldest to newest, starting at 0.

    :param list dates: List of dates

    :return: unique dates IDs
    :rtype: dict
    """

    dset = sorted(set(dates))
    return dict([(date_, i) for i, date_ in enumerate(dset)])


def factorise_integer(n, memo={}, left=2):
    """
    Returns two factors a and b of a supplied number n such that a * b = n.
    The two factors are evaluated to be as close to each other in size as possible

    :param int n: Number to factorise
    :param dict memo: dictionary of candidate factors
    :param int left: operation flag (default = 2)

    :return: a, factor one
    :rtype: int
    :return: b, factor two
    :rtype: int
    """
    n = int(n)
    if (n, left) in memo:
        return memo[(n, left)]
    if left == 1:
        return n, [n]
    i = 2
    best = n
    bestTuple = [n]
    while i * i <= n:
        if n % i == 0:
            rem = factorise_integer(n / i, memo, left - 1)
            if rem[0] + i < best:
                best = rem[0] + i
                bestTuple = [i] + rem[1]
        i += 1

    # handle edge case when only one processor is available
    if bestTuple == [4]:
        return 2, 2

    if len(bestTuple) == 1:
        bestTuple.append(1)

    return int(bestTuple[0]), int(bestTuple[1])
