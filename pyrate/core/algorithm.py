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
from numpy import unique, histogram

from core.ifgconstants import DAYS_PER_YEAR
from core.shared import EpochList, IfgException, PrereadIfg


def ifg_date_index_lookup(ifgs, date_pair):
    """Returns the Interferogram index which has the master and slave dates
    given in 'date_pair'.

    Args:
      ifgs(list): List of interferogram objects to search
      date_pair(tuple): A (datetime.date, datetime.date)

    Returns:
      int: interferogram index

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

    raise ValueError("Cannot find Ifg with " "master/slave of %s" % str(date_pair))


def get_epochs(ifgs):
    """Returns an EpochList derived from the given interferograms.

    Args:
      ifgs(list): List of interferogram objects

    Returns:
      list: EpochList

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
    """Returns a sequence of all master and slave dates in given interferograms.

    Args:
      ifgs(list): List of interferogram objects

    Returns:
      list: master and slave dates

    """

    return [ifg.master for ifg in ifgs] + [ifg.slave for ifg in ifgs]


def master_slave_ids(dates):
    """Returns a dictionary of 'date:unique ID' for each date in 'dates'. IDs
    are ordered from oldest to newest, starting at 0.

    Args:
      dates(list): List of dates

    Returns:
      dict: unique dates IDs

    """

    dset = sorted(set(dates))
    return dict([(date_, i) for i, date_ in enumerate(dset)])
