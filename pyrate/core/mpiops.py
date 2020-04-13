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
This Python module contains MPI convenience functions for PyRate
"""
# pylint: disable=no-member
# pylint: disable=invalid-name
import logging
import pickle
from typing import Callable, Any, Iterable
from mpi4py import MPI
import numpy as np

log = logging.getLogger(__name__)
# We're having trouble with the MPI pickling and 64bit integers
MPI.pickle.__init__(pickle.dumps, pickle.loads)

# module-level MPI 'world' object representing all connected nodes
comm = MPI.COMM_WORLD

# int: the total number of nodes in the MPI world
size = comm.Get_size()

# int: the index (from zero) of this node in the MPI world. Also known as
# the rank of the node.
rank = comm.Get_rank()


def run_once(f: Callable, *args, **kwargs) -> Any:
    """
    Run a function on one node and then broadcast result to all.

    :param Callable f: The function to be evaluated. Can take arbitrary arguments
                and return anything or nothing
    :param list args: Other positional arguments to pass on to f (optional)
    :param dict kwargs: Other named arguments to pass on to f (optional)

    :return: The value returned by f.
    :rtype: unknown
    """
    if rank == 0:
        f_result = f(*args, **kwargs)
    else:
        f_result = None
    result = comm.bcast(f_result, root=0)
    return result


def array_split(arr: Iterable, process: int = None) -> np.ndarray:
    """
    Convenience function for splitting array elements across MPI processes

    :param ndarray arr: Numpy array
    :param int process: Process for which array members are required.
                If None, MPI.comm.rank is used instead. (optional)

    :return List corresponding to array members in a process.
    :rtype: list
    """
    r = process if process else rank
    return np.array_split(arr, size)[r]


def sum_axis_0(x, y, dtype):
    s = np.sum([x, y], axis=0)
    return s

sum0_op = MPI.Op.Create(sum_axis_0, commute=True)
