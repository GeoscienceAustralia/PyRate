#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
import numpy as np


"""For MPI compatibility"""
try:
    from mpi4py import MPI
    MPI_INSTALLED = True
    # # We're having trouble with the MPI pickling and 64bit integers
    # MPI.pickle.__init__(pickle.dumps, pickle.loads)

    # module-level MPI 'world' object representing all connected nodes
    comm = MPI.COMM_WORLD

    # int: the total number of nodes in the MPI world
    size = comm.Get_size()

    # int: the index (from zero) of this node in the MPI world. Also known as
    # the rank of the node.
    rank = comm.Get_rank()
except ImportError:
    MPI_INSTALLED = False
    size = 1
    rank = 0

    class MPI:
        SUM = np.sum

    class comm:
        """
        the mpi simulators that are used in a non-mpi environment
        """

        @staticmethod
        def barrier():
            pass

        @staticmethod
        def reduce(arr, op, root=0):
            return op(arr)

        @staticmethod
        def Get_size():
            return 1

        @staticmethod
        def allgather(* args):
            return args

        @staticmethod
        def gather(* args, **kwargs):
            return args

        @staticmethod
        def allreduce(arr, op):
            return op(arr)

        @staticmethod
        def Bcast(arr, root=0):
            return

        @staticmethod
        def bcast(arr, root=0):
            return arr


class MPIException(Exception):
    pass


def validate_mpi():
    if not MPI_INSTALLED:
        raise MPIException("MPI needs to be installed in order to use this module")


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
    if MPI_INSTALLED:
        if rank == 0:
            # without the try except, any error in this step hangs the other mpi threads
            # and these processes never exit, i.e., this MPI call hangs in case process 0 errors.
            try:
                f_result = f(*args, **kwargs)
            except Exception as e:
                f_result = e
        else:
            f_result = None
        result = comm.bcast(f_result, root=0)
        if isinstance(result, Exception):
            raise result
        else:
            return result
    else:
        return f(*args, **kwargs)


def array_split(arr: Iterable, process: int = None) -> Iterable:
    """
    Convenience function for splitting array elements across MPI processes

    :param ndarray arr: Numpy array
    :param int process: Process for which array members are required.
                If None, MPI.comm.rank is used instead. (optional)

    :return List corresponding to array members in a process.
    :rtype: list
    """
    if MPI_INSTALLED:
        r = process if process else rank
        return np.array_split(np.array(arr, dtype=object), size)[r]
    else:
        return np.array(arr)


def sum_vars(x, y, dtype):
    s = np.sum([x, y], axis=0)
    return s


def sum_axis_0(x, y, dtype):
    s = np.sum(np.stack((x, y)), axis=0)
    return s


if MPI_INSTALLED:
    sum_op = MPI.Op.Create(sum_vars, commute=True)
    sum0_op = MPI.Op.Create(sum_axis_0, commute=True)
else:
    sum_op = lambda arr: np.sum(np.expand_dims(arr, 0), axis=0)
    sum0_op = lambda arr: np.sum(np.stack(np.expand_dims(arr, 0), axis=0), axis=0)
