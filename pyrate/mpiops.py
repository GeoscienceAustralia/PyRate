# pylint: disable=no-member
# pylint: disable=invalid-name
"""
MPI convenience functions
"""
import logging
import pickle
from mpi4py import MPI
import numpy as np

log = logging.getLogger(__name__)
# We're having trouble with the MPI pickling and 64bit integers
MPI.pickle.dumps = pickle.dumps
MPI.pickle.loads = pickle.loads

comm = MPI.COMM_WORLD
"""module-level MPI 'world' object representing all connected nodes
"""

size = comm.Get_size()
"""int: the total number of nodes in the MPI world
"""

rank = comm.Get_rank()
"""int: the index (from zero) of this node in the MPI world. Also known as
the rank of the node.
"""


def run_once(f, *args, **kwargs):
    """Run a function on one node, broadcast result to all
    This function evaluates a function on a single node in the MPI world,
    then broadcasts the result of that function to every node in the world.
    Parameters
    ----------
    f : callable
        The function to be evaluated. Can take arbitrary arguments and return
        anything or nothing
    args : optional
        Other positional arguments to pass on to f
    kwargs : optional
        Other named arguments to pass on to f
    Returns
    -------
    result
        The value returned by f
    """
    if rank == 0:
        f_result = f(*args, **kwargs)
    else:
        f_result = None
    result = comm.bcast(f_result, root=0)
    return result


def array_split(arr, process=None):
    """
    Convenience function for splitting array elements across MPI processes
    Parameters
    ----------
    arr: ndarray
        1-D array
    process: int, optional
        process for which array members are required.
        If None, MPI.comm.rank is used instead.

    Returns list corresponding to array members in a a process
    """
    r = process if process else rank
    return np.array_split(arr, size)[r]
