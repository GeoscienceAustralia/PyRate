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
import pickle

import numpy as np
from mpi4py import MPI

# We're having trouble with the MPI pickling and 64bit integers
MPI.pickle.__init__(pickle.dumps, pickle.loads)

# module-level MPI 'world' object representing all connected nodes
comm = MPI.COMM_WORLD

# int: the total number of nodes in the MPI world
size = comm.Get_size()

# int: the index (from zero) of this node in the MPI world. Also known as
# the rank of the node.
rank = comm.Get_rank()


def run_once(f, *args, **kwargs):
    """Run a function on one node and then broadcast result to all.
    
    Args:
        f (str): The function to be evaluated. Can take arbitrary arguments and
            return anything or nothing

    Args:
      kwargs: str
      f: param
      args: 
      kwargs: 
      *args: 
      **kwargs: 

    Returns:
      unknown: The value returned by f.

    """
    if rank == 0:
        f_result = f(*args, **kwargs)
    else:
        f_result = None
    result = comm.bcast(f_result, root=0)
    return result


def array_split(arr, process=None):
    """Convenience function for splitting array elements across MPI processes
    
    :return List corresponding to array members in a process. :rtype: list
    
    Args:
        arr (ndarray): Numpy array
        process (int): Process for which array members are required. If None,
            MPI.comm.rank is used instead. (optional)
    
    Args:
      arr: param process:  (Default value = None)
      process:  (Default value = None)
    
    Returns:

    Args:
      arr: 
      process:  (Default value = None)

    Returns:

    """
    r = process if process else rank
    return np.array_split(arr, size)[r]


def chunks(jobs, size):
    """

    Args:
      jobs: param size:
      size: 

    Returns:

    """
    n = int(round(len(jobs) / size, 0))
    # handle edge case: n <<< size
    if n == 0:
        n = 1
    jobs = [jobs[i * n: (i + 1) * n] for i in range((len(jobs) + n - 1) // n)]

    for i in range(size - len(jobs)):
        jobs.append([])

    return jobs
