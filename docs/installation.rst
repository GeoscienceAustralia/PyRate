Installation
===============

This is a quick guide to getting the PyRate software up and running in various platforms.
Step by step guide to install Pyrate and run a toy example.

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify PyRate has been successfully installed, run the workflow with the
included example config file and data:

::
    pyrate conv2tif -f input_parameters.conf
    pyrate prepifg -f input_parameters.conf
    pyrate process -f input_parameters.conf
    pyrate merge -f input_parameters.conf

On Raijin and other HPC systems, you can utilise MPI to run PyRate in parallel:

::

    # Modify 'n' based on the number of processors available.
    mpirun -n 4 pyrate conv2tif -f input_parameters.conf
    mpirun -n 4 pyrate prepifg -f input_parameters.conf
    mpirun -n 4 pyrate process -f input_parameters.conf -c 2 -r 2
    mpirun -n 4 pyrate merge -f input_parameters.conf -c 2 -r 2

If the installation has been successful, this workflow will complete without 
errors and results will be available in:

:: 
    
    ~/PyRate/out
