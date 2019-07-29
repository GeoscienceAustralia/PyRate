Installation
===============

This is a quick guide to getting the PyRate software up and running in various platforms.
Step by step guide to install Pyrate and run toy example.

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify PyRate has been successfully installed, run the workflow with the
included example config file and data:

::

    pyrate prepifg input_parameters.conf
    pyrate linrate input_parameters.conf
    pyrate postprocess input_parameters.conf

On Raijin and other HPC systems, you can utilise MPI to run PyRate in parallel:
    
::

    # Modify 'n' based on the number of processors available.
    mpirun -n 4 pyrate prepifg input_parameters.conf
    mpirun -n 4 pyrate linrate input_parameters.conf -c 2 -r 2
    mpirun -n 4 pyrate postprocess input_parameters.conf -c 2 -r 2

If the installation has been successful, this workflow will complete without 
errors and results will be available in:

:: 
    
    ~/PyRate/out
