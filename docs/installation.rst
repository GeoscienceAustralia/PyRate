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

    pyrate prepifg pyrate_gamma.conf
    pyrate linrate pyrate_gamma.conf -c 3 -r 4
    pyrate postprocess pyrate_gamma.conf -c 3 -r 4

On Raijin and other HPC systems, you can utilise MPI to run PyRate in parallel:
    
::

    # Modify 'n' based on the number of processors available.
    mpirun -n 4 pyrate prepifg pyrate_gamma.conf
    mpirun -n 4 pyrate linrate pyrate_gamma.conf -c 3 -r 4
    mpirun -n 4 pyrate postprocess pyrate_gamma.conf -c 3 -r 4

If the installation has been successful, this workflow will complete without 
errors and results will be available in:

:: 
    
    ~/PyRate/out
