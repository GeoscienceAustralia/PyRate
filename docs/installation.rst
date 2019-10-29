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

    pyrate conv2tif -f input_parameters.conf
    pyrate prepifg -f input_parameters.conf
    pyrate process -f input_parameters.conf
    pyrate merge -f input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and results will be available in:

:: 
    
    ~/PyRate/out
