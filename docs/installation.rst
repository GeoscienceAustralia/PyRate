Installation
===============

This is a quick guide to getting the PyRate software up and running in various platforms.
Step by step guide to install Pyrate and run toy example.

.. include:: ubuntu.rst
.. include:: windows.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify PyRate has been successfully installed, run the workflow with the
included example config file and data:

::

    pyrate conv2tif -f sample_data/input_parameters.conf
    pyrate prepifg -f sample_data/input_parameters.conf
    pyrate process -f sample_data/input_parameters.conf
    pyrate merge -f sample_data/input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and results will be available in:

:: 
    
    ~/PyRate/out
