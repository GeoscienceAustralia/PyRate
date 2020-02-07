Installation
===============

This is an installation guide to get PyRate up and running on various platforms.
Follow the instructions to install PyRate and run a small toy example.

.. include:: ubuntu.rst
.. include:: windows.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify PyRate has been successfully installed, run the workflow with the
example config file and data included in the repository:

::

    pyrate conv2tif -f sample_data/input_parameters.conf
    pyrate prepifg -f sample_data/input_parameters.conf
    pyrate process -f sample_data/input_parameters.conf
    pyrate merge -f sample_data/input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and geotiff files will be available in:

:: 
    
    ~/PyRate/out
