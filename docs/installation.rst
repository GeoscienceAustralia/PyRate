Installation
===============

This is an installation guide to get PyRate running on various platforms.
Follow the instructions to install PyRate and run a small toy example.

.. include:: dependencies.rst
.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify PyRate has been successfully installed, run the workflow with the
example config file and data included in the repository:

::

    pyrate conv2tif -f input_parameters.conf
    pyrate prepifg -f input_parameters.conf
    pyrate process -f input_parameters.conf
    pyrate merge -f input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and geotiff files will be available in:

:: 
    
    ~/PyRate/out
