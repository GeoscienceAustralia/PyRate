Installation
============

This is an installation guide to get `PyRate` running on various platforms.
Follow the instructions to install `PyRate` and run a small toy example.

.. include:: dependencies.rst

Platforms
---------

PyPI
^^^^

`PyRate` and its Python dependencies can be installed directly from the Python Package Index (PyPI_):

.. _PyPI: https://pypi.org/project/Py-Rate/

::

    pip install Py-Rate

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify `PyRate` has been successfully installed, run the workflow with the
example config file and the small dataset included in the repository:

::

    pyrate conv2tif -f input_parameters.conf
    pyrate prepifg -f input_parameters.conf
    pyrate correct -f input_parameters.conf
    pyrate timeseries -f input_parameters.conf
    pyrate stack -f input_parameters.conf
    pyrate merge -f input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and geotiff files will be available in:

:: 
    
    ~/PyRate/out
