Installation
============

This is an installation guide to get `PyRate` running on various platforms.
Follow the instructions to install `PyRate` and run a small toy example using
the test dataset included in the repository.

.. include:: dependencies.rst

Platforms
---------

PyPI
^^^^

`PyRate` and its Python dependencies can be installed directly from the
Python Package Index (PyPI_)::

    pip install Py-Rate

.. _PyPI: https://pypi.org/project/Py-Rate/

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To verify `PyRate` has been successfully installed, run the full processing
workflow with the example config file and the small dataset included
in the repository. If you compiled the ``pyrate`` executable program::

    pyrate workflow -f input_parameters.conf

If you installed from the Python Package Index (PyPI_)::

    python3 pyrate/main.py workflow -f input_parameters.conf

If the installation has been successful, this workflow will complete without 
errors and geotiff files will be available in::

    ~/PyRate/out
