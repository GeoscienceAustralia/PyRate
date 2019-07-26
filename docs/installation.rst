Installation
===============

This is a quick guide to getting the PyRate software up and running in various platforms.
Step by step guide to install Pyrate and run toy example.

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: hpc.rst


Verify Installation
-------------------

To run the verify PyRate has been sucessfully installed, run toy example workflow:

::

    pyrate prepifg pyrate_gamma.conf
    pyrate linrate pyrate_gamma.conf -c 3 -r 4
    pyrate postprocess pyrate_gamma.conf -c 3 -r 4

