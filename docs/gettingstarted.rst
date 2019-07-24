Installation
============

Step by step guide to install Pyrate and run toy example.

.. include:: ubuntu.rst
.. include:: docker.rst
.. include:: nci.rst
.. include:: conda.rst




Tests
-----

A suite of tests have been developed for use in testing PyRate
functionality and for further code development. The tests use
`pytest <http://doc.pytest.org/en/latest/>`__ and can be found in the
*tests/* directory. A small test dataset is included in the
*tests/test\_data/* directory.

To run the tests, use the following command inside the top level
*PyRate/* directory:

::

    cd PyRate
    pytest tests/

