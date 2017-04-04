# Installation

## Initial Setup

Before you start, you will need to have a number of packages installed on your Linux system. These can either be installed directly onto the system, or you can use a virtual environment.


### Direct Install

Make sure your system has the following packages installed. If not, run the command:


    sudo apt-get install gdal-bin libgdal-dev libpng12-dev libblas-dev liblapack-dev libatlas-dev libatlas-base-dev gfortran libproj-dev openmpi-bin libopenmpi-dev netcdf-bin libnetcdf11 libnetcdf-dev


### Virtual Environment

You can use one of the two virtual environment options below:

##### Virtualenv

1. Install [virtualenv](https://gist.github.com/basaks/b33ea9106c7d1d72ac3a79fdcea430eb).
2. You may need to install a slightly older python-daemon to install PyRate:


    pip install python-daemon==2.1.1

3. Install PyRate by using one of the following options.

- Run ``setup.py``:


    python setup.py install

- Install the latest version of PyRate with ``pip`` from github:


    pip install git+https://github.com/GeoscienceAustralia/PyRate

- PyRate is also on ``pypi``, the Python package manager. To install, run:


    pip install Py-Rate

The Python requirements should automatically be built and installed.


#### Anaconda

Install Anaconda [using this guide](https://github.com/GeoscienceAustralia/PyRate/blob/master/conda.md).


## Tests

A suite of tests have been developed for use in testing PyRate functionality. The tests use [pytest](http://doc.pytest.org/en/latest/) and can be found
in the *tests/* directory. A small test dataset is included in the
*tests/test_data/* directory.

To run the tests, use the following command inside the top level *PyRate/* directory:


    pip install pytest
    cd PyRate
    export PYRATEPATH=/path/to/PyRate
    pytest tests/



