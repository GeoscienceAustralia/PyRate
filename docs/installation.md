# Installation

## Quickstart

Before you start, make sure your system has the following packages installed:

    sudo apt-get install gdal-bin libgdal-dev libpng12-dev libblas-dev liblapack-dev libatlas-dev libatlas-base-dev gfortran libproj-dev openmpi-bin libopenmpi-dev netcdf-bin libnetcdf11 libnetcdf-dev


We strongly recommend using a [`virtualenv`](https://gist.github.com/basaks/b33ea9106c7d1d72ac3a79fdcea430eb).

One might need to install a slightly older `python-daemon` to install `PyRate`:

    pip install python-daemon==2.1.1

To install, simply run ``setup.py``:

    python setup.py install

or install with ``pip``:

    pip install git+https://github.com/GeoscienceAustralia/PyRate

The python requirements should automatically be built and installed.


## Anaconda setup for `PyRate`

For anaconda installation and `virtualenv` instruction, [see this guide](conda.md).

## Tests

Tests use [pytest](http://doc.pytest.org/en/latest/) and can be found in *tests*.

To run the tests, use the following command inside the `PyRate` directory:

	pip install pytest
	cd PyRate
	export PYRATEPATH=/path/to/PyRate
	pytest tests/