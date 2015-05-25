# PyRate - a Rate Estimation tool for Python

This is a partial Python port of [Pirate](http://homepages.see.leeds.ac.uk/~earhw/software/pi-rate/) (Poly-Interferogram Rate And Time-Series Estimator), a MATLAB tool written by Hua Wang and others.


## Dependencies:

- [SciPy](www.scipy.org)
- [GDAL](www.gdal.org)

It may be easier to use your Linux package manager to install the system versions of the two packages above (unless virtualenv can install them).

- [python-graph](https://pypi.python.org/pypi/python-graph) (easy_install python-graph-core)
- [pyproj](https://pypi.python.org/pypi/pyproj) (pip install pyproj)
- [nose](https://pypi.python.org/pypi/nose/) (pip install nose)
- [pytest](https://pypi.python.org/pypi/pytest) (pip install pytest)
- [sphinx](http://sphinx-doc.org/) for building the docs.


## Installation

Currently this application does not require any installation. You need to add the folder *pyrate* to your PYTHONPATH.

To run the tests, you need to set the envrironment variable *PYRATEPATH* to the folder where you put the test data... which you can get from the NCI at */g/data/dg9/PyRate_Project/data/*. There are several version of the data there found in date stampled zip files... you should probably use the most recent.


## Documentation

The [Sphinx](http://sphinx-doc.org/) documenation of the scripts can be found under the main documentation tree for the project, which is found in the top level folder *doc*. To build the documentation do

	apidoc
	make html

The documentation will be generated in *doc/build*. The main entry point is *index.html*.


## Tests

Tests use *unittest* and can be found in *pyrate/tests*. Note that you need to get the test data and set an environment variable as described in [Installation](#installation).


## Basic Usage Instructions

The runnable programs can be found in *pyrate/scripts*.

PyRate and associated tools are intended to be run TODO.

The main steps are:

1. Data formatting for PyRate
1. TODO: configuration options? (edit the config file)
1. Transformations (multi-looking, resampling etc)
1. PyRate workflow

Data formatting:
The utility scripts gamma.py & roipac.py are provided to convert ROI_PAC and GAMMA data into a GeoTIFF format usable by PyRate.

For command line options/help:
python gamma.py -h
python roipac.py -h

GAMMA needs some a DEM header to extract metadata required for the formatting. A ROIPAC translation automatically looks for some header files, and may need a DEM header file for projection data.


Data transformations/multilooking:
This separate step is handled with prepifg.py.

For command line options/help:
python prepifg.py -h

The script requires the pyrate runtime configuration file. If a config file is not provided as an arg, the script will look for 'pyrate.cfg'.  


PyRate workflow:
This is the core of the processing tools, handled with pyrate.py.

For command line options/help:
python pyrate.py -h

This script also requires runtime config file. As with the previous step, if a config file is not provided as an arg, the script will look for 'pyrate.cfg' by default.


## Todos

- NCI/GA contact details
- description of what the software does
- Licensing information
