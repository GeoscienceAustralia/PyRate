Ubuntu
------

Clone the PyRate repository, install required packages and install PyRate:

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    ./PyRate/utils/apt_install.sh

    # Add GDAL includes to C/CPLUS paths before building Python GDAL bindings.
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/gdal
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/gdal

    python3 setup.py install

PyRate can also be installed from the Python package index:

::

    pip install Py-Rate

