Linux
^^^^^

These instructions have been tested using Ubuntu 18.04. If using another
Linux distribution, you will have to install packages equivalent to those
listed in ``PyRate/scripts/apt_install.sh``.

Clone the repository, install required packages and install `PyRate`:

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    ./PyRate/scripts/apt_install.sh

    # Add GDAL includes to C/CPLUS paths before building Python GDAL bindings.
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/gdal
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/include/gdal

    python3 setup.py install
