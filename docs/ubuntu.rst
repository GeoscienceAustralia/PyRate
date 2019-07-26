Ubuntu
------

Setup
^^^^^

Build PyRate package

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    source PyRate/utils/apt_install.sh
    python3 setup.py install

Set environment variables

::

    export C_INCLUDE_PATH=/usr/include/gdal
    export CPLUS_INCLUDE_PATH=/usr/include/gdal
