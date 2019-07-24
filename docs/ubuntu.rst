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
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8


Usage
^^^^^

Run workflow

::

    pyrate prepifg pyrate_gamma.conf
    pyrate linrate pyrate_gamma.conf -c 3 -r 4
    pyrate postprocess pyrate_gamma.conf -c 3 -r 4