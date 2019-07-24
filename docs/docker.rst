Docker
------

Can be used to run PyRate on linux, Windows and MacOS.

Setup
^^^^^

::

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    cd PyRate
    docker build -t pyrate-image:1.0 .
    docker run -it --rm -v %cd%:/PyRate pyrate-image:1.0 /bin/bash

Build PyRate package

::

    cd PyRate
    python3 setup.py install

Usage
^^^^^

Run workflow

::

    pyrate prepifg pyrate_gamma.conf
    pyrate linrate pyrate_gamma.conf -c 3 -r 4
    pyrate postprocess pyrate_gamma.conf -c 3 -r 4