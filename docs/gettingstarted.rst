Installation
============

Before you start, you will need to have a number of packages installed
on your Linux system. These can either be installed directly onto the
system, or you can use a virtual environment.

Ubuntu
------

Setup
~~~~~
Usage
~~~~~

Docker
------

Setup
~~~~~
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
~~~~~

Run workflow

::

    pyrate prepifg pyrate_gamma.conf
    pyrate linrate pyrate_gamma.conf -c 3 -r 4
    pyrate postprocess pyrate_gamma.conf -c 3 -r 4

Build Sphinx docs

::

    cd /PyRate/docs
    make html


NCI
------

Setup
~~~~~

Login to Raijin

::

    ssh <user_name>@raijin.nci.org.au
    git clone git clone https://github.com/GeoscienceAustralia/PyRate.git

Load the required Raijin  modules and fix the PYTHONPATH

::

    source PyRate/utils/load_modules.sh

Create a python virtual environment

::


    python3 -m venv PyRateVM
    source PyRateVM/bin/activate

Install required libraries

::

    cd PyRate
    python setup.py install

Usage
~~~~~

Run interactive session

::


    qsub -I -q expressbw -l walltime=2:00:00,mem=256Gb,ncpus=28,wd

Wait for interactive session to start and then run workflow

::

    cd ~
    source PyRateVM/bin/activate
    source PyRate/utils/load_modules.sh

    cd PyRate
    mpirun -np 4 pyrate prepifg pyrate_gamma.conf
    mpirun -np 4 pyrate linrate pyrate_gamma.conf
    mpirun -np 4 pyrate postprocess pyrate_gamma.conf


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

