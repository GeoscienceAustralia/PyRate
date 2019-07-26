HPC
------

It is only applicable to the Raijin HPC system of the National Computational Infrastructure.
We don't know how applicable it is to other supercomputer systems.

Login to Raijin

::

    ssh <user_name>@raijin.nci.org.au
    git clone https://github.com/GeoscienceAustralia/PyRate.git

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


Run interactive session

::


    qsub -I -q expressbw -l walltime=2:00:00,mem=256Gb,ncpus=28,wd

Wait for interactive session to start and then run workflow

::

    cd ~
    source PyRateVM/bin/activate
    source PyRate/utils/load_modules.sh
