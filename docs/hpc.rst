HCP
------

Setup
^^^^^

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
^^^^^

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