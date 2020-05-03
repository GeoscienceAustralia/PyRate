HPC
------

These instructions have been tested with the Gadi supercomputer of the 
`National Computational Infrastructure (NCI) <https://nci.org.au/>`_. 
The process for other HPC platforms may differ. 

Login to Gadi:

::

    ssh <user_name>@gadi.nci.org.au
    cd ~
    git clone https://github.com/GeoscienceAustralia/PyRate.git

Load the required Gadi modules (this will also remove the default NCI GDAL
Python bindings so we can build and use our own):

::

    source ~/PyRate/utils/load_modules.sh

Create a Python virtual environment:

::

    python3 -m venv ~/PyRateVenv
    source ~/PyRateVenv/bin/activate

Install PyRate:

::

    cd ~/PyRate
    python setup.py install


Following this, PyRate will be available for PBS jobs. To verify the 
installation, first run an interactive session:

::

    qsub -I -q express -l walltime=01:00:00,mem=16Gb,ncpus=4,wd

Once the session has started, you will need to reactivate your virtual 
environment and reload the required modules:

::

    cd ~
    source PyRateVenv/bin/activate
    source PyRate/utils/load_modules.sh
