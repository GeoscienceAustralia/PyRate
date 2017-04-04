Running PyRate on a HPC System
==============================

This is a quick guide to getting the PyRate software up and
running in a Portable Batch System (PBS) batch environment with MPI
support. This setup is common in High Performance Compute (HPC) systems
such as the National Computational Infrastructure (NCI)'s `Raijin
<http://nci.org.au/systems-services/national-facility/peak-system/raijin/>`__
system.

The instructions below are tailored to the NCI Raijin system. Some
instructions may not be applicable depending on the setup of the HPC
system you are using. They should apply to both single-node and
multi-node jobs; just set the number of cpus in the PBS directives in
the job submission script accordingly (e.g. ncpus=32 for 2 nodes).

These instructions assume you are using bash shell.

----------------
Pre-installation
----------------

These instructions currently only work with gcc and not Intel compilers.

1. Clone the PyRate repository into your home directory, or another directory of your choice:

   .. code:: bash

       $ cd ~
       $ git clone git@github.com:GeoscienceAustralia/PyRate.git

2. Unload the icc compiler and default openmpi from the terminal:

   .. code:: bash

       $ module unload intel-cc
       $ module unload intel-fc
       $ module unload openmpi

3. Load the modules required for installation and running:

   .. code:: bash

       $ module load python3/3.4.3 python3/3.4.3-matplotlib
       $ module load hdf5/1.8.10 gdal/2.0.0 openmpi/1.8 netcdf/4.3.2

   (Alternatively, you may wish to add the above lines to your
   ``~/.profile`` file)

4. Now add the following lines to the end of your ``~/.profile`` file:

   .. code:: bash

       export PATH=$HOME/.local/bin:$PATH
       export PYTHONPATH=$HOME/.local/lib/python3.4/site-packages:$PYTHONPATH
       export PYRATEPATH=~/PyRate
       export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.4.3/bin/python3
       export LC_ALL=en_AU.UTF-8
       export LANG=en_AU.UTF-8
       source $HOME/.local/bin/virtualenvwrapper.sh

5. Install virtualenv and ``virtualenvwrapper`` from the terminal:

   .. code:: bash

       $ pip3 install  --user virtualenv virtualenvwrapper

6. Refresh your environment by reloading your ``~/.profile`` file:

   .. code:: bash

       $ source ~/.profile


------------
Installation
------------

1. Create a new virtualenv for PyRate:

   .. code:: bash

       $ mkvirtualenv --system-site-packages pyrate

2. Make sure the virtualenv is activated:

   .. code:: bash

       $ workon pyrate

3. Install ``pyrate``:

   .. code:: bash

       $ cd $PYRATEPATH
       $ pip install python-daemon==2.1.1  # the latest python-daemon had
       $ python setup.py install

4. Once installation is complete, you can run the tests to verify everything has gone correctly:

   .. code:: bash

       $ pip install pytest
       $ py.test ~/PyRate/tests/


-----------------
Updating the Code
-----------------

To update the PyRate code, first make sure you are in the ``pyrate`` virtual
environment:

.. code:: bash

    $ workon pyrate

Next, pull the latest commit from the master branch, and install:

.. code:: bash

    $ cd $PYRATEPATH
    $ git pull origin
    $ python setup.py install

If the pull and the installation complete successfully, the code is
ready to run!


------------------
Running Batch Jobs
------------------

In the ``pbs/`` subfolder of the ``PyRate`` repository there are some example
scripts to assist launching batch jobs over multiple nodes with PBS.


Batch testing
~~~~~~~~~~~~~

To check everything is working, submit the tests as a batch job:

.. code:: bash

    $ cd $PYRATEPATH/pbs
    $ qsub submit_tests.sh

MPIRun
~~~~~~

PyRate uses MPI internally for parallelization. To run a script or
demo, use the command:

.. code:: bash

    $ mpirun -n <num_procs> <command>

For example:

.. code:: bash

    $ mpirun -n 16 pyrate prepifg pyrate_pbs.conf

A PBS job submission script might look like this:

.. code:: bash

    #!/bin/bash
    #PBS -P <project>
    #PBS -q <queue>
    #PBS -l walltime=01:00:00,mem=128GB,ncpus=16,jobfs=20GB
    #PBS -l wd

    # setup environment
    module unload intel-cc
    module unload intel-fc
    module load python3/3.4.3 python3/3.4.3-matplotlib 
    module load load hdf5/1.8.10 gdal/2.0.0
    source $HOME/.profile

    # start the virtualenv
    workon pyrate

    # run PyRate commands
    mpirun -n 16 pyrate prepifg /path/to/config_file.conf
    mpirun -n 16 pyrate linrate /path/to/config_file.conf
    mpirun -n 16 pyrate postprocess /path/to/config_file.conf
