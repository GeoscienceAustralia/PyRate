# Running PyRate on the NCI

This README is a quick guide to getting `PyRate` up and running
in a PBS batch environment that has MPI support. This setup is common in
HPC systems such as the NCI (raijin).

The instructions below should apply to both single- and multi-node runs
on the NCI. Just set ncpus in the PBS  directives in the job submission
script accordingly (e.g. ncpus=32 for 2 nodes).

The instructions assume you are using bash shell.

## Pre-installation

These instructions currently only work with gcc and not the Intel compiler.
Note that on NCI it appears python is compiled against gcc anyway.


1. Clone the `PyRate` repo into your home directory:
```bash
$ cd ~
$ git clone git@github.com:GeoscienceAustralia/PyRate.git
```

2. Unload the icc compiler and default openmpi from the terminal:
```bash
$ module unload intel-cc
$ module unload intel-fc
$ module unload openmpi
```

3. Load the modules required for installation and running:
```bash
$ module load python3/3.4.3 python3/3.4.3-matplotlib 
$ module load hdf5/1.8.10 gdal/2.0.0 openmpi/1.8 netcdf/4.3.2
```
(Alternatively, you may wish to add the above lines to your ~/.profile)

4. Now add the following lines to the end of your ~/.profile:
```bash
export PATH=$HOME/.local/bin:$PATH
export PYTHONPATH=$HOME/.local/lib/python3.4/site-packages:$PYTHONPATH
export PYRATEPATH=~/PyRate
export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.4.3/bin/python3                 
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8
source $HOME/.local/bin/virtualenvwrapper.sh
``` 

5. Install virtualenv and `virtualenvwrapper` by running the following command
on the terminal:
```bash
$ pip3 install  --user virtualenv virtualenvwrapper
```

6. Refresh your environment by reloading your profile:
```bash
$ source ~/.profile
```

## Installation

1. Create a new virtualenv for pyrate:
```bash
$ mkvirtualenv --system-site-packages pyrate
```

2. Make sure the virtualenv is activated:
```bash
$ workon pyrate
```

3. Install `pyrate`:
```bash
$ cd $PYRATEPATH
$ pip install python-daemon==2.1.1  # the latest python-daemon had 
$ python setup.py install
```

4. Once installation has completed, you can run the tests to verify everything
has gone correctly:
```bash
$ pip install pytest
$ py.test ~/PyRate/tests/
```

## Updating the Code
To update the code, first make sure you are in the `pyrate` virtual environment:
```bash
$ workon pyrate
```
Next, pull the latest commit from the master branch, and install:
```bash
$ cd $PYRATEPATH
$ git pull origin
$ python setup.py install
```
If the pull and the installation complete successfully, the code is ready to run!

## Running Batch Jobs

in the `pbs` subfolder of `PyRate` there are some example scripts and a
helper function to assist launching batch jobs over multiple nodes with PBS

### Batch testing

To check everything is working, submit the tests as a batch job:
```bash
$ cd $PYRATEPATH/pbs
$ qsub submit_tests.sh
```

### MPIRun

`PyRate` uses MPI internally for parallelization. To run a script or demo
simply do

```bash
$ mpirun -n <num_procs> <command>
```

whilst a PBS job submission might look like this:

```bash
#!/bin/bash
#PBS -P ge3
#PBS -q normal
#PBS -l walltime=01:00:00,mem=128GB,ncpus=32,jobfs=20GB
#PBS -l wd

# setup environment
module unload intel-cc
module unload intel-fc
module load python3/3.4.3 python3/3.4.3-matplotlib 
module load load hdf5/1.8.10 gdal/2.0.0
source $HOME/.profile

# start the virtualenv
workon pyrate

# run command
mpirun --mca mpi_warn_on_fork 0 pyrate prepifg /path/to/config_file.txt
mpirun --mca mpi_warn_on_fork 0 pyrate linrate /path/to/config_file.txt
```

