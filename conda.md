Instructions for setting up an Anaconda virtual environment on a typical Red Hat linux server.
 
## Install Anaconda:

For 64bit linux download this version of the conda installer: https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86.sh

Run the installer using the linux bash interpreter:
    
    bash Anaconda2-4.1.1-Linux-x86.sh

Accept all default options. This will install anaconda in the `~/anaconda2` directory.
 
## Set up proxys for pip (only applicable if working behind proxy):
 
    export http_proxy=<proxy_server>:<port>
    export https_proxy=<proxy_server>:<port>
 
## Set up `~/.condarc` with proxy pass information (only applicable if working behind proxy):
 
For conda installer to work you need to create a `~/.condarc` file with the following proxy-pass information:
 
    proxy_servers:
        http: <proxy_server>:<port>
        https: <proxy_server>:<port>
        
## Clone the PyRate repo:

The rest of the instructions assume that you cloned the `PyRate` Git repository in your home directory: 

    cd ~
    git clone git@github.com:GeoscienceAustralia/PyRate.git
    
This will prompt for your rsa key passphrase. If you have access to `PyRate` repo and once you have entered the passphrase, `PyRate` will clone in your current (home) directory. 


## Set up envrionment variables

You need to add the directory containing the `pyrate/` directory to the `PYTHONPATH` environment variable, and the environment variable `PYRATEPATH` needs to point to the directory containing the `tests/` directory:

	export PYRATEPATH="~/PyRate"
	export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH
    

## Activate `anaconda` and install `conda-env`

    source ~/anaconda2/bin/activate ~/anaconda2/
    conda install -c conda conda-env        

The `conda-env` package enables YAML based installation.

Note that using the `source activate` command only works cleanly when using the bash shell

## Create the `pyrate` environment

The required packages for PyRate are listed in the environment.yml (YAML format) file:
    
    conda env create -f $PYRATEPATH/environment.yml

## Activate the `pyrate` environment

    source ~/anaconda2/bin/activate pyrate

## Support for `pygrib`, a package used by `PyAPS` (skip this step if not using PyAPS)

`pygrib` does not install properly in `redhat` systems. The work around is the following:

    cp ~/anaconda2/envs/pyrate/lib/libpng16.so.16 ~/anaconda2/envs/pyrate/lib/libpng16.so.16.bk
    conda install libpng=1.2.50
    mv ~/anaconda2/envs/pyrate/lib/libpng16.so.16.bk ~/anaconda2/envs/pyrate/lib/libpng16.so.16  

Explanation of the previous three steps:
    
1. The first copy command makes a temporary copy of the `.so`. This file is used by `matplotlib`, so we need to keep it. 

2. Then the second command  installs a much older version of `libpng`. This is used by `pygrib` on redhat systems. 

3. The third command just replaces the `.so` back so that `matplotlib` can find it. 
    
## Run `PyRate` tests

To run the full suite of tests, use the following command inside the `PyRate` directory:
		
	cd ~/PyRate
	py.test tests/

To run the tests for a particular module:

	py.test tests/test_orbital.py

## Deactivate

Once you are done using `PyRate` you could deactivate the conda environment using: 

    source deactivate

## Back to main Anaconda environment if you need to:
    
    source activate /home/user/anaconda2
