Instructions for setting up up an anaconda virtualenv on the linux server 'rhe-compute1' at GA.
 
## Download the conda installer

For `rhe-compute` [download this version] (https://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86.sh). 

## Install anaconda:

On 64bit linux:
    
    bash Anaconda2-4.1.1-Linux-x86.sh

Accept all default options. This will install anaconda in the `~/anaconda2` directory.
 
## Set up proxys for pip (only applicable if working behind proxy):
 
    export http_proxy=http://proxy.ga.gov.au:8080
    export https_proxy=http://proxy.ga.gov.au:8080
 
## Set up `~/.condarc` with proxy pass information (only applicable if working behind proxy):
 
For conda installer to work you need to create a `~/.condarc` file with the following proxy-pass information:
 
    proxy_servers:
        http: http://proxy.ga.gov.au:8080
        https: http://proxy.ga.gov.au:8080
        
## Clone the repo:

The rest of the instructions assume that you cloned the `PyRate` Git repository in your home directory. 

Clone the repo:    

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    
This will prompt for your rsa key passphrase. If you have access to `PyRate` repo and once you have entered the passphrase, `PyRate` will clone in your current directory. 


## Set up PYRATEPATH
You need to add the directory containing the `pyrate` folder to the `PYTHONPATH` environment variable.

The environment variable `PYRATEPATH` needs to point to the folder where you put the test data. This is how I set up my environment variable:

	export PYRATEPATH="/nas/users/u92456/unix/PyRate"
	export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH
    

## Activate `anaconda` and install `conda-env`

    source ~/anaconda2/bin/activate ~/anaconda2/
    conda install -c conda conda-env        

The `conda-env` package enables `yml` based installation.

Note that using the `source activate` command only works cleanly when using the bash shell

## Create the `pyrate` environment
    
    conda env create -f /path/to/PyRate/environment.yml

## Activate the `pyrate` environment

    source ~/anaconda2/bin/activate pyrate

## Support for `pygrib`, a package used by `PyAPS` (skip this step if not using PyAPS)

`pygrib` does not install properly in `redhat` systems, `rhe-compute1` being one of them. The way around is the following:

    cp ~/anaconda2/envs/pyrate/lib/libpng16.so.16 ~/anaconda2/envs/pyrate/lib/libpng16.so.16.bk
    conda install libpng=1.2.50
    mv ~/anaconda2/envs/pyrate/lib/libpng16.so.16.bk ~/anaconda2/envs/pyrate/lib/libpng16.so.16  

Explanation of the previous three steps:
    
1. The first copy command makes a temporary copy of the `.so`. This file is used by `matplotlib`, so we need to keep it. 

2. Then the second command  installs a much older version of `libpng`. This is used by `pygrib` on redhat systems. 

3. The third command just replaces the `.so` back so that `matplotlib` can find it. 
    
## Run `PyRate` tests

To run the tests, use the following command inside the `PyRate` directory:
		
	cd PyRate
	py.test tests/
	
## Deactivate
Once you are done using `PyRate` you could deactivate from the conda env using: 

    source deactivate

## Back to main anaconda if you need to:
    
    source activate /home/user/anaconda2