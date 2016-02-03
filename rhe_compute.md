Here is how I set up my anaconda virtualenv on rhe-compute.
 
## Download the conda installer

For `rhe-compute` [download this version] (https://repo.continuum.io/archive/Anaconda2-2.4.1-Linux-x86_64.sh). 

## Install anaconda:

On 64bit linux:
    
    bash Anaconda2-2.4.1-Linux-x86_64.sh

Please accept all defaults. This will install anaconda in the `~/anaconda2` directory.
 
## Set up proxys for pip:
 
    export http_proxy=http://proxy.agso.gov.au:8080
    export https_proxy=http://proxy.agso.gov.au:8080
 

## Clone the repo:

The rest of the instructions assume that you cloned `PyRate` in your home directory. 

Clone the repo:    

    git clone git@github.com:GeoscienceAustralia/PyRate.git
    
This will prompt for your rsa key passphrase. If you have access to `PyRate` repo and once you have entered the passphrase, `PyRate` will clone in your current directory. 


## Set up PYRATEPATH
You need to add the directory containing the `pyrate` folder to the `PYTHONPATH` environment variable.

The environment variable `PYRATEPATH` needs to point to the folder where you put the test data. This is how I set up my environment variable:

	export PYRATEPATH="/nas/users/u92456/unix/PyRate"
	export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH
    
 
# Set up `~/.condarc` with proxy pass information:
 
For conda installer to work you need to create a `~/.condarc` file with the following proxy-pass information:
 
    proxy_servers:
        http: http://proxy.agso.gov.au:8080
        https: http://proxy.agso.gov.au:8080
        

#### Activate `anaconda` and install `conda-env`

    source ~/anaconda2/bin/activate ~/anaconda2/
    conda install -c conda conda-env        

The `conda-env` package enables `yml` based installation.

## Create the `pyrate` environment
    
    conda env create -f /path/to/PyRate/environment.yml

## Activate the `pyrate` environment

    source ~/anaconda2/bin/activate pyrate
    
## Run `PyRate` tests

To run the tests, use the following command inside the `PyRate` directory:
		
	cd PyRate
	nosetests --nologcapture
	
Or you can use `unittest discover`:

	cd PyRate
	python -m unittest discover pyrate/

## Deactivate
Once you are done using `PyRate` you could deactivate from the conda env using: 

    source deactivate

## Back to main anaconda if you need to:
    
    source activate /home/user/anaconda2
