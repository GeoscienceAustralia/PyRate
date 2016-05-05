## How to run on NCI

### Prepare interferrograms: run_prepifg.py equivalent on NCI
To process unwrapped interferrograms on the NCI follow these steps:

1. Update the `pyrate_surat.conf` with the location of the data and output locations.
2. Make sure you have access to `/g/data/dg9/sudipta/site-packages`. This is where the packages that are not available in NCI lives.
3. Clone the PyRate repo:`git clone https://your_user_name@github.com/GeoscienceAustralia/PyRate.git`
It should prompt you for your `github` password. Enter it finish cloning. If you have two factor authentical setup, your `github` password is the long hash that you set up under `personal access tokens` sections under your `github` `accounts-->settings-->personal access tokens`.
4. Once you have cloned PyRate, you need update the `PYRATEPATH` variable in the `nci_script.sh` script to the location where you cloned `PyRate`:
    
    export PYRATEPATH=/g/data/dg9/sudipta/PyRate
    
5. Update the `nci_script.sh` with the name of your pyrate config file:

    - update the walltime, memory and cpu requirements
    - the `jobfs` requirement is not significant for `run_prepifg` and you can leave it unchanged: `#PBS -l walltime=00:10:00,mem=128GB,ncpus=64,jobfs=64GB`
    
     - change the `mpirun -np 64` to the desired number of mpi processes. Probably should use same number of processes as `ncpus`:
    
    `mpirun -np 128 python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_prepifg_pypar.py /g/data/dg9/sudipta/PyRate/pyrate_surat.conf`

6. Once you are happy with your configuration, just use `qsub nci_script.sh`.
    
    
    

    
    
    
    
    
