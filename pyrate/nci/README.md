## Initial set up to run on the NCI

### Generate new SSH key for GitHub access on the NCI
1. Log into NCI, navigate to your home directory.
2. Generate key: `ssh-keygen -t rsa -b 4096 -C "<github username>_github"`.
3. Add key to ssh-agent: `ssh-add /home/547/<nci_user_initials>547/.ssh/id_rsa`. 
4. Add key to GitHub:
    - copy contents of the `id_rsa.pub` file located in the `.ssh` directory of your home directory.
    - on GitHub, log in and navigate to: https://github.com/settings/keys.
    - click 'New SSH key' button.
    - add a title, such as 'nci' and paste the contents copied from the `id_rsa.pub` file.
    - click 'Add SSH key'.

### Clone PyRate and PyAPS repos
1. Navigate to your repo directory `/home/547/<nci_user_initials>547/repo/`
2. Clone PyRate: 
    - `git clone git@github.com:GeoscienceAustralia/PyRate.git` 
    - `cd PyRate`
    - create new branch: `git branch <branch_name>`
    - switch to new branch: `git checkout <branch_name>`
    - push new branch to GitHub: `git push origin <branch_name>`
3. Clone PyAPS (within `PyRate` directory):
    - `git clone git@github.com:GeoscienceAustralia/PyAPS.git`
    - `cd PyAPS`
    - create new branch: `git branch <branch_name>`
    - switch to new branch: `git checkout <branch_name>`
    - push new branch to GitHub: `git push origin <branch_name>`

### Update enivornment variables
1. Add the following lines to the `.profile` file in your home directory:
    - `export PYRATEPATH=/home/547/<nci_user_initials>547/repo/PyRate`
    - `export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH`
    - `export PYTHONPATH=/home/547/<nci_user_initials>547/repo/PyRate/PyAPS:$PYTHONPATH`
2. Add the following line to the `.bashrc` file in your home directory:
    - `alias loadpyrate='source /home/547/mcg547/dg9-apps/PYRATE/PYRATE_CONFIG'`
The `PYRATE_CONFIG` file details the modules and packages required for PyRate to work. They will automatically load when `loadpyrate` is run. The packages that are not available on the NCI are located in the directory: `/g/data/dg9/sudipta/site-packages`.
3. Activate these changes by running: `source .bashrc`



## How to run on NCI

### Prepare interferrograms: run_prepifg.py equivalent on NCI
To process unwrapped interferograms on the NCI follow these steps:

1. Load PyRate modules and packages: `loadpyrate`.
2. Copy the configuration file `pyrate_nci.conf` and rename it with the project name.
3. Update your configuration file with the location of the data and output locations.
4. Copy the PBS job file `process_pyrate_xxx` and rename it with the project name.
5. Update your PBS job file:
    - to point to your configuration file.
    - update the walltime, memory and cpu requirements.
    - the `jobfs` requirement is not significant for `run_prepifg` and you can leave it unchanged: `#PBS -l walltime=00:10:00,mem=128GB,ncpus=64,jobfs=64GB`.
    - change the `mpirun -np 64` to the desired number of mpi processes. Probably should use same number of processes as `ncpus`:
    
    `mpirun -np 128 python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_prepifg_pypar.py /g/data/dg9/sudipta/PyRate/pyrate_surat.conf`

6. Once you are happy with your configuration, just use `qsub process_pyrate_<project>`.
    
    
    
