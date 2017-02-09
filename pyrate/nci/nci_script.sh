#!/bin/bash
#PBS -P dg9
#PBS -q normal
#PBS -l walltime=00:10:00,mem=128GB,ncpus=64,jobfs=64GB
#PBS -l wd
module rm intel-fc intel-cc
module load intel-mkl
module load python/2.7.6
module load python/2.7.6-matplotlib
module rm openmpi/1.6.3
module load openmpi/1.8
module load gdal/1.11.1-python
module load pypar/30May-2.7.6-1.8
export GDAL_DISABLE_READDIR_ON_OPEN=YES
export PYRATEPATH=/g/data/dg9/sudipta/PyRate
PYTHONPATH=$PYRATEPATH:$PYTHONPATH
PYTHONPATH=/g/data/dg9/sudipta/site-packages:$PYTHONPATH
export PYTHONPATH
mpirun -np 64 python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_prepifg_pypar.py /g/data/dg9/sudipta/PyRate/pyrate_surat.conf