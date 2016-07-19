#!/bin/bash
#PBS -P dg9
#PBS -q normal
#PBS -l walltime=12:00:00,mem=256GB,ncpus=128,jobfs=200GB
#PBS -l wd
module rm intel-fc intel-cc
module load intel-mkl
module load python/2.7.6
module load python/2.7.6-matplotlib
module rm openmpi/1.6.3
module load openmpi/1.8
module load gdal/1.11.1-python
module load pypar/30May-2.7.6-1.8
export PYRATEPATH=/g/data/dg9/sudipta/PyRate
PYTHONPATH=$PYRATEPATH:$PYTHONPATH
PYTHONPATH=/g/data/dg9/sudipta/site-packages:$PYTHONPATH
export GDAL_DISABLE_READDIR_ON_OPEN=YES
export PYTHONPATH
mpirun -np 128 python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_pyrate_pypar.py /g/data/dg9/sudipta/PyRate/pyrate.conf
mpirun -map-by ppr:2:node python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_pyrate_pypar_2.py /g/data/dg9/sudipta/PyRate/pyrate.conf
mpirun -np 128 python /g/data/dg9/sudipta/PyRate/pyrate/nci/run_pyrate_pypar_3.py /g/data/dg9/sudipta/PyRate/pyrate.conf
mpirun -map-by ppr:8:node python /g/data/dg9/sudipta/PyRate/pyrate/nci/postprocessing.py /g/data/dg9/sudipta/PyRate/pyrate.conf
