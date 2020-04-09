#!/bin/bash
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/2.1.6
export PYTHONPATH=/apps/gdal/3.0.2/lib64/python3.7/site-packages:$PYTHONPATH
source ~/PyRateVenv/bin/activate
cd ~/PyRate/docs
make html
echo "Use following cmd to copy file to local machine:"
cd ~
echo "scp -r `whoami`@gadi.nci.org.au:`pwd ~`/PyRate/docs/_build/html C:/preview"
