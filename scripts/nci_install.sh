#!/bin/bash
rm -rf  ~/PyRateVenv
rm -rf  ~/PyRate
git clone https://github.com/GeoscienceAustralia/PyRate.git
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/4.0.2
export PYTHONPATH=/apps/gdal/3.0.2/lib64/python3.7/site-packages:$PYTHONPATH
python3 -m venv ~/PyRateVenv
source ~/PyRateVenv/bin/activate
cd ~/PyRate
python setup.py install
