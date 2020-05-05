#!/bin/bash
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/4.0.2
export PYTHONPATH=/apps/gdal/3.0.2/lib64/python3.7/site-packages:$PYTHONPATH
source ~/PyRateVenv/bin/activate
pip install -U pytest
cd ~/PyRate
pytest tests/
