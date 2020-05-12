#!bin/bash
module purge
module load python3/3.7.4
module load gdal/3.0.2
module load openmpi/4.0.2
# Remove the prebuilt GDAL Python bindings that get added by the NCI module.
export PYTHONPATH=/apps/gdal/3.0.2/lib64/python3.7/site-packages:$PYTHONPATH
# Required by Click
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8
