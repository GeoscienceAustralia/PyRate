module unload intel-cc
module unload intel-fc
module unload openmpi
module load python3/3.7.2
module load gdal/2.2.2
module load openmpi/2.1.1
# Remove the prebuilt GDAL Python bindings that get added by the NCI module.
PYTHONPATH=`echo $PYTHONPATH | sed -e 's|/apps/gdal/.*:||'`
# Required by Click
export LC_ALL=en_AU.UTF-8
export LANG=en_AU.UTF-8
module unload git
module load git