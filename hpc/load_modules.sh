module load python3/3.6.7
module load gdal/2.2.2
module load openmpi/2.1.1
# Remove the prebuilt GDAL Python bindings that get added by the NCI module.
PYTHONPATH=`echo $PYTHONPATH | sed -e 's|/apps/gdal/.*:||'`
