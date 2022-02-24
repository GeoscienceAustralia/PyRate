"""
This script generates input lists for PyRate baseed on the PyGAMMA workflow for descending frame S1 data in Australia.
Small modifications are necessary if using ascending frame data
Usage: . list_creator.sh
"""

# Provide path to the Gamma folder titled with the frame name
DIR="/path/to/frame/T045D"
ls -d "$DIR"/*/*/*unw.tif > ifgs.list
# "coh" is "cc" in asending frame data
ls -d "$DIR"/*/*/*flat*coh.tif > cohfiles.list
ls -d "$DIR"/*/*/*base.par > baseline.list
# Change "VV" depending on polarisation of data.
ls -d "$DIR"/*/*/*VV*mli.par > headers.list

wc -l *.list
