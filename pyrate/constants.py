import os
from pathlib import Path
import numpy as np

PYRATEPATH = Path(__file__).parent.parent


__version__ = "0.4.0"
CLI_DESCRIPTION = """
PyRate workflow: 

    Step 1: conv2tif
    Step 2: prepifg
    Step 3: process
    Step 4: merge 

Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for 
more details.
"""
from mpi4py import MPI
comm = MPI.COMM_WORLD
NO_OF_PARALLEL_PROCESSES = comm.Get_size()

CONV2TIF = 'conv2tif'
PREPIFG = 'prepifg'
PROCESS = 'process'
MERGE = 'merge'

REF_COLOR_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils", "colormap.txt")
# distance division factor of 1000 converts to km and is needed to match legacy output
DISTFACT = 1000
# mappings for metadata in header for interferogram
GAMMA_DATE = 'date'
GAMMA_TIME = 'center_time'
GAMMA_WIDTH = 'width'
GAMMA_NROWS = 'nlines'
GAMMA_CORNER_LAT = 'corner_lat'
GAMMA_CORNER_LONG = 'corner_lon'
GAMMA_Y_STEP = 'post_lat'
GAMMA_X_STEP = 'post_lon'
GAMMA_DATUM = 'ellipsoid_name'
GAMMA_FREQUENCY = 'radar_frequency'
GAMMA_INCIDENCE = 'incidence_angle'
# RADIANS = 'RADIANS'
# GAMMA = 'GAMMA'
# value assigned to no-data-value
LOW_FLOAT32 = np.finfo(np.float32).min*1e-10

# lookup keys for the metadata fields in PyRate GeoTIFF files
PYRATE_NCOLS = 'NCOLS'
PYRATE_NROWS = 'NROWS'
PYRATE_X_STEP = 'X_STEP'
PYRATE_Y_STEP = 'Y_STEP'
PYRATE_LAT = 'LAT'
PYRATE_LONG = 'LONG'
MASTER_DATE = 'MASTER_DATE'
MASTER_TIME = 'MASTER_TIME'
SLAVE_DATE = 'SLAVE_DATE'
SLAVE_TIME = 'SLAVE_TIME'
EPOCH_DATE = 'EPOCH_DATE'
PYRATE_DATUM = 'DATUM'
PYRATE_TIME_SPAN = 'TIME_SPAN_YEAR'
PYRATE_WAVELENGTH_METRES = 'WAVELENGTH_METRES'
PYRATE_INCIDENCE_DEGREES = 'INCIDENCE_DEGREES'
PYRATE_INSAR_PROCESSOR = 'INSAR_PROCESSOR'
PYRATE_WEATHER_ERROR = 'WEATHER_ERROR'
PYRATE_APS_ERROR = 'APS_ERROR'
PYRATE_MAXVAR = 'CVD_MAXVAR'
PYRATE_ALPHA = 'CVD_ALPHA'
COHERENCE = 'COHERENCE_MASKED_MULTILOOKED_IFG'
MULTILOOKED = 'MULTILOOKED_IFG'
ORIG = 'ORIGINAL_IFG'
DEM = 'ORIGINAL_DEM'
MLOOKED_DEM = 'MULTILOOKED_DEM'
INCIDENCE = 'INCIDENCE_ANGLE_MAP'
MLOOKED_INC = 'MULTILOOKED_INCIDENCE_ANGLE_MAP'
INCR = 'INCREMENTAL_TIME_SLICE'
CUML = 'CUMULATIVE_TIME_SLICE'
STACKRATE = 'STACKED_RATE_MAP'
STACKERROR = 'STACKED_RATE_ERROR'
STACKSAMP = 'STACKED_RATE_SAMPLES'
PYRATE_ORBITAL_ERROR = 'ORBITAL_ERROR'
ORB_REMOVED = 'REMOVED'
APS_REMOVED = 'REMOVED'
PYRATE_REF_PHASE = 'REFERENCE_PHASE'
REF_PHASE_REMOVED = 'REMOVED'
NAN_STATUS = 'NAN_STATUS'
NAN_CONVERTED = 'CONVERTED'
DATA_TYPE = 'DATA_TYPE'
DATA_UNITS = 'DATA_UNITS'
