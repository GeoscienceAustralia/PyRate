import os
import re
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

REF_COLOR_MAP_PATH = os.path.join(PYRATEPATH, "utils", "colourmap.txt")
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

SIXTEEN_DIGIT_EPOCH_PAIR = r'\d{8}-\d{8}'
sixteen_digits_pattern = re.compile(SIXTEEN_DIGIT_EPOCH_PAIR)
TWELVE_DIGIT_EPOCH_PAIR = r'\d{6}-\d{6}'
twelve_digits_pattern = re.compile(TWELVE_DIGIT_EPOCH_PAIR)
EIGHT_DIGIT_EPOCH = r'\d{8}'
PTN = re.compile(EIGHT_DIGIT_EPOCH)  # match 8 digits for the dates
