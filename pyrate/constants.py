import os
import re
from pathlib import Path
import numpy as np

PYRATEPATH = Path(__file__).parent.parent


__version__ = "0.5.0"
CLI_DESCRIPTION = """
PyRate workflow: 

    Step 1: conv2tif
    Step 2: prepifg
    Step 3: correct
    Step 4: timeseries
    Step 5: stack
    Step 6: merge

Refer to https://geoscienceaustralia.github.io/PyRate/usage.html for 
more details.
"""

from pyrate.core.mpiops import comm

NO_OF_PARALLEL_PROCESSES = comm.Get_size()

CONV2TIF = 'conv2tif'
PREPIFG = 'prepifg'
CORRECT = 'correct'
TIMESERIES = 'timeseries'
STACK = 'stack'
MERGE = 'merge'

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
GAMMA_HEADING = 'heading'
GAMMA_AZIMUTH = 'azimuth_angle'
GAMMA_RANGE_PIX = 'range_pixel_spacing'
GAMMA_RANGE_N = 'range_samples'
GAMMA_AZIMUTH_PIX = 'azimuth_pixel_spacing'
GAMMA_AZIMUTH_N = 'azimuth_lines'
GAMMA_AZIMUTH_LOOKS = 'azimuth_looks'
GAMMA_PRF = 'prf'
GAMMA_NEAR_RANGE = 'near_range_slc'
GAMMA_SAR_EARTH = 'sar_to_earth_center'
GAMMA_SEMI_MAJOR_AXIS = 'earth_semi_major_axis'
GAMMA_SEMI_MINOR_AXIS = 'earth_semi_minor_axis'
GAMMA_PRECISION_BASELINE = 'precision_baseline(TCN)'
GAMMA_PRECISION_BASELINE_RATE = 'precision_baseline_rate'
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
