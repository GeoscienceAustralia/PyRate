#   This Python module is part of the PyRate software package.
#
#   Copyright 2017 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""
This Python module contains a collection of constants used in
various components of the PyRate software
"""

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
PYRATE_APS_ERROR = 'APS_ERROR'
MULTILOOKED = 'MULTILOOKED_IFG'
ORIG = 'ORIGINAL_IFG'
DEM = 'DEM'
INCIDENCE = 'INCIDENCE_ANGLE_MAP'
INCR = 'INCREMENTAL_TIME_SLICE'
CUML = 'CUMULATIVE_TIME_SLICE'
LINRATE = 'LINEAR_RATE_MAP'
LINERROR = 'LINEAR_RATE_ERROR_MAP'
LINSAMP = 'LINEAR_RATE_SAMPLES'
PYRATE_ORBITAL_ERROR = 'ORBITAL_ERROR'
ORB_REMOVED = 'REMOVED'
REF_PHASE = 'REFERENCE_PHASE'
REF_PHASE_REMOVED = 'REMOVED'
NAN_STATUS = 'NAN_STATUS'
NAN_CONVERTED = 'CONVERTED'
DATA_TYPE = 'DATA_TYPE'
DATA_UNITS = 'DATA_UNITS'

DAYS_PER_YEAR = 365.25  # span of year, not a calendar year
SPEED_OF_LIGHT_METRES_PER_SECOND = 3e8
MM_PER_METRE = 1000
