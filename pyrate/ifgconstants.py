'''
Collection of constants for elements in ROI_PAC headers (\*.rsc)
+
Collection of constants used at different stages of ifg processing

Created on 14/09/2012

.. codeauthor:: Ben Davies
.. Sudipta Basak, GA
'''

# lookup keys for the fields in PyRate's custom GeoTIFF files
PYRATE_NCOLS = 'NCOLS'
PYRATE_NROWS = 'NROWS'
PYRATE_X_STEP = 'X_STEP'
PYRATE_Y_STEP = 'Y_STEP'
PYRATE_LAT = 'LAT'
PYRATE_LONG = 'LONG'
MASTER_DATE = 'MASTER_DATE'
SLAVE_DATE = 'SLAVE_DATE'
PYRATE_DATUM = 'DATUM'
PYRATE_TIME_SPAN = 'TIME_SPAN_YEAR'
PYRATE_WAVELENGTH_METRES = 'WAVELENGTH_METRES'
PYRATE_PHASE_UNITS = 'PHASE_UNITS'

PYRATE_ORBITAL_ERROR = 'ORBITAL_ERROR'
PYRATE_INSAR_PROCESSOR = 'INSAR_PROCESSOR'
PYRATE_APS_ERROR = 'APS_ERROR'
PROCESS_STEP = 'PR_TYPE'
GEOTIFF = 'GEOTIFF'
MULTILOOKED = 'MULTILOOKED'

DAYS_PER_YEAR = 365.25  # span of year, not a calendar year
