'''
This is used to create the dummy .lv_theta file for sydney_test_gamma
'''

import os
import numpy as np
from pyrate import shared
from pyrate.tests import common
from pyrate import gamma
from pyrate import config as cf
from pyrate import ifgconstants as ifc
from osgeo import gdal

lv_theta_file = os.path.join(common.SYD_TEST_GAMMA,
                             os.path.splitext(common.SYD_TEST_DEM_GAMMA)[0]
                             + '.lv_theta')
dest = os.path.splitext(lv_theta_file)[0] + '_lv_theta.tif'
dem_header_file = common.SYD_TEST_DEM_HDR_GAMMA

dem_header = gamma.parse_dem_header(dem_header_file)

header = gamma.parse_epoch_header(
    os.path.join(common.SYD_TEST_GAMMA, '20060828_slc.par'))


incidence_angle = header[ifc.INCIDENCE_ANGLE]
data = np.ones(shape=(dem_header[ifc.PYRATE_NROWS],
                      dem_header[ifc.PYRATE_NCOLS])
                      ) * incidence_angle

shared.write_unw_from_data_or_geotiff(geotif_or_data=data,
                                      dest_unw=lv_theta_file,
                                      ifg_proc=1)

header.update(dem_header)
header[ifc.PYRATE_TIME_SPAN] = 0
header[ifc.SLAVE_DATE] = 0
header[ifc.PYRATE_PHASE_UNITS] = 'degrees'
header[ifc.PROCESS_STEP] = ifc.GEOTIFF
shared.write_geotiff(header=header, data_path=lv_theta_file,
                     dest=dest, nodata=np.nan)


ds = gdal.Open(dest, gdal.GA_ReadOnly)
data_lv_theta = ds.ReadAsArray()
ds = None
np.testing.assert_array_almost_equal(data, data_lv_theta)

