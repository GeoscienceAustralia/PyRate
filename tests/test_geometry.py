import numpy as np
from os.path import join
import pyrate.core.config as cf
from pyrate.core.geometry import get_lonlat_coords, get_lonlat_coords_slow
from tests import common
from pyrate.configuration import Configuration
from subprocess import run, PIPE
from pyrate import prepifg


def test_get_lonlat_coords_vectorised(dem):
    lon, lat = get_lonlat_coords_slow(dem)
    lon_v, lat_v = get_lonlat_coords(dem)
    np.testing.assert_array_almost_equal(lon, lon_v)
    np.testing.assert_array_almost_equal(lat, lat_v)


def get_pyrate_angle(self, x0, y0, tif_file):
    """
    Get angle at particular pixel in the azimuth/incidence tif file
    """
    # get azimuth angle value of PyRate file azimuth_angle.tif
    temp = run(['gdallocationinfo', tif_file, str(x0), str(y0)], universal_newlines=True, stdout=PIPE)
    out = temp.stdout
    angle_value = float(out.split('Value:')[1].strip('\n'))

    return angle_value

class TestPyRateAngleFiles:

    @classmethod
    def setup_method(cls):
        cls.params = Configuration(common.MEXICO_CONF).__dict__
        # run prepifg
        prepifg.main(cls.params)

    @classmethod
    def teardown_method(cls):
        pass

    # the xy position in the original GAMMA geometry (before cropping and further multi-looking is required for
    # comparison of PyRate with GAMMA angle values. To get this position for a particular coordinate, the following
    # commands can be applied on Gadi:
    # cd /g/data/dg9/INSAR_ANALYSIS/MEXICO/S1/UNIT-TEST-DATA/CROP-A/GEOTIFF
    # gdallocationinfo -wgs84 cropA_20180106-20180130_VV_8rlks_eqa_unw.tif -99.06 19.37
    # -> outputs the cropped pixel coordinates: (94P,58L) -> x0 = 94, y0 = 58
    # cd /g/data/dg9/INSAR_ANALYSIS/MEXICO/S1/GAMMA/T005A/geotiff_files/unw_ifg_geotiffs/
    # gdallocationinfo -wgs84 20180106-20180130_VV_8rlks_eqa_unw.tif -99.06 19.37
    # -> outputs the original pixel coordinates: (940P,3271L)
    # cd /g/data/dg9/INSAR_ANALYSIS/MEXICO/S1/GAMMA/T005A/DEM/
    # gdallocationinfo 20180106_VV_8rlks_eqa_lv_phi.tif 940 3271
    # -> outputs the GAMMA azimuth angle: -2.94684600830078
    # gdallocationinfo 20180106_VV_8rlks_eqa_lv_theta.tif 940 3271
    # > outputs the GAMMA incidence angle: 1.0111095905304

    def test_pyrate_azimuth_matches_gamma_azimuth_x50_y29(self):
        x0 = 29; y0 = 50
        # azimuth angle validation value extracted from GAMMA azimuth file for pixel location corresponding
        # to (29,50) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_phi.tif 613 3228
        exp = -2.94634866714478 # GAMMMA azimuth angle at pixel location
        # GAMMA azimuth is defined towards satellite in an anti-clockwise angle system, with East being zero
        # PyRate azimuth angle is defined towards the satellite in a clockwise angle system with North being zero
        tif_file = join(self.params[cf.OUT_DIR], 'azimuth_angle.tif')
        azimuth_angle_pyrate = get_pyrate_angle(self, x0, y0, tif_file)
        # convert PyRate azimuth into GAMMA azimuth
        res = -(azimuth_angle_pyrate - np.pi / 2)
        np.testing.assert_array_almost_equal(exp, res, decimal=2) # max difference < 0.01 rad
        # screen output of difference if need be:
        #print(exp - res)

    def test_pyrate_incidence_matches_gamma_incidence_x50_y29(self):
        x0 = 29; y0 = 50
        # incidence angle validation value extracted from GAMMA incidence file for pixel location corresponding
        # to (29,50) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_theta.tif 613 3228
        exp = 1.02217936515808
        # GAMMA angle is defined from the horizontal plane with the zenith direction being pi / 2 radians (i.e.90 deg)
        # PyRate angle is defined from the vertical axis with the zenith direction being zero
        tif_file = join(self.params[cf.OUT_DIR], 'incidence_angle.tif')
        incidence_angle_pyrate = get_pyrate_angle(self, x0, y0, tif_file)
        # convert PyRate incidence into GAMMA incidence
        res = np.pi / 2 - incidence_angle_pyrate
        # convert PyRate incidence into GAMMA incidence
        np.testing.assert_array_almost_equal(exp, res, decimal=3) # max difference < 0.001 rad
        # screen output of difference if need be:
        #print(exp - res)

    def test_pyrate_azimuth_matches_gamma_azimuth_x94_y58(self):
        x0 = 94; y0 = 58
        # azimuth angle validation value extracted from GAMMA azimuth file for pixel location corresponding
        # to (94,58) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_phi.tif 940 3271
        exp = -2.94684600830078 # GAMMMA azimuth angle at pixel location
        # GAMMA azimuth is defined towards satellite in an anti-clockwise angle system, with East being zero
        # PyRate azimuth angle is defined towards the satellite in a clockwise angle system with North being zero
        tif_file = join(self.params[cf.OUT_DIR], 'azimuth_angle.tif')
        azimuth_angle_pyrate = get_pyrate_angle(self, x0, y0, tif_file)
        # convert PyRate azimuth into GAMMA azimuth
        res = -(azimuth_angle_pyrate - np.pi / 2)
        np.testing.assert_array_almost_equal(exp, res, decimal=2) # max difference < 0.01 rad
        # screen output of difference if need be:
        #print(exp - res)

    def test_pyrate_incidence_matches_gamma_incidence_x94_y58(self):
        x0 = 94; y0 = 58
        # incidence angle validation value extracted from GAMMA incidence file for pixel location corresponding
        # to (94,58) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_theta.tif 940 3271
        exp = 1.0111095905304
        # GAMMA angle is defined from the horizontal plane with the zenith direction being pi / 2 radians (i.e.90 deg)
        # PyRate angle is defined from the vertical axis with the zenith direction being zero
        tif_file = join(self.params[cf.OUT_DIR], 'incidence_angle.tif')
        incidence_angle_pyrate = get_pyrate_angle(self, x0, y0, tif_file)
        # convert PyRate incidence into GAMMA incidence
        res = np.pi / 2 - incidence_angle_pyrate
        # convert PyRate incidence into GAMMA incidence
        np.testing.assert_array_almost_equal(exp, res, decimal=3) # max difference < 0.001 rad
        # screen output of difference if need be:
        #print(exp - res)
