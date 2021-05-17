import shutil
from typing import Tuple
import numpy as np
from os.path import join
import pytest

import pyrate.constants as C
from pyrate.core import ifgconstants as ifc
from pyrate.core.geometry import get_lonlat_coords, get_sat_positions, vincinv
from pyrate.core.refpixel import convert_pixel_value_to_geographic_coordinate
from tests import common
from pyrate.configuration import Configuration
from subprocess import run, PIPE
from pyrate import prepifg, correct
from pyrate.core.shared import Ifg, Geometry


def get_lonlat_coords_slow(ifg: Ifg) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function to get longitude and latitude coordinates for each pixel in the multi-looked.
    interferogram dataset. Coordinates are identical for each interferogram in the stack.
    :param ifg: pyrate.core.shared.Ifg Class object.
    :return: lon: Longitude for each pixel (decimal degrees)
    :return: lat: Latitude for each pixel (decimal degrees)
    """
    # assume all interferograms have same projection and will share the same transform
    transform = ifg.dataset.GetGeoTransform()
    # number of rows and columns in dataset
    nrows, ncols = ifg.shape
    lon = np.zeros((nrows, ncols))  # pre-allocate 2D numpy array
    lat = np.zeros((nrows, ncols))  # pre-allocate 2D numpy array
    for i in range(0, nrows):  # rows are y-direction
        for j in range(0, ncols):  # cols are x-direction
            lon[i, j], lat[i, j] = convert_pixel_value_to_geographic_coordinate(j, i, transform)

    return lon, lat


def test_get_lonlat_coords_vectorised(dem):
    lon, lat = get_lonlat_coords_slow(dem)
    lon_v, lat_v = get_lonlat_coords(dem)
    np.testing.assert_array_almost_equal(lon, lon_v.data)
    np.testing.assert_array_almost_equal(lat, lat_v.data)


@pytest.fixture(params=[(29, 50, -2.94634866714478), (94, 58, -2.94684600830078)])
def point_azimuth(request):
    return request.param


@pytest.fixture(params=[(29, 50, 1.02217936515808), (94, 58, 1.0111095905304)])
def point_incidence(request):
    return request.param


class TestPyRateAngleFiles:

    @classmethod
    def setup_class(cls):
        cls.params = Configuration(common.MEXICO_CROPA_CONF).__dict__
        # run prepifg
        prepifg.main(cls.params)
        # copy IFGs to temp folder
        correct._copy_mlooked(cls.params)

    @classmethod
    def teardown_class(cls):
        shutil.rmtree(cls.params[C.OUT_DIR], ignore_errors=True)

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

    def get_pyrate_angle(self, x0, y0, tif_file):
        """
        Get angle at particular pixel in the azimuth/incidence tif file
        """
        # get azimuth angle value of PyRate file azimuth_angle.tif
        temp = run(['gdallocationinfo', tif_file, str(x0), str(y0)], universal_newlines=True, stdout=PIPE)
        out = temp.stdout
        angle_value = float(out.split('Value:')[1].strip('\n'))

        return angle_value

    def test_pyrate_azimuth_matches_gamma_azimuth(self, point_azimuth):
        x0, y0, exp = point_azimuth
        # azimuth angle validation value extracted from GAMMA azimuth file for pixel location corresponding
        # to (29,50) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_phi.tif 613 3228
        # exp = -2.94634866714478  # GAMMMA azimuth angle at pixel location
        # azimuth angle validation value extracted from GAMMA azimuth file for pixel location corresponding
        # to (94,58) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_phi.tif 940 3271
        # exp = -2.94684600830078
        # GAMMA azimuth is defined towards satellite in an anti-clockwise angle system, with East being zero
        # PyRate azimuth angle is defined towards the satellite in a clockwise angle system with North being zero
        tif_file = join(self.params[C.GEOMETRY_DIR], 'azimuth_angle.tif')
        azimuth_angle_pyrate = self.get_pyrate_angle(x0, y0, tif_file)
        # convert PyRate azimuth into GAMMA azimuth
        res = -(azimuth_angle_pyrate - np.pi / 2)
        np.testing.assert_array_almost_equal(exp, res, decimal=2)  # max difference < 0.01 rad
        # screen output of difference if need be:
        # print(exp - res)

    def test_pyrate_incidence_matches_gamma_incidence(self, point_incidence):
        x0, y0, exp = point_incidence
        # incidence angle validation value extracted from GAMMA incidence file for pixel location corresponding
        # to (29,50) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_theta.tif 613 3228
        # exp = 1.02217936515808
        # incidence angle validation value extracted from GAMMA incidence file for pixel location corresponding
        # to (94,58) using gdallocationinfo 20180106_VV_8rlks_eqa_lv_theta.tif 940 3271
        # exp = 1.0111095905304
        # GAMMA angle is defined from the horizontal plane with the zenith direction being pi / 2 radians (i.e.90 deg)
        # PyRate angle is defined from the vertical axis with the zenith direction being zero
        tif_file = join(self.params[C.GEOMETRY_DIR], 'incidence_angle.tif')
        incidence_angle_pyrate = self.get_pyrate_angle(x0, y0, tif_file)
        # convert PyRate incidence into GAMMA incidence
        res = np.pi / 2 - incidence_angle_pyrate
        # convert PyRate incidence into GAMMA incidence
        np.testing.assert_array_almost_equal(exp, res, decimal=3)  # max difference < 0.001 rad
        # screen output of difference if need be:
        # print(exp - res)

    def test_azimuth_angle_calculation(self):
        """
         Calculate local azimuth angle using a spherical model and compare to result using Vincenty's equations
        """
        # get first IFG in stack to calculate lon/lat values
        multi_paths = self.params[C.INTERFEROGRAM_FILES]
        tmp_paths = [ifg_path.tmp_sampled_path for ifg_path in multi_paths]
        # keep only ifg files in path list (i.e. remove coherence and dem files)
        ifg_paths = [item for item in tmp_paths if 'ifg.tif' in item]
        # read and open the first IFG in list
        ifg0_path = ifg_paths[0]
        ifg0 = Ifg(ifg0_path)
        ifg0.open(readonly=True)
        lon, lat = get_lonlat_coords(ifg0)

        # read incidence and look angle files
        tif_file = join(self.params[C.GEOMETRY_DIR], 'incidence_angle.tif')
        geom = Geometry(tif_file)
        incidence_angle = geom.data
        tif_file = join(self.params[C.GEOMETRY_DIR], 'look_angle.tif')
        geom = Geometry(tif_file)
        look_angle = geom.data
        # get metadata
        a = float(ifg0.meta_data[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES])
        b = float(ifg0.meta_data[ifc.PYRATE_SEMI_MINOR_AXIS_METRES])
        heading = float(ifg0.meta_data[ifc.PYRATE_HEADING_DEGREES])
        azimuth = float(ifg0.meta_data[ifc.PYRATE_AZIMUTH_DEGREES])

        # convert all angles from deg to radians
        lon = np.radians(lon.data)
        lat = np.radians(lat.data)
        heading = np.radians(heading)
        azimuth = np.radians(azimuth)

        # calculate satellite positions
        sat_lat, sat_lon = get_sat_positions(lat, lon, look_angle, incidence_angle, heading, azimuth)

        # calculate azimuth angle from pixel to satellite using spherical trigonometry
        # see Eq. 86 on page 4-11 in EARTH-REFERENCED AIRCRAFT NAVIGATION AND SURVEILLANCE ANALYSIS
        # (https://ntlrepository.blob.core.windows.net/lib/59000/59300/59358/DOT-VNTSC-FAA-16-12.pdf)
        azimuth_angle_spherical = np.arctan2(np.cos(sat_lat) * np.sin(sat_lon - lon), \
                                             np.sin(sat_lat) * np.cos(lat) - np.cos(sat_lat) * np.sin(lat) * \
                                             np.cos(sat_lon - lon))
        # add 2 pi in case an angle is below zero
        for azi in np.nditer(azimuth_angle_spherical, op_flags=['readwrite']):
            if azi < 0:
                azi[...] = azi + 2 * np.pi

        azimuth_angle_elliposidal = vincinv(lat, lon, sat_lat, sat_lon, a, b)

        # the difference between Vincenty's azimuth calculation and the spherical approximation is ~0.001 radians
        azimuth_angle_diff = azimuth_angle_spherical - azimuth_angle_elliposidal
        # max difference < 0.01 rad
        np.testing.assert_array_almost_equal(azimuth_angle_spherical, azimuth_angle_elliposidal, decimal=2)
