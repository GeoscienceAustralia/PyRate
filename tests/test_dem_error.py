import numpy as np
import glob
import os
from scipy.interpolate import griddata
import math
import pytest
import shutil
from tests import common
from pyrate.configuration import Configuration
from pyrate import prepifg, correct
import pyrate.core.config as cf
from pyrate.core import shared
import pyrate.core.geometry as geom
from pyrate.core.shared import Ifg, Geometry


geometry_path = common.MEXICO_TEST_DIR_GEOMETRY


@pytest.fixture
def gamma_bperp():
    # interpolation location (azimuth and range coordinate of crop A pixel (50,29)
    az0 = 2636.9912
    rg0 = 103.4759

    # round azimuth and range coordinates to closest step (500 for az, 200 for rg)
    azstep = 500
    rgstep = 200
    az1 = azstep * math.floor(az0 / azstep)
    rg1 = rgstep * math.floor(rg0 / rgstep)
    az2 = azstep * math.ceil(az0 / azstep)
    rg2 = rgstep * math.ceil(rg0 / rgstep)

    # four coordinates for bi-linear interpolation
    teststr1 = str(az1).rjust(6, ' ') + str(rg1).rjust(7, ' ')
    teststr2 = str(az1).rjust(6, ' ') + str(rg2).rjust(7, ' ')
    teststr3 = str(az2).rjust(6, ' ') + str(rg1).rjust(7, ' ')
    teststr4 = str(az2).rjust(6, ' ') + str(rg2).rjust(7, ' ')

    # loop through all corresponding bperp.par files in base.par list
    bperp_files = sorted(list(glob.glob(os.path.join(geometry_path, '*_bperp.par'))))
    bperp_int = np.empty(shape=(len(bperp_files)))
    for i, bperp_file in enumerate(bperp_files):
        # read Mexico city bperp file
        with open(bperp_file, 'r') as f:
            for line in f.readlines():
                if teststr1 in line:
                    bperp1 = line.split()[7]
                if teststr2 in line:
                    bperp2 = line.split()[7]
                if teststr3 in line:
                    bperp3 = line.split()[7]
                if teststr4 in line:
                    bperp4 = line.split()[7]

        # setup numpy array for bi-linear interpolation
        n = np.array([(az1, rg1, bperp1),
                      (az1, rg2, bperp2),
                      (az2, rg1, bperp3),
                      (az2, rg2, bperp4)])
        # interpolate using scipy function "griddata"
        bperp_int[i] = griddata(n[:, 0:2], n[:, 2], [(az0, rg0)], method='linear')

    return bperp_int


# TODO: calculate this with pyrate dem_error
pyrate_bperp = np.array([33.48592183, 3.44669685, -75.37369399, -26.88597679, -33.25298942, -108.84360354, 3.74075472,
                         -3.19700977, -14.58390611, 10.06920291, -51.12649599, -5.74544068, -17.36872483, -30.4772929,
                         7.12691256, -37.68943916, -73.14248882, -11.45674522, -24.64851804, 12.69928323, -32.16248418,
                         -20.86746046, 61.514626, 48.30928659, -13.17640207, 24.28126177, -36.84111057, -20.5870326,
                         77.8291117, -8.66115426])


# compare the GAMMA and PyRate Bperp estimates
def test_pyrate_bperp_matches_gamma_bperp(gamma_bperp):
    np.testing.assert_array_almost_equal(pyrate_bperp / 1e6, gamma_bperp / 1e6)  # max difference < 10mm

# diff = bperp_int - pyrate_bperp
# print(1000 * np.mean(diff))  # mean difference in mm
# print(1000 * np.max(diff))  # max difference in mm
# Mean difference: 0.19 mm, maximum difference: 1.29 mm

class TestDEMerror:

    @classmethod
    @pytest.fixture(autouse=True)
    def setup_method(cls):
        cls.conf = common.MEXICO_CONF
        params = Configuration(cls.conf).__dict__
        prepifg.main(params)
        cls.params = Configuration(cls.conf).__dict__
        correct._copy_mlooked(cls.params)
        correct._update_params_with_tiles(cls.params)
        correct._create_ifg_dict(cls.params)
        multi_paths = cls.params[cf.INTERFEROGRAM_FILES]
        cls.ifg_paths = [p.tmp_sampled_path for p in multi_paths]
        cls.ifgs = [shared.Ifg(i) for i in cls.ifg_paths]
        for i in cls.ifgs:
            i.open()
        shared.save_numpy_phase(cls.ifg_paths, cls.params)
        correct.mst_calc_wrapper(cls.params)

    @classmethod
    def teardown_method(cls):
        pass

    def calc_pyrate_bperp(self, row, col):

        print(self)
        # read radar azimuth and range tif files
        rdc_az_file = join(self.params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        az = geom_az.geometry_data
        rdc_rg_file = join(self.params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        rg = geom_rg.geometry_data

        ifg0 = self.ifgs[0]
        lon, lat = geom.get_lonlat_coords(ifg0)

        # values for test pixel at 50, 29
        az_pix = az[50, 29]
        rg_pix = rg[50, 29]
        lon_pix = lon[50, 29]
        lat_pix = lat[50, 29]

        ifg_paths = self.ifg_paths
        nifgs = len(ifg_paths)
        bperp_pyrate = np.empty((nifgs, 1, 1)) * np.nan
        # calculate per-pixel perpendicular baseline for each IFG
        for ifg_num, ifg_path in enumerate(
                ifg_paths):  # loop could be avoided by approximating the look angle for the first Ifg
            ifg = Ifg(ifg_path)
            ifg.open(readonly=True)
            # calculate look angle for interferograms (using the Near Range of the primary SLC)
            look_angle, range_dist = geom.write_local_geometry_files(ifg, None, rg_pix, lon_pix,
                                                                         lat_pix, self.params)
            bperp_pyrate[ifg_num, :, :] = geom.calc_local_baseline(ifg, az_pix, look_angle)

