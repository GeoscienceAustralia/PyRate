import numpy as np
import glob
import os
from os.path import join
from scipy.interpolate import griddata
import math
import pytest
from tests import common
from pyrate.configuration import Configuration
from pyrate import prepifg, correct
import pyrate.core.config as cf
import pyrate.core.geometry as geom
from pyrate.core.shared import Ifg, Geometry


geometry_path = common.MEXICO_TEST_DIR_GEOMETRY


def gamma_bperp(self, x0, y0):
    """
    Calculate Bperp for specified pixel from GAMMA out files (interpolation required)
    x0, y0 is the interpolation location in azimuth and range
    """
    az0 = self.az[x0, y0]
    rg0 = self.rg[x0, y0]

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


def pyrate_bperp(self):
    """
    Calculate Bperp image for each ifg using PyRate functions
    """
    multi_paths = self.params[cf.INTERFEROGRAM_FILES]
    tmp_paths = [ifg_path.tmp_sampled_path for ifg_path in multi_paths]
    # keep only ifg files in path list (i.e. remove coherence and dem files)
    ifg_paths = [item for item in tmp_paths if 'ifg.tif' in item]
    # read and open the first IFG in list
    ifg0_path = ifg_paths[0]
    ifg0 = Ifg(ifg0_path)
    ifg0.open(readonly=True)
    # size of ifg dataset
    nrows, ncols = ifg0.shape
    nifgs = len(ifg_paths)
    bperp = np.empty(shape=(nrows, ncols, nifgs)) * np.nan

    # calculate per-pixel perpendicular baseline for each IFG
    for ifg_num, ifg_path in enumerate(
            ifg_paths):  # loop could be avoided by approximating the look angle for the first Ifg
        ifg = Ifg(ifg_path)
        ifg.open(readonly=True)
        # calculate look angle for interferograms (using the Near Range of the primary SLC)
        look_angle, _, _, _ = geom.calc_pixel_geometry(ifg, self.rg, self.params)
        bperp[:, :, ifg_num] = geom.calc_local_baseline(ifg, self.az, look_angle)

    return bperp

class TestPyRateGammaBperp:

    @classmethod
    def setup_method(cls):
        cls.params = Configuration(common.MEXICO_CONF).__dict__
        # run prepifg
        prepifg.main(cls.params)
        # copy IFGs to temp folder
        correct._copy_mlooked(cls.params)
        # read radar azimuth and range tif files
        rdc_az_file = join(cls.params[cf.OUT_DIR], 'rdc_azimuth.tif')
        geom_az = Geometry(rdc_az_file)
        geom_az.open(readonly=True)
        cls.az = geom_az.geometry_data
        rdc_rg_file = join(cls.params[cf.OUT_DIR], 'rdc_range.tif')
        geom_rg = Geometry(rdc_rg_file)
        geom_rg.open(readonly=True)
        cls.rg = geom_rg.geometry_data
        # calc bperp using pyrate funcs
        cls.pbperp = pyrate_bperp(cls)

    @classmethod
    def teardown_method(cls):
        pass

    def test_pyrate_bperp_matches_gamma_bperp_x50_y29(self):
        x0 = 50; y0 = 29
        res = self.pbperp[x0, y0, :]
        exp = gamma_bperp(self, x0, y0)
        np.testing.assert_array_almost_equal(exp, res, 3)  # max difference < 10mm
        # screen output if need be:
        #diff = exp - res
        #print(1000 * np.mean(diff))  # mean difference in mm
        #print(1000 * np.max(diff))  # max difference in mm
        # Mean difference: 0.19 mm, maximum difference: 1.29 mm

    def test_pyrate_bperp_matches_gamma_bperp_x59_y99(self):
        x0 = 59; y0 = 99
        res = self.pbperp[x0, y0, :]
        exp = gamma_bperp(self, x0, y0)
        np.testing.assert_array_almost_equal(exp, res, 3)

    def test_pyrate_bperp_matches_gamma_bperp_x3_y2(self):
        x0 = 3; y0 = 2
        res = self.pbperp[x0, y0, :]
        exp = gamma_bperp(self, x0, y0)
        np.testing.assert_array_almost_equal(exp, res, 2) # Max absolute difference: 0.00239929

