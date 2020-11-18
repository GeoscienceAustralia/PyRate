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


def setup():
    params = Configuration(common.MEXICO_CONF).__dict__
    # run prepifg
    prepifg.main(params)
    # copy IFGs to temp folder
    correct._copy_mlooked(params)
    # read radar azimuth and range tif files
    rdc_az_file = join(params[cf.OUT_DIR], 'rdc_azimuth.tif')
    geom_az = Geometry(rdc_az_file)
    geom_az.open(readonly=True)
    az = geom_az.geometry_data
    rdc_rg_file = join(params[cf.OUT_DIR], 'rdc_range.tif')
    geom_rg = Geometry(rdc_rg_file)
    geom_rg.open(readonly=True)
    rg = geom_rg.geometry_data

    return az, rg, params


@pytest.fixture
def gamma_bperp(x0=50, y0=29):
    # calculate Bperp from GAMMA out files (interpolation required)
    # interpolation location (azimuth and range coordinate of crop A pixel (50,29)
    az, rg, params = setup()
    az0 = az[x0, y0]
    rg0 = rg[x0, y0]

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


@pytest.fixture
def pyrate_bperp(x0=50, y0=29):
    # calculate Bperp using PyRate functions

    az, rg, params = setup()
    multi_paths = params[cf.INTERFEROGRAM_FILES]
    tmp_paths = [ifg_path.tmp_sampled_path for ifg_path in multi_paths]
    # keep only ifg files in path list (i.e. remove coherence and dem files)
    ifg_paths = [item for item in tmp_paths if 'ifg.tif' in item]
    # read and open the first IFG in list
    ifg0_path = ifg_paths[0]
    ifg0 = Ifg(ifg0_path)
    ifg0.open(readonly=True)
    # geometry information needed to calculate Bperp for each pixel using first IFG in list
    lon, lat = geom.get_lonlat_coords(ifg0)

    # values for test pixel at 50, 29
    az0 = az[x0, y0]
    rg0 = rg[x0, y0]
    lon0 = lon[x0, y0]
    lat0 = lat[x0, y0]

    # calculate Bperp
    nifgs = len(ifg_paths)
    bperp = np.empty((nifgs)) * np.nan
    # calculate per-pixel perpendicular baseline for each IFG
    for ifg_num, ifg_path in enumerate(
            ifg_paths):  # loop could be avoided by approximating the look angle for the first Ifg
        ifg = Ifg(ifg_path)
        ifg.open(readonly=True)
        # calculate look angle for interferograms (using the Near Range of the primary SLC)
        look_angle, range_dist = geom.write_local_geometry_files(ifg, None, rg0, lon0,
                                                                 lat0, params)
        bperp[ifg_num] = geom.calc_local_baseline(ifg, az0, look_angle)

    return bperp


def test_pyrate_bperp_matches_gamma_bperp(gamma_bperp, pyrate_bperp):
    # compare the GAMMA and PyRate Bperp estimates
    np.testing.assert_array_almost_equal(gamma_bperp / 1e6, pyrate_bperp / 1e6)  # max difference < 10mm
    # screen output if need be:
    diff = gamma_bperp - pyrate_bperp
    print(1000 * np.mean(diff))  # mean difference in mm
    print(1000 * np.max(diff))  # max difference in mm
    # Mean difference: 0.19 mm, maximum difference: 1.29 mm

