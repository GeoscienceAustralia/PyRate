import numpy as np
from scipy.interpolate import griddata
import math
import pytest
from pyrate.constants import PYRATEPATH

test_data = PYRATEPATH.joinpath('tests', 'test_data')
geotiffs = test_data.joinpath('geotiffs')
geometry = test_data.joinpath('geometry')


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
    bperp_files = sorted(list(geometry.glob('*_bperp.par')))
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
pyrate_bperp = np.array([33.48592183, 3.44669685, -75.37369399, -26.88597679, -33.25298942, -108.84360354, 3.74075472, \
                         -3.19700977, -14.58390611, 10.06920291, -51.12649599, -5.74544068, -17.36872483, -30.4772929, \
                         7.12691256, -37.68943916, -73.14248882, -11.45674522, -24.64851804, 12.69928323, -32.16248418, \
                         -20.86746046, 61.514626, 48.30928659, -13.17640207, 24.28126177, -36.84111057, -20.5870326, \
                         77.8291117, -8.66115426])


# compare the GAMMA and PyRate Bperp estimates
def test_pyrate_bperp_matches_gamma_bperp(gamma_bperp):
    np.testing.assert_array_almost_equal(pyrate_bperp/1e6, gamma_bperp/1e6)  # max difference < 10mm


# diff = bperp_int - pyrate_bperp
# print(1000 * np.mean(diff))  # mean difference in mm
# print(1000 * np.max(diff))  # max difference in mm
# Mean difference: 0.19 mm, maximum difference: 1.29 mm
