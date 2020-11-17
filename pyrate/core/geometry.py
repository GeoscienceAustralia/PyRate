#   This Python module is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
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
This Python module implements the calculation and output of the per-pixel vector of the radar viewing geometry
 (i.e. local angles, incidence angles and azimuth angles) as well as the calculation of per-pixel baseline values
 used for correcting interferograms for residual topographic effects (a.k.a. as DEM errors).
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import numpy as np
import os
from math import sqrt, sin, cos, tan, asin, atan, atan2, isnan, pi
from pyrate.core import ifgconstants as ifc, config as cf
from pyrate.core.refpixel import convert_pixel_value_to_geographic_coordinate
from pyrate.core.gamma import read_lookup_table


def get_lonlat_coords_slow(ifg):
    """
    Function to get longitude and latitude coordinates for each pixel in the multi-looked.
    interferogram dataset. Coordinates are identical for each interferogram in the stack.
    """
    # assume all interferograms have same projection and will share the same transform
    transform = ifg.dataset.GetGeoTransform()
    # number of rows and columns in dataset
    nrows, ncols = ifg.shape
    lon = np.zeros((nrows, ncols))  # pre-allocate 2D numpy array
    lat = np.zeros((nrows, ncols))  # pre-allocate 2D numpy array
    for i in range(0, nrows): # rows are y-direction
        for j in range(0, ncols): # cols are x-direction
            lon[i, j], lat[i, j] = convert_pixel_value_to_geographic_coordinate(j, i, transform)

    return lon, lat


def get_lonlat_coords(ifg):
    """
    Function to get longitude and latitude coordinates for each pixel in the multi-looked.
    interferogram dataset. Coordinates are identical for each interferogram in the stack.
    """
    # assume all interferograms have same projection and will share the same transform
    transform = ifg.dataset.GetGeoTransform()
    # number of rows and columns in dataset
    nrows, ncols = ifg.shape
    yOrigin = transform[3]
    pixelHeight = -transform[5]
    xOrigin = transform[0]
    pixelWidth = transform[1]

    lons = np.arange(0, ncols) * pixelWidth + xOrigin
    lats = yOrigin - np.arange(0, nrows) * pixelHeight
    lon, lat = np.meshgrid(lons, lats)

    return lon, lat


def calc_radar_coords(ifg, params, xmin, xmax, ymin, ymax):
    """
    Function to calculate radar coordinates for each pixel in the multi-looked
    interferogram dataset. Radar coordinates are identical for each interferogram
    in the stack.
    """
    # lookup table file:
    lookup_table = params[cf.LT_FILE]
    # PyRate IFG multi-looking factors
    ifglksx = params[cf.IFG_LKSX]
    ifglksy = params[cf.IFG_LKSY]
     # transform float lookup table file to np array, min/max pixel coordinates are required for cropping
    lt_az, lt_rg = read_lookup_table(ifg, lookup_table, ifglksx, ifglksy, xmin, xmax, ymin, ymax)
    # replace 0.0 with NaN
    lt_az[lt_az == 0.0] = np.nan
    lt_rg[lt_rg == 0.0] = np.nan

    return lt_az, lt_rg


def calc_pixel_geometry(ifg, rg, params):
    """
    Function to calculate local look angle, incidence angle and geodetic azimuth for each pixel.
    """
    # read relevant metadata from first IFG
    a = float(ifg.meta_data[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES])
    b = float(ifg.meta_data[ifc.PYRATE_SEMI_MINOR_AXIS_METRES])
    se = float(ifg.meta_data[ifc.PYRATE_SAR_EARTH_METRES])
    # near range of primary image is stored in the interferogram metadata
    near_range = float(ifg.meta_data[ifc.PYRATE_NEAR_RANGE_METRES])
    rps = float(ifg.meta_data[ifc.PYRATE_RANGE_PIX_METRES])
    heading = float(ifg.meta_data[ifc.PYRATE_HEADING_DEGREES])
    azimuth = float(ifg.meta_data[ifc.PYRATE_AZIMUTH_DEGREES])

    # calculate per-pixel lon/lat
    lon, lat = get_lonlat_coords(ifg)

    # convert angles to radians
    lon = np.radians(lon)
    lat = np.radians(lat)
    heading = np.radians(heading)
    azimuth = np.radians(azimuth)

    # Earth radius at given latitude
    re = np.sqrt(np.divide(np.square(a**2 * np.cos(lat)) + np.square(b**2 * np.sin(lat)), \
                           np.square(a * np.cos(lat)) + np.square(b * np.sin(lat))))

    # range measurement at pixel ij
    range_dist = near_range + rps * rg

    # look angle at pixel ij -> law of cosines in "satellite - Earth centre - ground pixel" triangle
    # see e.g. Section 2 in https://www.cs.uaf.edu/~olawlor/ref/asf/sar_equations_2006_08_17.pdf
    look_angle = np.arccos(np.divide(se**2 + np.square(range_dist) - np.square(re), 2 * se * range_dist))

    # incidence angle at pixel ij -> law of cosines in "satellite - Earth centre - ground pixel" triangle
    # see e.g. Section 2 in https://www.cs.uaf.edu/~olawlor/ref/asf/sar_equations_2006_08_17.pdf
    incidence_angle = np.pi - np.arccos(np.divide(np.square(range_dist) + np.square(re) - se**2, \
                                                  2 * np.multiply(range_dist, re)))

    # todo (once new test data is ready): move next line into test for validation with GAMMA output
    #incidence_angle_gamma = np.pi / 2 - incidence_angle
    # maximum differences to the GAMMA-derived local incidence angles for Sentinel-1 test data are within +/-0.1 deg
    # to improve the accuracy further one would need to consider the local slope observed from the DEM

    # local azimuth angle at pixel ij using constant satellite heading angle and spherical approximations
    epsilon = np.pi - look_angle - (np.pi - incidence_angle) # angle at the Earth's center between se and re
    # azimuth of satellite look vector (satellite heading + azimuth of look direction (+90 deg for right-looking SAR)
    sat_azi = heading + azimuth
    # the following equations are adapted from Section 4.4 (page 4-16) in EARTH-REFERENCED AIRCRAFT NAVIGATION AND
    # SURVEILLANCE ANALYSIS (https://ntlrepository.blob.core.windows.net/lib/59000/59300/59358/DOT-VNTSC-FAA-16-12.pdf)
    sat_lon = np.divide(np.arcsin(-(np.multiply(np.sin(epsilon), np.sin(sat_azi)))), np.cos(lat)) + lon # Eq. 103
    temp = np.multiply(np.divide(np.cos(0.5 * (sat_azi + sat_lon - lon)), np.cos(0.5 * (sat_azi - sat_lon + lon))), \
                       np.tan(0.5 * (np.pi / 2 + lat - epsilon))) # Eq. 104
    sat_lat = -np.pi / 2 + 2 * np.arctan(temp)
    # the above code could be improved by calculating satellite positions for each azimuth row
    # from orbital state vectors given in .par file, using the following workflow:
    # 1. read orbital state vectors and start/stop times from mli.par
    # 2. for each pixel get the corresponding radar row (from matrix az)
    # 3. get the corresponding radar time for that row (using linear interpolation)
    # 4. calculate the satellite XYZ position for that time by interpolating the time and velocity state vectors

    # calc azimuth angle using Vincenty's equations
    if np.isscalar(lat): # function works also for scalar input instead of numpy array
        azimuth_angle = vincinv(lat, lon, sat_lat, sat_lon, a, b)
    else:
        azimuth_angle = np.empty(lat.shape) * np.nan # pre-allocate 2D numpy array
        for ix, iy in np.ndindex(lat.shape):
            if not isnan(sat_lat[ix, iy]):
                az12 = vincinv(lat[ix, iy], lon[ix, iy], sat_lat[ix, iy], sat_lon[ix, iy], a, b)
                azimuth_angle[ix, iy] = az12
        np.reshape(azimuth_angle, lat.shape)

    # todo: move this old code for azimuth angle calculation using a spherical Earth model to tests
    #azimuth_angle2 = np.arccos(np.divide(np.multiply(np.sin(sat_lat), np.cos(lat)) - \
    #                                np.multiply(np.multiply(np.cos(sat_lat), np.sin(lat)), np.cos(sat_lon - lon)), \
    #                                np.sin(epsilon)))
    #azimuth_angle_diff = azimuth_angle2 - azimuth_angle
    #print(np.nanmin(azimuth_angle_diff), np.nanmax(azimuth_angle_diff))
    # the difference between Vincenty's azimuth calculation and the spherical approximation is ~0.001 radians

    # todo (once new test data is ready): move next line into test for validation with GAMMA output
    #azimuth_angle_gamma = -(azimuth_angle - np.pi / 2) # local azimuth towards satellite as output by GAMMA
    # maximum differences to the GAMMA-derived local azimuth angles for Sentinel-1 test data are within +/-0.5 deg
    # this could be improved by using orbital state vectors to calculate that satellite positions (see above comment)

    return look_angle, incidence_angle, azimuth_angle, range_dist


def calc_local_baseline(ifg, az, look_angle):
    """
    Function to calculate perpendicular baseline values for each pixel.
    """
    # open ifg object
    if not ifg.is_open:
        ifg.open()

    # read relevant metadata from IFG
    prf = float(ifg.meta_data[ifc.PYRATE_PRF_HERTZ])
    az_looks = int(ifg.meta_data[ifc.PYRATE_AZIMUTH_LOOKS])
    az_n = int(ifg.meta_data[ifc.PYRATE_AZIMUTH_N])
    base_C = float(ifg.meta_data[ifc.PYRATE_BASELINE_C])
    base_N = float(ifg.meta_data[ifc.PYRATE_BASELINE_N])
    baserate_C = float(ifg.meta_data[ifc.PYRATE_BASELINE_RATE_C])
    baserate_N = float(ifg.meta_data[ifc.PYRATE_BASELINE_RATE_N])

    # calculate per pixel baseline vectors across track (C) and normal to the track (N)
    mean_az = az_n / 2 - 0.5 # mean azimuth line
    prf = prf / az_looks # Pulse Repetition Frequency needs to be adjusted according to GAMMA azimuth looks
    base_C_local = base_C + baserate_C * (az - mean_az) / prf
    base_N_local = base_N + baserate_N * (az - mean_az) / prf

    # calculate the per-pixel perpendicular baseline (see Eq. 3.5 in Baehr, 2012 available here:
    # http://www.dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/c-719.pdf)
    bperp = np.multiply(base_C_local, np.cos(look_angle)) - np.multiply(base_N_local, np.sin(look_angle))

    return bperp


def vincinv(lat1, lon1, lat2, lon2, semimaj, semimin):
    """
    Vincenty's Inverse Formula, adapted from GeodePy function vincinv
    (see https://github.com/GeoscienceAustralia/GeodePy/blob/master/geodepy/geodesy.py)
    :param lat1: Latitude of Point 1 (radians)
    :param lon1: Longitude of Point 1 (radians)
    :param lat2: Latitude of Point 2 (radians)
    :param lon2: Longitude of Point 2 (radians)
    :param semimaj: semi-major axis of ellipsoid
    :param semimin: semi-minor axis of ellipsoid
    :return: azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees)
    """
    # Exit if the two input points are the same
    if lat1 == lat2 and lon1 == lon2:
        return 0, 0, 0
    # calculate flattening
    f = (semimaj-semimin)/semimaj
    # Equation numbering is from the GDA2020 Tech Manual v1.0
    # Eq. 71
    u1 = atan((1 - f) * tan(lat1))
    # Eq. 72
    u2 = atan((1 - f) * tan(lat2))
    # Eq. 73; initial approximation
    lon = lon2 - lon1
    omega = lon
    # Iterate until the change in lambda, lambda_sigma, is insignificant
    # (< 1e-12) or after 1000 iterations have been completed
    for i in range(1000):
        # Eq. 74
        sin_sigma = sqrt((cos(u2)*sin(lon))**2 + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lon))**2)
        # Eq. 75
        cos_sigma = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lon)
        # Eq. 76
        sigma = atan2(sin_sigma, cos_sigma)
        # Eq. 77
        alpha = asin((cos(u1)*cos(u2)*sin(lon)) / sin_sigma)
        # Eq. 78
        cos_two_sigma_m = cos(sigma) - (2*sin(u1)*sin(u2) / cos(alpha)**2)
        # Eq. 79
        c = (f / 16) * cos(alpha)**2 * (4 + f * (4 - 3*cos(alpha)**2))
        # Eq. 80
        new_lon = omega + (1 - c) * f * sin(alpha) * (
                sigma + c*sin(sigma) * (cos_two_sigma_m + c * cos(sigma) * (-1 + 2*(cos_two_sigma_m**2)))
        )
        delta_lon = new_lon - lon
        lon = new_lon
        if abs(delta_lon) < 1e-12:
            break
    # Calculate the azimuth from point 1 to point 2
    azimuth1to2 = atan2((cos(u2)*sin(lon)), (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lon)))
    if azimuth1to2 < 0:
        azimuth1to2 = azimuth1to2 + 2*pi

    return round(azimuth1to2, 9)
