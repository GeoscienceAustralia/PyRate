#   This Python module is part of the PyRate software package.
#
#   Copyright 2022 Geoscience Australia
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
used for correcting interferograms for residual topographic effects (a.k.a. DEM errors).
"""
# pylint: disable=invalid-name, too-many-locals, too-many-arguments
import numpy as np
from typing import Tuple, Union

import pyrate.constants as C
from pyrate.core import ifgconstants as ifc
from pyrate.core.gamma import read_lookup_table
from pyrate.core.shared import Ifg, IfgPart, MemGeometry


def get_lonlat_coords(ifg: Ifg) -> Tuple[MemGeometry, MemGeometry]:
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
    yOrigin = transform[3]
    pixelHeight = -transform[5]
    xOrigin = transform[0]
    pixelWidth = transform[1]

    lons = np.arange(0, ncols) * pixelWidth + xOrigin
    lats = yOrigin - np.arange(0, nrows) * pixelHeight
    lon, lat = np.meshgrid(lons, lats)

    return MemGeometry(lon), MemGeometry(lat)


def calc_radar_coords(ifg: Ifg, params: dict, xmin: int, xmax: int,
                      ymin: int, ymax: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function to calculate radar coordinates for each pixel in the multi-looked
    interferogram dataset. Radar coordinates are identical for each interferogram
    in the stack. Uses the Gamma lookup table defined in the configuration file.
    :param ifg: pyrate.core.shared.Ifg Class object.
    :param params: Dictionary of PyRate configuration parameters.
    :param xmin: Minimum longitude of cropped image (decimal degrees).
    :param xmax: Maximum longitude of cropped image (decimal degrees)
    :param ymin: Minimum latitude of cropped image (decimal degrees).
    :param ymax: Maximum latitude of cropped image (decimal degrees)
    :return: lt_az: Radar geometry azimuth coordinate for each pixel.
    :return: lt_rg: Radar geometry range coordinate for each pixel.
    """
    # lookup table file:
    lookup_table = params[C.LT_FILE]

    if lookup_table is None:
        msg = f"No lookup table file supplied: Geometry cannot be computed"
        raise FileNotFoundError(msg)

    # PyRate IFG multi-looking factors
    ifglksx = params[C.IFG_LKSX]
    ifglksy = params[C.IFG_LKSY]
    # transform float lookup table file to np array, min/max pixel coordinates are required for cropping
    lt_az, lt_rg = read_lookup_table(ifg, lookup_table, ifglksx, ifglksy, xmin, xmax, ymin, ymax)
    # replace 0.0 with NaN
    lt_az[lt_az == 0.0] = np.nan
    lt_rg[lt_rg == 0.0] = np.nan

    return lt_az, lt_rg


def get_sat_positions(lat: np.ndarray, lon: np.ndarray, look_angle: np.ndarray, inc_angle: np.ndarray,
                      heading: np.float64, look_dir: np.float64) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function to calculate the lon/lat position of the satellite for each pixel.
    :param lat: Ground latitude for each pixel (decimal degrees).
    :param lon: Ground point longitude for each pixel (decimal degrees).
    :param look_angle: Look angle (between nadir and look vector) for each pixel (radians).
    :param inc_angle: Local incidence angle (between vertical and look vector) for each pixel (radians).
    :param heading: Satellite flight heading (radians).
    :param look_dir: Look direction w.r.t. satellite heading; +ve = right looking (radians).
    :return: sat_lat: Satellite position latitude for each pixel (decimal degrees).
    :return: sat_lon: Satellite position longitude for each pixel (decimal degrees).
    """
    # note that the accuracy of satellite lat/lon positions code could be improved by calculating satellite positions
    # for each azimuth row from orbital state vectors given in .par file, using the following workflow:
    # 1. read orbital state vectors and start/stop times from mli.par
    # 2. for each pixel get the corresponding radar row (from matrix az)
    # 3. get the corresponding radar time for that row (using linear interpolation)
    # 4. calculate the satellite XYZ position for that time by interpolating the time and velocity state vectors

    # angle at the Earth's center between se and re
    epsilon = np.pi - look_angle - (np.pi - inc_angle)
    # azimuth of satellite look vector (satellite heading + look direction (+90 deg for right-looking SAR)
    sat_azi = heading + look_dir
    # the following equations are adapted from Section 4.4 (page 4-16) in EARTH-REFERENCED AIRCRAFT NAVIGATION AND
    # SURVEILLANCE ANALYSIS (https://ntlrepository.blob.core.windows.net/lib/59000/59300/59358/DOT-VNTSC-FAA-16-12.pdf)
    sat_lon = np.divide(np.arcsin(-(np.multiply(np.sin(epsilon), np.sin(sat_azi)))), np.cos(lat)) + lon  # Eq. 103
    temp = np.multiply(np.divide(np.cos(0.5 * (sat_azi + sat_lon - lon)), np.cos(0.5 * (sat_azi - sat_lon + lon))), \
                       np.tan(0.5 * (np.pi / 2 + lat - epsilon)))  # Eq. 104
    sat_lat = -np.pi / 2 + 2 * np.arctan(temp)

    return np.real(sat_lat), np.real(sat_lon)


def calc_pixel_geometry(ifg: Union[Ifg, IfgPart], rg: np.ndarray, lon: np.ndarray, lat: np.ndarray,
                        dem_height: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Function to calculate angular satellite to ground geometries and distance for each pixel.
    :param ifg: pyrate.core.shared.Ifg Class object.
    :param rg: Range image coordinate for each pixel.
    :param lon: Longitude for each pixel (decimal degrees).
    :param lat: Latitude for each pixel (decimal degrees).
    :param dem_height: Height from DEM for each pixel (metres).
    :return: look_angle: Look angle (between nadir and look vector) for each pixel (radians).
    :return: incidence_angle: Local incidence angle (between vertical and look vector) for each pixel (radians).
    :return: azimuth_angle: Geodetic azimuth for each pixel (radians).
    :return: range_dist: Distance from satellite to ground for each pixel (metres).
    """
    # read relevant metadata from first IFG
    a = float(ifg.meta_data[ifc.PYRATE_SEMI_MAJOR_AXIS_METRES])
    b = float(ifg.meta_data[ifc.PYRATE_SEMI_MINOR_AXIS_METRES])
    se = float(ifg.meta_data[ifc.PYRATE_SAR_EARTH_METRES])
    # near range of primary image is stored in the interferogram metadata
    near_range = float(ifg.meta_data[ifc.PYRATE_NEAR_RANGE_METRES])
    rps = float(ifg.meta_data[ifc.PYRATE_RANGE_PIX_METRES])
    heading = float(ifg.meta_data[ifc.PYRATE_HEADING_DEGREES])
    # direction of look vector w.r.t. satellite heading. 
    # Gamma convention: +ve = right; -ve = left.
    look_dir = float(ifg.meta_data[ifc.PYRATE_AZIMUTH_DEGREES])

    # convert to radians
    lon = np.radians(lon)
    lat = np.radians(lat)
    heading = np.radians(heading)
    look_dir = np.radians(look_dir)

    # Earth radius at given latitude
    re = np.sqrt(np.divide(np.square(a ** 2 * np.cos(lat)) + np.square(b ** 2 * np.sin(lat)),
                           np.square(a * np.cos(lat)) + np.square(b * np.sin(lat))))

    # range measurement at pixel ij
    range_dist = near_range + rps * rg

    # look angle at pixel ij -> law of cosines in "satellite - Earth centre - ground pixel" triangle
    # see e.g. Section 2 in https://www.cs.uaf.edu/~olawlor/ref/asf/sar_equations_2006_08_17.pdf
    look_angle = np.arccos(np.divide(se ** 2 + np.square(range_dist) - np.square(re), 2 * se * range_dist))

    # add per-pixel height to the earth radius to obtain a more accurate ground pixel position for
    # incidence angle calculation
    re = re + dem_height

    # incidence angle at pixel ij -> law of cosines in "satellite - Earth centre - ground pixel" triangle
    # see e.g. Section 2 in https://www.cs.uaf.edu/~olawlor/ref/asf/sar_equations_2006_08_17.pdf
    incidence_angle = np.pi - np.arccos(np.divide(np.square(range_dist) + np.square(re) - se ** 2,
                                                  2 * np.multiply(range_dist, re)))

    # calculate satellite positions for each pixel
    sat_lat, sat_lon = get_sat_positions(lat, lon, look_angle, incidence_angle, heading, look_dir)

    # # calc azimuth angle using Vincenty's equations
    azimuth_angle = vincinv(lat, lon, sat_lat, sat_lon, a, b)

    return look_angle, incidence_angle, azimuth_angle, range_dist


def calc_local_baseline(ifg: Ifg, az: np.ndarray, look_angle: np.ndarray) -> np.ndarray:
    """
    Function to calculate perpendicular baseline values for each pixel.
    :param ifg: pyrate.core.shared.Ifg Class object
    :param az: Azimuth image coordinate for each pixel
    :param look_angle: Look angle for each pixel (radians)
    :return: bperp: Perpendicular baseline for each pixel (metres)
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
    mean_az = az_n / 2 - 0.5  # mean azimuth line
    prf = prf / az_looks  # Pulse Repetition Frequency needs to be adjusted according to GAMMA azimuth looks
    base_C_local = base_C + baserate_C * (az - mean_az) / prf
    base_N_local = base_N + baserate_N * (az - mean_az) / prf

    # calculate the per-pixel perpendicular baseline (see Eq. 3.5 in Baehr, 2012 available here:
    # http://www.dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/c-719.pdf)
    bperp = np.multiply(base_C_local, np.cos(look_angle)) - np.multiply(base_N_local, np.sin(look_angle))

    return bperp


def vincinv(lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray,
            semimaj: float, semimin: float) -> np.ndarray:
    """
    Vincenty's Inverse Formula, adapted from GeodePy function vincinv
    (see https://github.com/GeoscienceAustralia/GeodePy/blob/master/geodepy/geodesy.py)
    - only relevant parts of the Geodepy to retrieve the azimuth angle have been used
    - vectorised the function for use with numpy arrays
    :param lat1: Latitude of Point 1 (radians)
    :param lon1: Longitude of Point 1 (radians)
    :param lat2: Latitude of Point 2 (radians)
    :param lon2: Longitude of Point 2 (radians)
    :param semimaj: semi-major axis of ellipsoid
    :param semimin: semi-minor axis of ellipsoid
    :return: azimuth1to2: Azimuth from Point 1 to 2 (Decimal Degrees)
    """
    # Exit if any of the two sets of input points are the same
    if np.any(lat1 == lat2) and np.any(lon1 == lon2):
        return 0
    # calculate flattening
    f = (semimaj - semimin) / semimaj
    # Equation numbering is from the GDA2020 Tech Manual v1.0
    # Eq. 71
    u1 = np.arctan((1 - f) * np.tan(lat1))
    # Eq. 72
    u2 = np.arctan((1 - f) * np.tan(lat2))
    # Eq. 73; initial approximation
    lon = lon2 - lon1
    omega = lon
    # Iterate until the change in lambda, lambda_sigma, is insignificant
    # (< 1e-12) or after 1000 iterations have been completed
    for i in range(1000):
        # Eq. 74
        sin_sigma = np.sqrt(
            (np.cos(u2) * np.sin(lon)) ** 2 + (np.cos(u1) * np.sin(u2) - np.sin(u1) * np.cos(u2) * np.cos(lon)) ** 2)
        # Eq. 75
        cos_sigma = np.sin(u1) * np.sin(u2) + np.cos(u1) * np.cos(u2) * np.cos(lon)
        # Eq. 76
        sigma = np.arctan2(sin_sigma, cos_sigma)
        # Eq. 77
        alpha = np.arcsin((np.cos(u1) * np.cos(u2) * np.sin(lon)) / sin_sigma)
        # Eq. 78
        cos_two_sigma_m = np.cos(sigma) - (2 * np.sin(u1) * np.sin(u2) / np.cos(alpha) ** 2)
        # Eq. 79
        c = (f / 16) * np.cos(alpha) ** 2 * (4 + f * (4 - 3 * np.cos(alpha) ** 2))
        # Eq. 80
        new_lon = omega + (1 - c) * f * np.sin(alpha) * (
                sigma + c * np.sin(sigma) * (cos_two_sigma_m + c * np.cos(sigma) * (-1 + 2 * (cos_two_sigma_m ** 2)))
        )
        delta_lon = new_lon - lon
        lon = new_lon
        if np.all(np.absolute(delta_lon) < 1e-12):
            break
    # Calculate the azimuth from point 1 to point 2
    azimuth1to2 = np.arctan2((np.cos(u2) * np.sin(lon)),
                             (np.cos(u1) * np.sin(u2) - np.sin(u1) * np.cos(u2) * np.cos(lon)))

    # add 2 pi in case an angle is below zero
    for azi in np.nditer(azimuth1to2, op_flags=['readwrite']):
        if azi < 0:
            azi[...] = azi + 2 * np.pi

    return np.round(azimuth1to2, 9)
