'''
Collection of geodesy/pyproj algorithms for PyRate.

This module depends on PyProj/PROJ4 to replace llh2local.m in MATLAB Pirate.

Created: 13/3/13

.. codeauthor:: Ben Davies
'''


from math import floor

import pyproj


def utm_zone(longitude):
    '''
    Returns basic UTM zone for given longitude in degrees. Currently does NOT
    handle the sub-zoning around Scandanavian countries.
    See http://www.dmap.co.uk/utmworld.htm
    '''

    if longitude == 180:
        return 60.0
    return floor((longitude + 180) / 6.0) + 1


def cell_size(lat, lon, x_step, y_step):
    '''
    Converts X|Y_STEP in degrees to X & Y cell length/width in metres.
    lat: latitude in degrees
    lon: longitude in degrees
    x_step: horizontal step size in degrees
    y_step: vertical step size in degrees
    '''

    if lat > 84.0 or lat < -80:
        msg = "No UTM zone for polar region: > 84 degrees N or < 80 degrees S"
        raise ValueError(msg)

    zone = utm_zone(lon)
    p0 = pyproj.Proj(proj='latlong', ellps='WGS84')
    p1 = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')
    assert p0.is_latlong()
    assert not p1.is_latlong()

    x0, y0 = pyproj.transform(p0, p1, lon, lat)
    x1, y1 = pyproj.transform(p0, p1, lon + x_step, lat + y_step)
    return tuple(abs(e) for e in (x1 - x0, y1 - y0))
