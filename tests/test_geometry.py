import numpy as np
from pyrate.core.geometry import get_lonlat_coords, get_lonlat_coords_slow, get_and_write_radar_coords


def test_get_lonlat_coords_vectorised(dem):
    lon, lat = get_lonlat_coords_slow(dem)
    lon_v, lat_v = get_lonlat_coords(dem)
    np.testing.assert_array_almost_equal(lon, lon_v)
    np.testing.assert_array_almost_equal(lat, lat_v)
