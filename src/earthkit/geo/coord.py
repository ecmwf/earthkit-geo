import numpy as np


def _xyz_to_latlon(x, y, z):
    return (
        np.rad2deg(np.arcsin(np.minimum(1.0, np.maximum(-1.0, z)))),
        np.rad2deg(np.arctan2(y, x)),
    )


def _latlon_to_xyz(lat, lon):
    """Works on the unit sphere."""
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    lat = np.radians(lat)
    lon = np.radians(lon)

    s_lon = np.sin(lon)
    c_lon = np.cos(lon)
    s_lat = np.sin(lat)
    c_lat = np.cos(lat)

    return c_lat * c_lon, c_lat * s_lon, s_lat
