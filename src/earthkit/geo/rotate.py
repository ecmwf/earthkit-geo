# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


def rotate(lat, lon, south_pole_lat, south_pole_lon):
    """
    Rotate geographical coordinates with respect to a south pole.

    Parameters
    ----------
    lat: ndarray
        Latitudes (degrees).
    lon: ndarray
        Longitudes (degrees).
    south_pole_lat: float
        Latitude of the south pole defining the rotation (degrees).
    south_pole_lon: float
        Longitude of the south pole defining the rotation (degrees).

    Returns
    -------
    ndarray
        Rotated latitudes (degrees).
    ndarray
        Rotated longitudes (degrees).

    """
    import numpy as np

    def from_xyz(x, y, z):
        return (
            np.rad2deg(np.arcsin(np.minimum(1.0, np.maximum(-1.0, z)))),
            np.rad2deg(np.arctan2(y, x)),
        )

    def to_xyz(lat, lon):
        lam = np.deg2rad(lon)
        phi = np.deg2rad(lat)

        sp = np.sin(phi)
        cp = np.cos(phi)
        sl = np.sin(lam)
        cl = np.cos(lam)

        return (cp * cl, cp * sl, sp)

    theta = -np.deg2rad(south_pole_lat + 90.0)
    phi = -np.deg2rad(south_pole_lon)

    ct = np.cos(theta)
    st = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    matrix = np.array(
        [[cp * ct, sp, cp * st], [-ct * sp, cp, -sp * st], [-st, 0.0, ct]]
    )

    return from_xyz(*np.dot(matrix, to_xyz(lat, lon)))


def unrotate(lat, lon, south_pole_lat, south_pole_lon):
    """
    Unrotate geographical coordinates with respect to a south pole.

    Parameters
    ----------
    lat: ndarray
        Latitudes (degrees).
    lon: ndarray
        Longitudes (degrees).
    south_pole_lat: float
        Latitude of the south pole defining the rotation (degrees).
    south_pole_lon: float
        Longitude of the south pole defining the rotation (degrees).

    Returns
    -------
    ndarray
        Unrotated latitudes (degrees).
    ndarray
        Unrotated longitudes (degrees).

    """
    import numpy as np

    def from_xyz(x, y, z):
        return (
            np.rad2deg(np.arcsin(np.minimum(1.0, np.maximum(-1.0, z)))),
            np.rad2deg(np.arctan2(y, x)),
        )

    def to_xyz(lat, lon):
        lam = np.deg2rad(lon)
        phi = np.deg2rad(lat)

        sp = np.sin(phi)
        cp = np.cos(phi)
        sl = np.sin(lam)
        cl = np.cos(lam)

        return (cp * cl, cp * sl, sp)

    theta = -np.deg2rad(south_pole_lat + 90.0)
    phi = -np.deg2rad(south_pole_lon)

    ct = np.cos(theta)
    st = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    matrix = np.array(
        [[cp * ct, sp, cp * st], [-ct * sp, cp, -sp * st], [-st, 0.0, ct]]
    )

    matrix = matrix.T
    return from_xyz(*np.dot(matrix, to_xyz(lat, lon)))
