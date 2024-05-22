# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

from .coord import latlon_to_xyz, xyz_to_latlon


def _rotation_matrix(south_pole_lat, south_pole_lon):
    """
    Define matrix for spherical rotation.
    """
    import numpy as np

    theta = -np.deg2rad(south_pole_lat + 90.0)
    phi = -np.deg2rad(south_pole_lon)

    ct = np.cos(theta)
    st = np.sin(theta)
    cp = np.cos(phi)
    sp = np.sin(phi)

    return np.array([[cp * ct, sp, cp * st], [-ct * sp, cp, -sp * st], [-st, 0.0, ct]])


def rotate(lat, lon, south_pole_lat, south_pole_lon):
    """
    Rotate geographical coordinates on a sphere.

    Parameters
    ----------
    lat: ndarray
        Latitudes (degrees).
    lon: ndarray
        Longitudes (degrees).
    south_pole_lat: float
        Latitude of the rotated south pole (degrees).
    south_pole_lon: float
        Longitude of the rotated south pole (degrees).

    Returns
    -------
    ndarray
        Rotated latitudes (degrees).
    ndarray
        Rotated longitudes (degrees).


    The rotation is specified by the position
    where the south pole is rotated to, i.e. by (``south_pole_lat``, ``south_pole_lon``).


    See Also
    --------
    :obj:`unrotate`

    Examples
    --------
    >>> from earthkit.geo.rotate import rotate
    >>> lat = [-90]
    >>> lon = [0]
    >>> south_pole_lat = -20
    >>> south_pole_lon = -40
    >>> rotate(lat, lon, south_pole_lat, south_pole_lon)
    (array([-20.]), array([-40.]))
    """
    import numpy as np

    matrix = _rotation_matrix(south_pole_lat, south_pole_lon)
    return xyz_to_latlon(*np.dot(matrix, latlon_to_xyz(lat, lon)))


def unrotate(lat, lon, south_pole_lat, south_pole_lon):
    """
    Unrotate geographical coordinates on a sphere.

    Parameters
    ----------
    lat: ndarray
        Latitudes of the rotated points (degrees).
    lon: ndarray
        Longitudes of the rotated points (degrees).
    south_pole_lat: float
        Latitude of the rotated south pole (degrees).
    south_pole_lon: float
        Longitude of the rotated south pole (degrees).

    Returns
    -------
    ndarray
        Unrotated latitudes (degrees).
    ndarray
        Unrotated longitudes (degrees).


    :func:`unrotate` operates on points already rotated with a rotation specified
    by (``south_pole_lat``, ``south_pole_lon``). The function rotates the points
    back to the original locations.

    See Also
    --------
    :obj:`rotate`

    Examples
    --------
    >>> from earthkit.geo.rotate import unrotate
    >>> lat = [-18]
    >>> lon = [-41]
    >>> south_pole_lat = -20
    >>> south_pole_lon = -40
    >>> unrotate(lat, lon, south_pole_lat, south_pole_lon)
    (array([-87.78778846]), array([-25.4673765]))
    """
    import numpy as np

    matrix = _rotation_matrix(south_pole_lat, south_pole_lon)
    matrix = matrix.T
    return xyz_to_latlon(*np.dot(matrix, latlon_to_xyz(lat, lon)))
