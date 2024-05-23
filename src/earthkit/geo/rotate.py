# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

from . import constants
from .coord import _normalise_longitude, latlon_to_xyz, xyz_to_latlon


def _normalise(x):
    return max(min(x, 1.0), -1.0)


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


def rotate_wind(lat, lon, x_wind, y_wind, source_projection, target_projection):
    """Rotate wind vectors from source to target projection coordinates.

    Parameters
    ----------
    lat: ndarray
        Latitude coordinates of the points (degrees).
    lon: ndarray
        Longitude coordinates of the points (degrees).
    x_wind: ndarray
        x-component of the wind in the source projection (m/s).
    y_wind: ndarray
        y-component of the wind in the source projection (m/s).
    source_projection: str, dict
        Projection that the wind components are defined in. It is either a proj string
        or a dict of map projection control parameter key/value pairs.
    target_projection: str, dict
        Projection that the wind components should be transformed to. It is either a proj string
        or a dict of map projection control parameter key/value pairs.

    Returns
    -------
    ndarray
        x-component of the wind in the target projection (m/s).
    ndarray
        y-component of the wind in the target projection (m/s).


    The wind vector is returned at the same locations (``lat``, ``lon``), but rotated
    into ``target_projection`` coordinates. The magnitude of the wind is preserved.
    """
    import numpy as np
    import pyproj

    if source_projection == target_projection:
        return x_wind, x_wind

    source_projection = pyproj.Proj(source_projection)
    target_projection = pyproj.Proj(target_projection)

    transformer = pyproj.transformer.Transformer.from_proj(
        source_projection, target_projection
    )

    # To compute the new vector components:
    # 1) perturb each position in the direction of the winds
    # 2) convert the perturbed positions into the new coordinate system
    # 3) measure the new x/y components.
    #
    # A complication occurs when using the longlat "projections", since this is not a cartesian grid
    # (i.e. distances in each direction is not consistent), we need to deal with the fact that the
    # width of a longitude varies with latitude
    orig_speed = np.sqrt(x_wind**2 + y_wind**2)

    x0, y0 = source_projection(lon, lat)

    if source_projection.name != "longlat":
        x1 = x0 + x_wind
        y1 = y0 + y_wind
    else:
        # Reduce the perturbation, since x_wind and y_wind are in meters, which would create
        # large perturbations in lat, lon. Also, deal with the fact that the width of longitude
        # varies with latitude.
        factor = 3600000.0
        x1 = x0 + x_wind / factor / np.cos(np.deg2rad(lat))
        y1 = y0 + y_wind / factor

    X0, Y0 = transformer.transform(x0, y0)
    X1, Y1 = transformer.transform(x1, y1)

    new_x_wind = X1 - X0
    new_y_wind = Y1 - Y0
    if target_projection.name == "longlat":
        new_x_wind *= np.cos(np.deg2rad(lat))

    if target_projection.name == "longlat" or source_projection.name == "longlat":
        # Ensure the wind speed is not changed (which might not be the case since the
        # units in longlat is degrees, not meters)
        curr_speed = np.sqrt(new_x_wind**2 + new_y_wind**2)
        new_x_wind *= orig_speed / curr_speed
        new_y_wind *= orig_speed / curr_speed

    return new_x_wind, new_y_wind


def unrotate_wind(
    lat,
    lon,
    lat_unrotated,
    lon_unrotated,
    x_wind,
    y_wind,
    south_pole_latitude,
    south_pole_longitude,
    south_pole_rotation_angle=0,
):
    """
    Rotate wind vectors on a rotated grid back into eastward and northward components.

    Parameters
    ----------
    lat: ndarray
        Latitude coordinates of the rotated points (degrees).
    lon: ndarray
        Longitude coordinates of the rotated points (degrees).
    lat_unrotated: ndarray
        Latitude coordinates of the points before rotation (degrees).
    lon_unrotated: ndarray
        Longitude coordinates of the points before rotation (degrees).
    x_wind: ndarray
        x-component of the wind in the rotated
        coordinate system at the rotated points (m/s).
    y_wind: ndarray
        y-component of the wind in the rotated
        coordinate system at the rotated points (m/s).
    south_pole_latitude: float
        Latitude of the south pole defining the rotation (degrees).
    south_pole_longitude: float
        Longitude of the south pole defining the rotation (degrees).
    south_pole_rotation_angle: float, optional
        Rotation angle around the south pole (degrees). Currently not supported.

    Returns
    -------
    ndarray
        x-component of the wind vector rotated back to
        eastward and northward components (m/s).
    ndarray
        y-component of the wind vector rotated back to
        eastward and northward components (m/s).


    When a grid is rotated spherically the wind components are locally rotated
    at each target grid point into the new local (rotated) coordinate
    system. :func:`unrotate_wind` performs the inverse operation, and rotates back the wind
    vectors to their original directions at each point. The wind vector is returned at
    the same locations (``lat``, ``lon``) and its magnitude is preserved.

    """
    import numpy as np

    assert south_pole_rotation_angle == 0

    C = np.deg2rad(constants.NORTH - south_pole_latitude)
    cos_C = np.cos(C)
    sin_C = np.sin(C)

    new_x = np.zeros_like(x_wind)
    new_y = np.zeros_like(y_wind)

    for i, (vx, vy, lon_r, lon_ur) in enumerate(
        zip(x_wind, y_wind, lon, lon_unrotated)
    ):
        lon_r = south_pole_longitude - lon_r
        lon_r = _normalise_longitude(lon_r, -constants.STRAIGHT_ANGLE)

        a = np.deg2rad(lon_r)
        b = np.deg2rad(lon_ur)
        q = 1.0 if (sin_C * lon_r < 0.0) else -1.0  # correct quadrant

        cos_c = _normalise(np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * cos_C)
        sin_c = q * np.sqrt(1.0 - cos_c * cos_c)

        new_x[i] = cos_c * vx + sin_c * vy
        new_y[i] = -sin_c * vx + cos_c * vy

    return new_x, new_y
