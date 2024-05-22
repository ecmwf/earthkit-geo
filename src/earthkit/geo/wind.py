# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

from . import constants


def _normalise(x):
    return max(min(x, 1.0), -1.0)


def _normalise_longitude(lon, minimum):
    while lon < minimum:
        lon += constants.FULL_ANGLE

    while lon >= minimum + constants.FULL_ANGLE:
        lon -= constants.FULL_ANGLE

    return lon


def rotate_wind(lats, lons, x_wind, y_wind, source_projection, target_projection):
    """Rotate wind vectors from source to target projection coordinates.

    Parameters
    ----------
    lats: ndarray
        Latitude coordinates of the points (degrees).
    lons: ndarray
        Longitude coordinates of the points (degrees).
    x_wind: ndarray
        x-component of the wind in the source projection (m/s).
    y_wind: ndarray
        y-component of the wind in the source projection (m/s).
    source_projection: str
        Projection that the wind components are defined in.
    target_projection: str
        Projection that the wind components should be transformed to.

    Returns
    -------
    ndarray
        x-component of the wind in the target projection (m/s).
    ndarray
        y-component of the wind in the target projection (m/s).


    The wind vector is returned at the same locations (``lats``, ``lons``), but rotated
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

    x0, y0 = source_projection(lons, lats)

    if source_projection.name != "longlat":
        x1 = x0 + x_wind
        y1 = y0 + y_wind
    else:
        # Reduce the perturbation, since x_wind and y_wind are in meters, which would create
        # large perturbations in lat, lon. Also, deal with the fact that the width of longitude
        # varies with latitude.
        factor = 3600000.0
        x1 = x0 + x_wind / factor / np.cos(np.deg2rad(lats))
        y1 = y0 + y_wind / factor

    X0, Y0 = transformer.transform(x0, y0)
    X1, Y1 = transformer.transform(x1, y1)

    new_x_wind = X1 - X0
    new_y_wind = Y1 - Y0
    if target_projection.name == "longlat":
        new_x_wind *= np.cos(np.deg2rad(lats))

    if target_projection.name == "longlat" or source_projection.name == "longlat":
        # Ensure the wind speed is not changed (which might not be the case since the
        # units in longlat is degrees, not meters)
        curr_speed = np.sqrt(new_x_wind**2 + new_y_wind**2)
        new_x_wind *= orig_speed / curr_speed
        new_y_wind *= orig_speed / curr_speed

    return new_x_wind, new_y_wind


def unrotate_wind(
    lats,
    lons,
    raw_lats,
    raw_lons,
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
    lats: ndarray
        Latitude coordinates of the rotated points (degrees).
    lons: ndarray
        Longitude coordinates of the rotated points (degrees).
    raw_lats: ndarray
        Latitude coordinates of the points before rotation (degrees).
    raw_lons: ndarray
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
    the same locations (``lats``, ``lons``) and its magnitude is preserved.

    """
    import numpy as np

    assert south_pole_rotation_angle == 0

    C = np.deg2rad(constants.NORTH - south_pole_latitude)
    cos_C = np.cos(C)
    sin_C = np.sin(C)

    new_x = np.zeros_like(x_wind)
    new_y = np.zeros_like(y_wind)

    for i, (vx, vy, lat, lon, raw_lat, raw_lon) in enumerate(
        zip(x_wind, y_wind, lats, lons, raw_lats, raw_lons)
    ):
        lonRotated = south_pole_longitude - lon
        lon_rotated = _normalise_longitude(lonRotated, -constants.STRAIGHT_ANGLE)
        lon_unrotated = raw_lon

        a = np.deg2rad(lon_rotated)
        b = np.deg2rad(lon_unrotated)
        q = 1 if (sin_C * lon_rotated < 0.0) else -1.0  # correct quadrant

        cos_c = _normalise(np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * cos_C)
        sin_c = q * np.sqrt(1.0 - cos_c * cos_c)

        new_x[i] = cos_c * vx + sin_c * vy
        new_y[i] = -sin_c * vx + cos_c * vy

    return new_x, new_y
