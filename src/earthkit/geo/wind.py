# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


def _normalise(x):
    return max(min(x, 1.0), -1.0)


def _normalise_longitude(lon, minimum):
    while lon < minimum:
        lon += 360

    while lon >= minimum + 360:
        lon -= 360

    return lon


def rotate_winds(lats, lons, x_wind, y_wind, source_projection, target_projection):
    """Code provided by MetNO"""
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
        # Ensure the wind speed is not changed (which might not the case since the units in longlat
        # is degrees, not meters)
        curr_speed = np.sqrt(new_x_wind**2 + new_y_wind**2)
        new_x_wind *= orig_speed / curr_speed
        new_y_wind *= orig_speed / curr_speed

    return new_x_wind, new_y_wind


def unrotate_winds(
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
    import numpy as np

    # Code from MIR
    assert south_pole_rotation_angle == 0
    C = np.deg2rad(90 - south_pole_latitude)
    cos_C = np.cos(C)
    sin_C = np.sin(C)

    new_x = np.zeros_like(x_wind)
    new_y = np.zeros_like(y_wind)

    for i, (vx, vy, lat, lon, raw_lat, raw_lon) in enumerate(
        zip(x_wind, y_wind, lats, lons, raw_lats, raw_lons)
    ):
        lonRotated = south_pole_longitude - lon
        lon_rotated = _normalise_longitude(lonRotated, -180)
        lon_unrotated = raw_lon

        a = np.deg2rad(lon_rotated)
        b = np.deg2rad(lon_unrotated)
        q = 1 if (sin_C * lon_rotated < 0.0) else -1.0  # correct quadrant

        cos_c = _normalise(np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * cos_C)
        sin_c = q * np.sqrt(1.0 - cos_c * cos_c)

        new_x[i] = cos_c * vx + sin_c * vy
        new_y[i] = -sin_c * vx + cos_c * vy

    return new_x, new_y
