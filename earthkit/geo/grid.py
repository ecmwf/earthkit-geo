# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


from . import constants


def to_latlon(grid):
    """Generate latitude and longitude arrays for a given ``grid``.

    Parameters
    ----------
    grid: dict
        Grid specification. At the moment, only global regular latitude-longitude
        grids are supported. The required format is as follows::

        {"grid": [dlon, dlat]}

    Returns
    -------
    ndarray
        Latitudes (degrees) with the shape of the grid
    ndarray
        Longitudes (degrees) with the shape of the grid


    The global regular latitude-longitude grid is defined as follows:

    - the points start at the north pole and go down to the south pole from
      west (0°) to east in each latitude band
    - the **shape** of the grid is (ny, nx) where ny is the number of latitudes
      and nx is the number of longitudes

    Examples
    --------
    >>> import earthkit.geo
    >>> lat, lon = earthkit.geo.grid.to_latlon({"grid": [120, 45]})
    >>> lat
    array([[ 90.,  90.,  90.],
        [ 45.,  45.,  45.],
        [  0.,   0.,   0.],
        [-45., -45., -45.],
        [-90., -90., -90.]])
    >>> lon
    array([[  0., 120., 240.],
        [  0., 120., 240.],
        [  0., 120., 240.],
        [  0., 120., 240.],
        [  0., 120., 240.]])

    """
    if not isinstance(grid, dict):
        raise TypeError(f"Invalid grid type {type(grid)}, must be a dict.")

    grid = dict(grid)
    dlon, dlat = grid.pop("grid")
    if grid:
        raise KeyError(f"Unexpected keys in grid = {list(grid.keys())}")

    if dlat <= 0:
        raise ValueError(f"Increment in latitude must be positive: {dlat=}")
    if dlon <= 0:
        raise ValueError(f"Increment in longitude must be positive: {dlon=}")

    import numpy as np

    eps_lat = dlat / 1e2
    eps_lon = dlon / 1e2

    north = constants.north
    south = constants.south
    west = 0.0
    east = constants.full_circle

    # north-south
    ny2 = round(constants.quarter_circle / dlat)
    half_lat_range = ny2 * dlat
    delta = constants.quarter_circle - half_lat_range
    if delta < -eps_lat:
        ny2 -= 1
        half_lat_range = ny2 * dlat
        delta = constants.quarter_circle - half_lat_range

    ny = ny2 * 2 + 1

    if abs(delta) < eps_lat:
        north = constants.north
        south = constants.south
    else:
        north = half_lat_range
        south = -north

    assert ny % 2 == 1
    lats = np.array([north - i * dlat for i in range(ny)])
    lats[-1] = south

    # west-east
    west = 0.0
    east = constants.full_circle

    nx = round(constants.full_circle / dlon)
    lon_range = nx * dlon
    delta = constants.full_circle - lon_range
    if abs(delta) < eps_lon or delta < 0:
        nx = nx - 1
    nx = nx + 1
    east = (nx - 1) * dlon

    lons = np.array([west + i * dlon for i in range(nx)])
    lons[-1] = east

    lon, lat = np.meshgrid(lons, lats)
    return lat, lon
