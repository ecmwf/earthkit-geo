# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

from . import constants


def generate_latlon(dx, dy):
    """Generate latitude and longitude arrays for a global regular latitude-longitude grid.

    Parameters
    ----------
    dx: number
        Increment in longitudes (degrees)
    dy: number
        Increment in latitudes (degrees)

    Returns
    -------
    ndarray
        Latitudes (degrees) with the shape of the grid
    ndarray
        Longitudes (degrees) with the shape of the grid


    The global regular latitude-longitude grid is formed as follows:

    - the points start at the north pole and go down to the south pole from west to east in each latitude band
    - the shape of the grid is (ny, nx) where ny is the number of latitudes and nx is the number of longitudes

    Examples
    --------
    >>> import earthkit.geo
    >>> lat, lon = earthkit.geo.generate_latlon(120, 45)
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
    import numpy as np

    lat_v = np.linspace(
        constants.north,
        constants.south,
        int((constants.north - constants.south) / dy) + 1,
    )
    lon_v = np.linspace(0, constants.full_circle - dx, int(constants.full_circle / dx))
    lon, lat = np.meshgrid(lon_v, lat_v)
    return lat, lon
