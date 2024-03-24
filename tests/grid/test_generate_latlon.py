#!/usr/bin/env python3

# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import numpy as np
import pytest

from earthkit.geo import generate_latlon


@pytest.mark.parametrize(
    "dx,dy,shape,vals_x,vals_y",
    [
        (20, 20, (10, 18), np.arange(0, 360, 20), np.arange(90, -90 - 20, -20)),
        (30, 20, (10, 12), np.arange(0, 360, 30), np.arange(90, -90 - 20, -20)),
    ],
)
def test_generate_latlon(dx, dy, shape, vals_x, vals_y):
    lat, lon = generate_latlon(dx, dy)
    assert lat.shape == shape
    assert lon.shape == shape

    ny = lat.shape[0]
    nx = lat.shape[1]
    for i in range(ny):
        first = lat[i, 0]
        assert np.allclose(lat[i, :], np.ones(nx) * first)

    for i in range(nx):
        first = lon[0, i]
        assert np.allclose(lon[:, i], np.ones(ny) * first)

    assert np.allclose(lat[:, 0], vals_y)
    assert np.allclose(lon[0, :], vals_x)
