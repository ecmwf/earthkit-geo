#!/usr/bin/env python3

# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import os

import numpy as np
import pytest

from earthkit.geo.grid import to_latlon

DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def file_in_testdir(filename):
    return os.path.join(DATA_PATH, filename)


def _global_grids():
    import json

    with open(file_in_testdir("global_ll_ref.json"), "r") as f:
        d = json.load(f)

    for _, v in d.items():
        yield v


# these are all the cases for the earthkit-regrid target latlon grids
@pytest.mark.parametrize("grid", _global_grids())
def test_to_latlon_1(grid):
    dx, dy = grid["grid"]
    lat, lon = to_latlon({"grid": [dx, dy]})

    assert list(lat.shape) == grid["shape"]
    assert list(lon.shape) == grid["shape"]

    # bounding box
    north = lat[0, 0]
    south = lat[-1, 0]
    west = lon[0, 0]
    east = lon[0, -1]
    ref_north, ref_west, ref_south, ref_east = grid["area"]
    assert np.isclose(north, ref_north)
    assert np.isclose(south, ref_south)
    assert np.isclose(west, ref_west)
    assert np.isclose(east, ref_east)

    # equator
    ny = lat.shape[0]
    assert ny % 2 == 1
    assert np.isclose(lat[int((ny - 1) / 2), 0], 0.0)


@pytest.mark.parametrize(
    "dx,dy,shape,vals_x,vals_y",
    [
        (10, 10, (19, 36), np.arange(0, 360, 10), np.arange(90, -90 - 10, -10)),
        (30, 20, (9, 12), np.arange(0, 360, 30), np.arange(80, -80 - 20, -20)),
    ],
)
def test_to_latlon_2(dx, dy, shape, vals_x, vals_y):
    gs = {"grid": [dx, dy]}
    lat, lon = to_latlon(gs)
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


@pytest.mark.parametrize(
    "grid,error",
    [
        (1, TypeError),
        ({"grid": 1}, TypeError),
        ({"grid": [-1, 1]}, ValueError),
        ({"grid": [1, 1], "area": [70, 0, 40, 25]}, KeyError),
    ],
)
def test_to_latlon_invalid(grid, error):
    with pytest.raises(error):
        to_latlon(grid)
