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

from earthkit.geo.coord import latlon_to_xyz
from earthkit.geo.coord import xyz_to_latlon

from .testing import normalise_lon

R = 1.0
L = R * np.sqrt(2.0) / 2

# latlon -> xyz
REF_DATA = [
    ((90, 0), (0, 0, R)),
    ((-90, 0), (0, 0, -R)),
    ((0, 0), (R, 0, 0)),
    ((0, -360), (R, 0, 0)),
    ((0, 90), (0, R, 0)),
    ((0, -270), (0, R, 0)),
    ((0, 180), (-R, 0, 0)),
    ((0, -180), (-R, 0, 0)),
    ((0, 270), (0, -R, 0)),
    ((0, -90), (0, -R, 0)),
    ((0, 45), (L, L, 0)),
    ((0, -315), (L, L, 0)),
    ((0, 135), (-L, L, 0)),
    ((0, -225), (-L, L, 0)),
    ((0, 225), (-L, -L, 0)),
    ((0, -135), (-L, -L, 0)),
    ((0, 315), (L, -L, 0)),
    ((0, -45), (L, -L, 0)),
]


def _make_latlon_input():
    for r in REF_DATA:
        yield r[0], r[1]


def _make_xyz_input():
    for r in REF_DATA:
        yield r[1], r[0]


@pytest.mark.parametrize("latlon,expected_result", _make_latlon_input())
def test_latlon_to_xyz(latlon, expected_result):
    res = latlon_to_xyz(*latlon)
    assert np.allclose(res, expected_result), f"latlon={latlon}"


@pytest.mark.parametrize("xyz,expected_result", _make_xyz_input())
def test_xyz_to_latlon(xyz, expected_result):
    res = xyz_to_latlon(*xyz)
    assert np.allclose(res[0], expected_result[0]), f"xyz={xyz}"
    assert np.allclose(normalise_lon(res[1]), normalise_lon(expected_result[1])), f"xyz={xyz}"
