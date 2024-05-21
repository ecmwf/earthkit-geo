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

from earthkit.geo.rotate import rotate, unrotate
from earthkit.geo.wind import _normalise_longitude


def normalise_lon(x):
    return np.array([_normalise_longitude(lon, -180) for lon in x])


@pytest.mark.parametrize(
    "data,expected_result",
    [
        (
            {
                "lat": np.array([-89, -45, 0, 45, 89]),
                "lon": np.array([0] * 5),
                "south_pole_lat": -20,
                "south_pole_lon": -40,
            },
            ([-19.0, 25.0, 70.0, 65.0, 21.0], [-40.0, -40.0, -40.0, 140.0, 140.0]),
        ),
        (
            {
                "lat": np.array([-89] * 10),
                "lon": np.linspace(-180, 180, 10),
                "south_pole_lat": -20,
                "south_pole_lon": -40,
            },
            (
                [
                    -21.0,
                    -20.76470947,
                    -20.17055585,
                    -19.49764434,
                    -19.05994355,
                    -19.05994355,
                    -19.49764434,
                    -20.17055585,
                    -20.76470947,
                    -21.0,
                ],
                [
                    -40.0,
                    -40.68742246,
                    -41.04915724,
                    -40.91870127,
                    -40.36184216,
                    -39.63815784,
                    -39.08129873,
                    -38.95084276,
                    -39.31257754,
                    -40.0,
                ],
            ),
        ),
    ],
)
def test_rotate_points(data, expected_result):
    lat_r, lon_r = rotate(*data.values())
    assert np.allclose(lat_r, expected_result[0], atol=1e-5)
    assert np.allclose(
        normalise_lon(lon_r), normalise_lon(expected_result[1]), atol=1e-5
    )

    lat_ur, lon_ur = unrotate(
        lat_r, lon_r, data["south_pole_lat"], data["south_pole_lon"]
    )
    assert np.allclose(lat_ur, data["lat"], atol=1e-5)
    assert np.allclose(normalise_lon(lon_ur), normalise_lon(data["lon"]), atol=1e-5)
