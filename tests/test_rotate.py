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

from earthkit.geo.rotate import rotate, rotate_wind, unrotate, unrotate_wind

from .testing import earthkit_test_data_file, normalise_lon


def _make_proj(name):
    from pyproj import CRS

    if name == "lambert_conformal_conic":
        grid_mapping = {
            "grid_mapping_name": "lambert_conformal_conic",
            "standard_parallel": [63.3, 63.3],
            "longitude_of_central_meridian": 15.0,
            "latitude_of_projection_origin": 63.3,
            "earth_radius": 6371000.0,
            "proj4": "+proj=lcc +lat_0=63.3 +lon_0=15 +lat_1=63.3 +lat_2=63.3 +no_defs +R=6.371e+06",
        }
    else:
        raise ValueError(f"Unknown projection: {name}")

    return CRS.from_cf(grid_mapping)


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


@pytest.mark.parametrize(
    "data,expected_result",
    [
        (
            # central meridian of the Lambert grid
            {
                "lats": np.array([63.3] * 5),
                "lons": np.array([15.0] * 5),
                "x_wind": np.array([0, 0, 1, -1, 1]),
                "y_wind": np.array([1, -1, 0, 0, 1]),
                "source_projection": _make_proj("lambert_conformal_conic"),
                "target_projection": "+proj=longlat",
            },
            (
                [0, 0, 1, -1, 1.00000023],
                [1, -1, -1.55647232e-07, -1.55647232e-07, 9.99999766e-01],
            ),
        ),
        (
            # south-west corner of the Lambert grid
            {
                "lats": np.array([50.31961636316951] * 5),
                "lons": np.array([0.2782806572089653] * 5),
                "x_wind": np.array([0, 0, 1, -1, 1]),
                "y_wind": np.array([1, -1, 0, 0, 1]),
                "source_projection": _make_proj("lambert_conformal_conic"),
                "target_projection": "+proj=longlat",
            },
            (
                [-0.22753458, 0.22753449, 0.97377004, -0.97376999, 0.74623563],
                [0.97377, -0.97377002, 0.22753443, -0.22753464, 1.20130445],
            ),
        ),
    ],
)
def test_rotate_wind(data, expected_result):
    res = rotate_wind(*data.values())
    assert np.allclose(res[0], expected_result[0])
    assert np.allclose(res[1], expected_result[1])


@pytest.mark.parametrize(
    "data,expected_result",
    [
        (
            {
                "lats": [70] * 3,
                "lons": [0] * 3,
                "raw_lats": [80] * 3,
                "raw_lons": [0] * 3,
                "x_wind": [1, -1, 0],
                "y_wind": [0, 0, 1],
                "south_pole_latitude": -80,
                "south_pole_longitude": 0,
            },
            ([1, -1, 0], [0, 0, 1]),
        ),
    ],
)
def test_unrotate_wind_points(data, expected_result):
    res = unrotate_wind(*data.values())
    assert np.allclose(res[0], expected_result[0])
    assert np.allclose(res[1], expected_result[1])


def test_unrotate_wind_global():
    """The data was generated by combining the results of the following
    ECMWF MARS retrievals:

        # global latlon grid
        param    = ["10u","10v"],
        levtype  = "sfc",
        date     = 20240514,
        grid     = [20, 20]

        # rotated global latlon grid
        param    = ["10u","10v"],
        levtype  = "sfc",
        date     = 20240514,
        rotation = [-20,-40],
        grid     = [20, 20]

    """
    data = np.load(earthkit_test_data_file("rotated_wind_20x20_input.npz"))
    ref = np.load(earthkit_test_data_file("rotated_wind_20x20_ref.npz"))

    res = unrotate_wind(
        data["lats_rot"],
        data["lons_rot"],
        data["lats_ori"],
        data["lons_ori"],
        data["x_wind_rot"],
        data["y_wind_rot"],
        data["rotation"][0],
        data["rotation"][1],
    )

    assert np.allclose(res[0], ref["x_wind"])
    assert np.allclose(res[1], ref["y_wind"])
