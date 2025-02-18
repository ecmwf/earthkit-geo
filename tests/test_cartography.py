# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


import pytest
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from earthkit.geo.cartography import _closest_resolution
from earthkit.geo.cartography import country_polygons
from earthkit.geo.cartography import multipolygon_to_coordinates


@pytest.mark.parametrize(
    "distance_m, expected_resolution",
    [
        (0, "10m"),  # Closest to 10m
        (20e6, "10m"),  # Closest to 10m
        (80e6, "110m"),  # Closest to 110m
        (79e6, "50m"),  # Closest to 50m
        (45e6, "50m"),  # Closest to 50m
        (110e6, "110m"),  # Exactly 110m
    ],
)
def test_closest_resolution(distance_m, expected_resolution):
    assert _closest_resolution(distance_m) == expected_resolution


def test_multipolygon_to_coordinates():
    poly1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    poly2 = Polygon([(2, 2), (3, 2), (3, 3), (2, 3), (2, 2)])
    multipolygon = MultiPolygon([poly1, poly2])

    coordinates = multipolygon_to_coordinates(multipolygon)

    expected_coordinates = [
        [[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]],
        [[2, 2], [2, 3], [3, 3], [3, 2], [2, 2]],
    ]

    assert coordinates == expected_coordinates


def test_country_polygons_invalid_country():
    with pytest.raises(ValueError):
        country_polygons("InvalidCountryName")


def test_country_polygons_single_country():
    coordinates = country_polygons("France")

    # Verify the output format
    assert isinstance(coordinates, list)
    assert all(isinstance(polygon, list) for polygon in coordinates)
    assert all(isinstance(coord, list) and len(coord) == 2 for polygon in coordinates for coord in polygon)


def test_country_polygons_multiple_countries():
    coordinates = country_polygons(["Poland", "Germany"])

    # Verify the output format
    assert isinstance(coordinates, list)
    assert all(isinstance(polygon, list) for polygon in coordinates)
    assert all(isinstance(coord, list) and len(coord) == 2 for polygon in coordinates for coord in polygon)


def test_country_polygons_output_format():
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    multipolygon = MultiPolygon([poly])

    coordinates = multipolygon_to_coordinates(multipolygon)

    # Verify the output format
    assert isinstance(coordinates, list)
    assert all(isinstance(polygon, list) for polygon in coordinates)
    assert all(isinstance(coord, list) and len(coord) == 2 for polygon in coordinates for coord in polygon)
