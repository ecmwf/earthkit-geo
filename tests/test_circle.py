#!/usr/bin/env python3

# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import decimal

import numpy as np
import pytest

from earthkit.geo.coord import latlon_to_xyz
from earthkit.geo.distance import haversine_distance
from earthkit.geo.figure import IFS_SPHERE, UNIT_SPHERE
from earthkit.geo.polyline import circle, polyline_normal

R = IFS_SPHERE.radius  # m

# radii (km): 1 m, 1 km, 1000 km, 10000 km, 20000 km and 1 m short of the antipode
RADII_KM = [0.001, 1.0, 1000.0, 10000.0, 20000.0, (np.pi * R - 1.0) / 1000.0]

# angular radii (rad) for the unit sphere (radius in km = 0.001 * angle)
UNIT_ANGLES = [1e-6, 0.1, 1.0, 3.0, np.pi - 1e-6]


def _random_centres():
    rng = np.random.default_rng(20261022)
    lats = rng.uniform(-90, 90, 30)
    lons = rng.uniform(-180, 180, 30)
    return [(float(a), float(b)) for a, b in zip(lats, lons)]


CENTRES = [
    # poles and next to them
    (90.0, 0.0),
    (-90.0, 0.0),
    (90.0, 45.0),
    (-90.0, -120.0),
    (89.999999, 10.0),
    (-89.9999, -170.0),
    # next to the dateline
    (0.0, 179.99),
    (0.0, -179.99),
    (45.0, -179.5),
    (-45.0, 179.5),
    (30.0, 180.0),
    (-30.0, -180.0),
    # others
    (0.0, 0.0),
    (51.5, -0.1),
    (-33.9, 151.2),
] + _random_centres()


def _vec(lat, lon):
    return np.stack(latlon_to_xyz(lat, lon), axis=-1)


def _angle(a, b):
    """Independent angular distance (rad) between the unit vectors ``a`` and ``b`` (last axis of size 3)."""
    return np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1))


def _angdiff(a, b):
    """Absolute difference of the angles (rad) modulo a full turn."""
    return np.abs((np.asarray(a) - np.asarray(b) + np.pi) % (2 * np.pi) - np.pi)


def _delta_deg(radius_km, figure=IFS_SPHERE):
    return float(np.rad2deg(radius_km * 1000.0 / figure.radius))


def _near_pole(lat, radius_km, limit=1.0):
    """True if the circle itself passes within ``limit`` degrees of a pole.

    The distance from the circle to the north (south) pole is the difference between the angular radius and the
    angular distance from the centre to that pole, which is ``90 - lat`` (``90 + lat``).
    """
    d = _delta_deg(radius_km)
    return abs(d - (90 - lat)) <= limit or abs(d - (90 + lat)) <= limit


# ---------------------------------------------------------------------------
# structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [4, 5, 7, 100, 361])
def test_circle_length_and_closed(n):
    lats, lons = circle(51.5, -0.1, 500, n=n)
    for a in (lats, lons):
        assert isinstance(a, np.ndarray)
        assert a.shape == (n,)
        assert a.dtype == np.float64
    assert lats[-1] == lats[0]
    assert lons[-1] == lons[0]


def test_circle_default_n():
    lats, lons = circle(10, 20, 100)
    assert lats.shape == (100,)
    assert lons.shape == (100,)


@pytest.mark.parametrize("radius", RADII_KM)
@pytest.mark.parametrize("centre", [(0.0, 0.0), (51.5, -0.1), (90.0, 0.0), (-89.9, 100.0), (-30.0, 179.9)])
def test_circle_distinct_points(centre, radius):
    lats, lons = circle(*centre, radius, n=100)
    pts = set(zip(lats[:-1].tolist(), lons[:-1].tolist()))
    assert len(pts) == 99


# ---------------------------------------------------------------------------
# distance from the centre
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("radius", RADII_KM)
def test_circle_distance_to_centre(radius):
    for lat, lon in CENTRES:
        lats, lons = circle(lat, lon, radius, n=100)
        c = _vec(lat, lon)
        p = _vec(lats, lons)
        dist = _angle(c, p) * R
        np.testing.assert_allclose(
            dist, radius * 1000.0, rtol=1e-9, atol=1e-6, err_msg=f"centre=({lat}, {lon}) radius={radius}"
        )


@pytest.mark.parametrize("angle", UNIT_ANGLES)
def test_circle_distance_to_centre_unit_sphere(angle):
    radius = 0.001 * angle
    for lat, lon in CENTRES:
        lats, lons = circle(lat, lon, radius, n=100, figure=UNIT_SPHERE)
        c = _vec(lat, lon)
        p = _vec(lats, lons)
        dist = _angle(c, p) * UNIT_SPHERE.radius
        np.testing.assert_allclose(dist, angle, rtol=1e-9, atol=1e-6, err_msg=f"centre=({lat}, {lon}) angle={angle}")


def test_circle_figure_gives_the_same_circle_for_the_same_angle():
    # an angular radius of 0.1 rad on both spheres: radius km = 0.1 * radius_m / 1000
    lats1, lons1 = circle(10, 20, 0.1 * R / 1000.0, n=10)
    lats2, lons2 = circle(10, 20, 0.1 / 1000.0, n=10, figure=UNIT_SPHERE)
    np.testing.assert_allclose(lats1, lats2, rtol=0, atol=1e-12)
    np.testing.assert_allclose(lons1, lons2, rtol=0, atol=1e-12)


@pytest.mark.parametrize("radius", [r for r in RADII_KM if r <= 10000])
def test_circle_haversine_distance(radius):
    for lat, lon in CENTRES:
        lats, lons = circle(lat, lon, radius, n=100)
        dist = haversine_distance((lat, lon), (lats, lons))
        np.testing.assert_allclose(
            dist, radius * 1000.0, rtol=1e-9, atol=1e-6, err_msg=f"centre=({lat}, {lon}) radius={radius}"
        )


def test_circle_haversine_smoke():
    lats, lons = circle(51.5, -0.1, 500)
    dist = haversine_distance((51.5, -0.1), (lats, lons))
    np.testing.assert_allclose(dist, 500000.0, rtol=1e-9, atol=1e-6)


# ---------------------------------------------------------------------------
# order of the points
# ---------------------------------------------------------------------------


def test_circle_cardinal_points():
    lats, lons = circle(0, 0, 1000, n=5)
    d = _delta_deg(1000)
    np.testing.assert_allclose(lats, [d, 0, -d, 0, d], rtol=0, atol=1e-9)
    # the centre longitude is >= 0: the western point is at 360 - d
    np.testing.assert_allclose(lons, [0, d, 0, 360 - d, 0], rtol=0, atol=1e-9)


def test_circle_cardinal_points_negative_lon():
    lats, lons = circle(0, -0.0001, 1000, n=5)
    d = _delta_deg(1000)
    np.testing.assert_allclose(lats, [d, 0, -d, 0, d], rtol=0, atol=1e-9)
    np.testing.assert_allclose(lons, [-0.0001, d - 0.0001, -0.0001, -d - 0.0001, -0.0001], rtol=0, atol=1e-9)


def _bearing(lat1, lon1, lat2, lon2):
    """Initial bearing (rad, from north, clockwise) from point 1 to the points 2, written independently."""
    p1, l1 = np.deg2rad(lat1), np.deg2rad(lon1)
    p2, l2 = np.deg2rad(lat2), np.deg2rad(lon2)
    dl = l2 - l1
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.arctan2(y, x)


@pytest.mark.parametrize("radius", [1.0, 10.0, 500.0, 1000.0, 2500.0, 5000.0, 10000.0])
@pytest.mark.parametrize("n", [5, 37, 100])
def test_circle_bearing(radius, n):
    count = 0
    for lat, lon in CENTRES + [(80.0, 10.0), (-80.0, -10.0), (79.0, 0.0), (0.0, 0.0), (-45.0, 170.0)]:
        if abs(lat) > 80 or _near_pole(lat, radius):
            continue
        count += 1
        lats, lons = circle(lat, lon, radius, n=n)
        theta = np.deg2rad(360.0 * np.arange(n - 1) / (n - 1))
        b = _bearing(lat, lon, lats[:-1], lons[:-1])
        np.testing.assert_allclose(
            _angdiff(b, theta), 0.0, rtol=0, atol=1e-9, err_msg=f"centre=({lat}, {lon}) radius={radius}"
        )
    assert count > 0


def test_circle_bearing_is_clockwise():
    lats, lons = circle(40, 10, 300, n=9)
    # first point north of the centre, then to the east (higher longitude), then south
    assert lats[0] > 40
    assert lons[0] == pytest.approx(10.0, abs=1e-9)
    assert lons[1] > 10 and lats[1] > 40
    assert lons[2] > 10
    assert lats[4] < 40  # 180 degrees is index 4 of 8 distinct points, i.e. due south
    assert lons[4] == pytest.approx(10.0, abs=1e-9)
    assert lons[6] < 10


# ---------------------------------------------------------------------------
# spacing
# ---------------------------------------------------------------------------


def _check_spacing(lat, lon, radius, figure, n=100):
    lats, lons = circle(lat, lon, radius, n=n, figure=figure)
    p = _vec(lats[:-1], lons[:-1])
    steps = _angle(p, np.roll(p, -1, axis=0))
    np.testing.assert_allclose(
        steps, np.full_like(steps, steps.mean()), rtol=1e-9, atol=1e-14, err_msg=f"centre=({lat}, {lon}) r={radius}"
    )


@pytest.mark.parametrize("radius", RADII_KM)
def test_circle_equal_spacing(radius):
    for lat, lon in CENTRES:
        _check_spacing(lat, lon, radius, IFS_SPHERE)


@pytest.mark.parametrize("angle", UNIT_ANGLES)
def test_circle_equal_spacing_unit_sphere(angle):
    for lat, lon in CENTRES:
        _check_spacing(lat, lon, 0.001 * angle, UNIT_SPHERE)


# ---------------------------------------------------------------------------
# poles
# ---------------------------------------------------------------------------


def test_circle_containing_pole():
    lats, lons = circle(80, 0, 2000, n=101)  # the 100 distinct points include the bearing 180
    assert np.all(np.isfinite(lats))
    assert np.all(np.isfinite(lons))
    assert np.all(np.abs(lats) <= 90)
    c = _vec(80.0, 0.0)
    dist = _angle(c, _vec(lats, lons)) * R
    np.testing.assert_allclose(dist, 2000000.0, rtol=1e-9, atol=1e-6)
    # 2000 km is 17.98 degrees: the circle goes beyond the pole, to the other side (lon 180).
    # The first point (due north) is therefore at latitude 180 - (80 + 17.98) = 82.02.
    d = _delta_deg(2000)
    assert d > 10
    assert lats[0] == pytest.approx(90 - (80 + d - 90), abs=1e-9)
    assert np.min(lats) == pytest.approx(80 - d, abs=1e-9)
    # the circle encloses the pole: its points surround it on both sides in longitude
    assert np.any(np.cos(np.deg2rad(lons)) > 0.9)
    assert np.any(np.cos(np.deg2rad(lons)) < -0.9)


@pytest.mark.parametrize("pole", [90.0, -90.0])
@pytest.mark.parametrize("lon", [0.0, 123.0, -45.0])
@pytest.mark.parametrize("radius", [0.001, 1.0, 1000.0, 5000.0, 10000.0])
def test_circle_pole_centre(pole, lon, radius):
    lats, lons = circle(pole, lon, radius)
    d = _delta_deg(radius)
    expected = 90 - d if pole > 0 else -90 + d
    np.testing.assert_allclose(lats, expected, rtol=0, atol=1e-9)
    assert np.all(np.isfinite(lons))
    # the bearings are those of the limit of approaching the pole along the meridian ``lon``: the
    # direction north, which at the north pole continues beyond it, is on the meridian lon + 180
    first = lon + 180.0 if pole > 0 else lon
    assert _angdiff(np.deg2rad(lons[0]), np.deg2rad(first)) < 1e-9


# ---------------------------------------------------------------------------
# longitude range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("lon", [0.0, 1.0, 179.0, 180.0, 200.0, 359.0, 720.5])
@pytest.mark.parametrize("lat", [0.0, 60.0, -85.0, 90.0])
@pytest.mark.parametrize("radius", [10.0, 500.0, 5000.0])
def test_circle_lon_range_non_negative_centre(lat, lon, radius):
    _, lons = circle(lat, lon, radius)
    assert np.all(lons >= 0)
    assert np.all(lons < 360)


@pytest.mark.parametrize("lon", [-1e-9, -0.5, -179.0, -180.0, -181.0, -540.5])
@pytest.mark.parametrize("lat", [0.0, 60.0, -85.0, -90.0])
@pytest.mark.parametrize("radius", [10.0, 500.0, 5000.0])
def test_circle_lon_range_negative_centre(lat, lon, radius):
    _, lons = circle(lat, lon, radius)
    assert np.all(lons >= -180)
    assert np.all(lons < 180)


def test_circle_lon_range_dateline():
    # centre east of the dateline: points cross it and are given as 180..360
    _, lons = circle(0, 179, 500)
    assert np.all((lons >= 0) & (lons < 360))
    assert np.any(lons > 180)
    assert np.any(lons < 180)

    # centre west of the dateline: points cross it and are given as -180..180
    _, lons = circle(0, -179, 500)
    assert np.all((lons >= -180) & (lons < 180))
    assert np.any(lons > 0)
    assert np.any(lons < 0)


def test_circle_lon_range_negative_zero():
    _, lons = circle(0.0, -0.0, 500)
    assert np.all((lons >= 0) & (lons < 360))


# ---------------------------------------------------------------------------
# outward normals
# ---------------------------------------------------------------------------


def _basis(lat, lon):
    phi = np.deg2rad(lat)
    lam = np.deg2rad(lon)
    east = np.stack([-np.sin(lam), np.cos(lam), 0 * lam], axis=-1)
    north = np.stack([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)], axis=-1)
    return east, north


@pytest.mark.parametrize("radius", [10.0, 100.0, 1000.0, 2500.0, 5000.0])
def test_circle_polyline_normal_outwards(radius):
    count = 0
    centres = [c for c in CENTRES if abs(c[0]) <= 70] + [(70.0, 0.0), (-70.0, 90.0), (0.0, 180.0)]
    for lat, lon in centres:
        if _near_pole(lat, radius):
            continue
        count += 1
        lats, lons = circle(lat, lon, radius, n=100)
        east, north = polyline_normal(lats, lons)
        e, u = _basis(lats, lons)
        normal = east[:, None] * e + north[:, None] * u

        c = _vec(lat, lon)
        p = _vec(lats, lons)
        out = (p @ c)[:, None] * p - c
        out /= np.linalg.norm(out, axis=-1, keepdims=True)
        dots = np.sum(normal * out, axis=-1)
        assert np.all(dots > 0.99), f"centre=({lat}, {lon}) radius={radius}: min={dots.min()}"
    assert count > 0


def test_circle_polyline_normal_undefined_at_pole_vertex():
    # a quarter of the circumference around (0, 0) starts at the north pole, where polyline_normal is undefined
    lats, lons = circle(0, 0, np.pi / 2 * R / 1000.0)
    assert lats[0] == 90
    with pytest.raises(ValueError, match="pole"):
        polyline_normal(lats, lons)


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        True,
        False,
        np.True_,
        1j,
        np.complex128(1 + 0j),
        np.array(1j),
        np.array(True),
        "10",
        "a",
        b"1",
        None,
        object(),
        decimal.Decimal("1"),
        [10.0],
        (10.0,),
        np.array([10.0]),
        np.array([[10.0]]),
        np.array(["10"]),
        np.array("10"),
        np.array([]),
    ],
)
@pytest.mark.parametrize("arg", ["lat", "lon", "radius"])
def test_circle_invalid_scalar_type(arg, value):
    kwargs = {"lat": 10.0, "lon": 20.0, "radius": 100.0}
    kwargs[arg] = value
    with pytest.raises(ValueError, match=f"^circle: {arg} must be a real-valued numeric scalar"):
        circle(**kwargs)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, np.float32("nan"), np.array(np.nan)])
@pytest.mark.parametrize("arg", ["lat", "lon", "radius"])
def test_circle_not_finite(arg, value):
    kwargs = {"lat": 10.0, "lon": 20.0, "radius": 100.0}
    kwargs[arg] = value
    with pytest.raises(ValueError, match=f"^circle: {arg} must be finite"):
        circle(**kwargs)


@pytest.mark.parametrize("lat", [90.0001, -90.0001, 91, -100, 1e10])
def test_circle_lat_out_of_range(lat):
    with pytest.raises(ValueError, match="^circle: lat must be in the range"):
        circle(lat, 0, 10)


@pytest.mark.parametrize("radius", [0, 0.0, -0.0, -1, -1e-9, np.float64(-5)])
def test_circle_radius_not_positive(radius):
    with pytest.raises(ValueError, match="^circle: radius must be positive"):
        circle(0, 0, radius)


@pytest.mark.parametrize("radius", [20015.9, 20016, 30000, 1e9])
def test_circle_radius_reaches_antipode(radius):
    with pytest.raises(ValueError, match="^circle: radius=.* reaches the antipode"):
        circle(0, 0, radius)


def test_circle_radius_reaches_antipode_unit_sphere():
    with pytest.raises(ValueError, match="^circle: radius=.* reaches the antipode"):
        circle(0, 0, 0.001 * np.pi + 1e-6, figure=UNIT_SPHERE)
    # but just under it works
    lats, lons = circle(0, 0, 0.001 * (np.pi - 1e-6), figure=UNIT_SPHERE)
    assert np.all(np.isfinite(lats)) and np.all(np.isfinite(lons))


def test_circle_radius_equal_to_antipode():
    # delta == pi exactly, whatever the rounding: use the unit sphere, where delta = radius * 1000
    with pytest.raises(ValueError, match="^circle: "):
        circle(0, 0, np.pi / 1000.0 * (1 + 1e-12), figure=UNIT_SPHERE)


@pytest.mark.parametrize("n", [5.0, 5.5, np.float64(5), "5", None, True, False, np.bool_(True), [5], np.array([5]), 1j])
def test_circle_n_not_int(n):
    with pytest.raises(ValueError, match="^circle: n must be an int"):
        circle(0, 0, 10, n=n)


@pytest.mark.parametrize("n", [3, 2, 1, 0, -1, -100, np.int64(3), np.int32(0)])
def test_circle_n_too_small(n):
    with pytest.raises(ValueError, match="^circle: n must be at least 4"):
        circle(0, 0, 10, n=n)


# ---------------------------------------------------------------------------
# accepted input types
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lat, lon, radius, n",
    [
        (10, 20, 100, 5),
        (np.int64(10), np.int32(20), np.uint8(100), np.int64(5)),
        (np.float32(10), np.float64(20), np.float16(100), np.int8(5)),
        (np.array(10.0), np.array(20), np.array(100.0), 5),
    ],
)
def test_circle_accepted_types(lat, lon, radius, n):
    lats, lons = circle(lat, lon, radius, n=n)
    assert lats.shape == lons.shape == (5,)
    dist = _angle(_vec(float(lat), float(lon)), _vec(lats, lons)) * R
    np.testing.assert_allclose(dist, float(radius) * 1000.0, rtol=1e-9, atol=1e-6)


def test_circle_n_min():
    lats, lons = circle(0, 0, 100, n=4)
    assert lats.shape == (4,)
    assert lats[0] == lats[-1]


def test_circle_does_not_change_inputs_for_numpy_scalars():
    lat, lon, radius = np.float64(10.0), np.float64(20.0), np.float64(100.0)
    a = circle(lat, lon, radius)
    b = circle(10.0, 20.0, 100.0)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
