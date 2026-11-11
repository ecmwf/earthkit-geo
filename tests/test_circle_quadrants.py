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
import re

import numpy as np
import pytest

from earthkit.geo.coord import latlon_to_xyz
from earthkit.geo.figure import IFS_SPHERE, UNIT_SPHERE
from earthkit.geo.polyline import circle, circle_quadrants

R = IFS_SPHERE.radius  # m

# sets of 4 radii (m) in the order NE, SE, SW, NW: all equal, all different (both orders),
# mixed sizes and mixed with 1 m short of the antipode
RADIUS_SETS = [
    (1_000_000.0, 1_000_000.0, 1_000_000.0, 1_000_000.0),
    (1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0),
    (4_000_000.0, 3_000_000.0, 2_000_000.0, 1_000_000.0),
    (1.0, 1_000.0, 1_000_000.0, 10_000_000.0),
    (np.pi * R - 1.0, 1.0, 15_000_000.0, 1_000.0),
    (10_000_000.0, 1.0, 10_000_000.0, 1.0),
]

# radii (m) with every pair of adjacent ones (including NW/NE) different: the boundary
# points are at different places, so all the n - 1 points are distinct
DISTINCT_ADJACENT = [r for r in RADIUS_SETS[1:]]

# angular radii (rad) sets for the unit sphere (radius in m = angle, as the radius of the sphere is 1 m)
UNIT_ANGLE_SETS = [
    (1e-6, 0.1, 1.0, 3.0),
    (np.pi - 1e-6, 1e-6, 2.0, np.pi - 0.5),
]


def _random_centres():
    rng = np.random.default_rng(20261104)
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

# values that are not real-valued numeric scalars
INVALID_SCALARS = [
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
]

NOT_FINITE = [np.nan, np.inf, -np.inf, np.float32("nan"), np.array(np.nan)]


def _vec(lat, lon):
    return np.stack(latlon_to_xyz(lat, lon), axis=-1)


def _angle(a, b):
    """Independent angular distance (rad) between the unit vectors ``a`` and ``b`` (last axis of size 3)."""
    return np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, axis=-1))


def _angdiff(a, b):
    """Absolute difference of the angles (rad) modulo a full turn."""
    return np.abs((np.asarray(a) - np.asarray(b) + np.pi) % (2 * np.pi) - np.pi)


def _delta_deg(radius_m, figure=IFS_SPHERE):
    return float(np.rad2deg(radius_m / figure.radius))


def _basis(lat, lon):
    phi = np.deg2rad(lat)
    lam = np.deg2rad(lon)
    east = np.stack([-np.sin(lam), np.cos(lam), 0 * lam], axis=-1)
    north = np.stack([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)], axis=-1)
    return east, north


def _bearing(lat1, lon1, lat2, lon2):
    """Initial bearing (rad, from north, clockwise) from point 1 to the points 2, written independently."""
    p1, l1 = np.deg2rad(lat1), np.deg2rad(lon1)
    p2, l2 = np.deg2rad(lat2), np.deg2rad(lon2)
    dl = l2 - l1
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.arctan2(y, x)


def _bearing_vec(lat0, lon0, lat, lon):
    """Bearings (rad, clockwise from north) from a centre to points, computed with vectors.

    More robust than :func:`_bearing` next to the poles; still noisy next to the antipode.
    """
    c = _vec(lat0, lon0)
    east, north = _basis(lat0, lon0)
    p = _vec(lat, lon)
    rel = p - (p @ c)[:, None] * c
    return np.arctan2(rel @ east, rel @ north)


def _split(n):
    """Expected number of the ``n - 1`` distinct points in each quadrant (NE, SE, SW, NW)."""
    q, extra = divmod(n - 1, 4)
    return [q + 1 if i < extra else q for i in range(4)]


def _near_pole(lat, radii, limit=1.0):
    """True if any quadrant arc passes within ``limit`` degrees of a pole."""
    return any(abs(_delta_deg(r) - d) <= limit for r in radii for d in (90 - lat, 90 + lat))


# ---------------------------------------------------------------------------
# structure, shape, closure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [9, 10, 12, 13, 100, 361])
def test_quadrants_length_and_closed(n):
    lats, lons = circle_quadrants(51.5, -0.1, (500_000.0, 600_000.0, 700_000.0, 800_000.0), n=n)
    for a in (lats, lons):
        assert isinstance(a, np.ndarray)
        assert a.shape == (n,)
        assert a.dtype == np.float64
    # the last point is an exact copy of the first one
    assert lats[-1] == lats[0]
    assert lons[-1] == lons[0]


def test_quadrants_default_n():
    lats, lons = circle_quadrants(10, 20, (100_000.0, 200_000.0, 300_000.0, 400_000.0))
    assert lats.shape == (100,)
    assert lons.shape == (100,)


def test_quadrants_n_min():
    lats, lons = circle_quadrants(0, 0, (100_000.0, 200_000.0, 300_000.0, 400_000.0), n=9)
    assert lats.shape == (9,)
    assert lats[0] == lats[-1]


@pytest.mark.parametrize("radii", DISTINCT_ADJACENT)
def test_quadrants_distinct_points(radii):
    lats, lons = circle_quadrants(51.5, -0.1, radii, n=100)
    pts = set(zip(lats[:-1].tolist(), lons[:-1].tolist()))
    assert len(pts) == 99


@pytest.mark.parametrize("n", [9, 10, 11, 12, 13, 14, 15, 16, 17, 100, 101, 102, 103, 361])
def test_quadrants_split(n):
    # every remainder of (n - 1) mod 4; the earlier quadrants get the extra points
    expected = _split(n)
    assert all(c >= 2 for c in expected)
    radii = np.array([1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0])
    lats, lons = circle_quadrants(45.0, 10.0, radii, n=n)
    dist = _angle(_vec(45.0, 10.0), _vec(lats[:-1], lons[:-1])) * R
    # the radii are far apart, the quadrant of a point is the one with the nearest radius
    group = np.abs(dist[:, None] - radii[None, :]).argmin(axis=1)
    counts = np.bincount(group, minlength=4)
    np.testing.assert_array_equal(counts, expected)
    # the groups are in order: the first c0 points are NE, then SE, SW and NW
    np.testing.assert_array_equal(group, np.repeat(np.arange(4), np.asarray(expected)))


# ---------------------------------------------------------------------------
# distance from the centre, per quadrant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("radii", RADIUS_SETS)
def test_quadrants_distance_to_centre(radii):
    deltas = np.repeat(np.asarray(radii) / R, _split(100))
    for lat, lon in CENTRES:
        lats, lons = circle_quadrants(lat, lon, radii, n=100)
        c = _vec(lat, lon)
        p = _vec(lats, lons)
        dist = _angle(c, p) * R
        # the closing copy is at the place of the first point, i.e. at the NE radius
        np.testing.assert_allclose(
            dist, np.append(deltas, deltas[0]) * R, rtol=1e-9, atol=1e-6, err_msg=f"centre=({lat}, {lon})"
        )


@pytest.mark.parametrize("angles", UNIT_ANGLE_SETS)
def test_quadrants_distance_to_centre_unit_sphere(angles):
    radii = list(angles)
    deltas = np.repeat(np.asarray(radii) / UNIT_SPHERE.radius, _split(100))
    for lat, lon in CENTRES:
        lats, lons = circle_quadrants(lat, lon, radii, n=100, figure=UNIT_SPHERE)
        dist = _angle(_vec(lat, lon), _vec(lats, lons)) * UNIT_SPHERE.radius
        np.testing.assert_allclose(
            dist, np.append(deltas, deltas[0]), rtol=1e-9, atol=1e-6, err_msg=f"centre=({lat}, {lon})"
        )


def test_quadrants_figure_gives_the_same_outline_for_the_same_angle():
    # an angular radius pattern of (0.1, 0.2, 0.3, 0.4) rad on both spheres
    a = [0.1, 0.2, 0.3, 0.4]
    lats1, lons1 = circle_quadrants(10, 20, [x * R for x in a], n=10)
    lats2, lons2 = circle_quadrants(10, 20, a, n=10, figure=UNIT_SPHERE)
    np.testing.assert_allclose(lats1, lats2, rtol=0, atol=1e-12)
    np.testing.assert_allclose(lons1, lons2, rtol=0, atol=1e-12)


# ---------------------------------------------------------------------------
# bearings: quadrants, spacing, start due north, clockwise
# ---------------------------------------------------------------------------


def test_quadrants_cardinal_example():
    # the n=9 example of the docstring: the points are at the cardinal bearings, on the
    # equator or on the meridian of the centre, at the radius of their quadrant
    radii = [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0]
    d = [_delta_deg(r) for r in radii]
    lats, lons = circle_quadrants(0, 0, radii, n=9)
    np.testing.assert_allclose(lats, [d[0], 0, 0, -d[1], -d[2], 0, 0, d[3], d[0]], rtol=0, atol=1e-9)
    np.testing.assert_allclose(lons, [0, d[0], d[1], 0, 0, 360 - d[2], 360 - d[3], 0, 0], rtol=0, atol=1e-9)


@pytest.mark.parametrize("n", [9, 13, 100])
def test_quadrants_expected_bearings(n):
    # each quadrant runs from its start bearing to its end bearing inclusive, equally spaced
    radii = [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0]
    theta = np.deg2rad(np.concatenate([np.linspace(90.0 * i, 90.0 * (i + 1), c) for i, c in enumerate(_split(n))]))
    assert theta[0] == 0.0  # starts due north
    assert theta[-1] == np.deg2rad(360.0)  # the NW quadrant ends due north
    lats, lons = circle_quadrants(45.0, 10.0, radii, n=n)
    b = _bearing(45.0, 10.0, lats[:-1], lons[:-1])
    np.testing.assert_allclose(_angdiff(b, theta), 0.0, rtol=0, atol=1e-9)


@pytest.mark.parametrize(
    "radii",
    [(1_000.0, 10_000.0, 500_000.0, 1_000_000.0), (1_000_000.0, 2_500_000.0, 3_000_000.0, 5_000_000.0)],
)
@pytest.mark.parametrize("n", [9, 37, 100])
def test_quadrants_bearings(radii, n):
    count = 0
    theta = np.deg2rad(np.concatenate([np.linspace(90.0 * i, 90.0 * (i + 1), c) for i, c in enumerate(_split(n))]))
    for lat, lon in CENTRES + [(80.0, 10.0), (-80.0, -10.0), (79.0, 0.0), (0.0, 0.0), (-45.0, 170.0)]:
        if abs(lat) > 80 or _near_pole(lat, radii):
            continue
        count += 1
        lats, lons = circle_quadrants(lat, lon, radii, n=n)
        b = _bearing(lat, lon, lats[:-1], lons[:-1])
        np.testing.assert_allclose(
            _angdiff(b, theta), 0.0, rtol=0, atol=1e-9, err_msg=f"centre=({lat}, {lon}) radii={radii}"
        )
    assert count > 0


@pytest.mark.parametrize(
    "radii",
    [(10_000.0, 2_000_000.0, 5_000_000.0, 9_000_000.0), (9_000_000.0, 5_000_000.0, 2_000_000.0, 10_000.0)],
)
def test_quadrants_bearings_stay_in_quadrant(radii):
    # the bearing of a point is inside the bearing range of its own quadrant
    count = 0
    radii = np.asarray(radii)
    for lat, lon in CENTRES:
        if abs(lat) > 80 or _near_pole(lat, radii):
            continue
        count += 1
        lats, lons = circle_quadrants(lat, lon, radii, n=100)
        dist = _angle(_vec(lat, lon), _vec(lats[:-1], lons[:-1])) * R
        group = np.abs(dist[:, None] - radii[None, :]).argmin(axis=1)
        b = _bearing_vec(lat, lon, lats[:-1], lons[:-1])
        # distance of the bearing from the quadrant interval [90g, 90g + 90] on the circle:
        # measured from the middle of the quadrant, it must not exceed a half quadrant
        half = np.deg2rad(45.0) + 1e-7
        np.testing.assert_array_less(
            _angdiff(b, np.deg2rad(45.0 + 90.0 * group)), half, err_msg=f"centre=({lat}, {lon})"
        )
    assert count > 0


def test_quadrants_first_point_due_north_at_ne_radius():
    lat, lon = 40.0, 10.0
    radii = [300_000.0, 500_000.0, 700_000.0, 900_000.0]
    lats, lons = circle_quadrants(lat, lon, radii, n=100)
    dist = _angle(_vec(lat, lon), _vec(lats, lons)) * R
    assert lons[0] == pytest.approx(lon, abs=1e-9)
    assert lats[0] > lat
    assert dist[0] == pytest.approx(radii[0], rel=1e-9, abs=1e-6)
    # the last distinct point (index n - 2, before the closing copy) is also due north, at the NW radius
    assert lons[-2] == pytest.approx(lon, abs=1e-9)
    assert lats[-2] > lat
    assert dist[-2] == pytest.approx(radii[3], rel=1e-9, abs=1e-6)


def test_quadrants_order_is_clockwise():
    lat, lon = 40.0, 10.0
    n = 13
    lats, lons = circle_quadrants(lat, lon, [300_000.0] * 4, n=n)
    # the quadrants have 3 points each, at the bearings 0, 45, 90 | 90, 135, 180 | ...:
    # the points go north, east, south, west and back to north, clockwise seen from above
    assert lats[0] > lat and lons[0] == pytest.approx(lon, abs=1e-9)  # 0 degrees, due north
    assert lats[1] > lat and lons[1] > lon  # 45 degrees
    assert lons[2] > lon  # 90 degrees, due east
    assert lats[5] < lat and lons[5] == pytest.approx(lon, abs=1e-9)  # 180 degrees, due south
    assert lats[7] < lat and lons[7] < lon  # 225 degrees
    assert lats[10] > lat and lons[10] < lon  # 315 degrees
    assert lats[-2] > lat and lons[-2] == pytest.approx(lon, abs=1e-9)  # 360 degrees, due north again
    b = _bearing(lat, lon, lats[:-1], lons[:-1]) % (2 * np.pi)
    b[-1] += 2 * np.pi  # the last point is due north: count it as 360 degrees, not 0
    diffs = np.diff(b)
    assert np.all(diffs > -1e-9)  # the bearings never turn back
    assert np.all(diffs < np.deg2rad(90.0) + 1e-9)  # and never jump by more than a quadrant


# ---------------------------------------------------------------------------
# radial step at the quadrant boundaries
# ---------------------------------------------------------------------------


def _boundary_pairs(n):
    """Indices of the pairs of points at the boundaries 90, 180, 270 and 0 (360) degrees."""
    c = _split(n)
    return [
        (c[0] - 1, c[0]),
        (c[0] + c[1] - 1, c[0] + c[1]),
        (sum(c[:3]) - 1, sum(c[:3])),
        (n - 2, 0),
    ]


def test_quadrants_radial_step_at_boundaries():
    lat, lon = 30.0, 20.0
    radii = [1_000_000.0, 2_500_000.0, 1_500_000.0, 4_000_000.0]
    n = 100
    lats, lons = circle_quadrants(lat, lon, radii, n=n)
    p = _vec(lats, lons)
    dist = _angle(_vec(lat, lon), p) * R
    b = _bearing(lat, lon, lats, lons)
    for (i, j), boundary in zip(_boundary_pairs(n), [90, 180, 270, 360]):
        ra, rb = radii[(boundary // 90 - 1) % 4], radii[(boundary // 90) % 4]
        # both points are exactly on the boundary bearing (no tangential step)
        assert _angdiff(b[i], np.deg2rad(boundary)) < 1e-9, boundary
        assert _angdiff(b[j], np.deg2rad(boundary)) < 1e-9, boundary
        # each is at the radius of its own quadrant
        np.testing.assert_allclose(dist[i], ra, rtol=1e-9, atol=1e-6)
        np.testing.assert_allclose(dist[j], rb, rtol=1e-9, atol=1e-6)
        # the distance between them is the radial step
        np.testing.assert_allclose(_angle(p[i], p[j]) * R, abs(ra - rb), rtol=1e-9, atol=1e-6)


def test_quadrants_no_step_between_equal_adjacent_radii():
    radii = [2_000_000.0, 2_000_000.0, 3_000_000.0, 3_000_000.0]
    n = 100
    lats, lons = circle_quadrants(45.0, 10.0, radii, n=n)
    (i_ne, i_se), (i_se2, i_sw), (i_sw2, i_nw), (i_nw2, i_ne2) = _boundary_pairs(n)
    # NE/SE and SW/NW have the same radius: the boundary points are exact copies of each other
    assert lats[i_ne] == lats[i_se] and lons[i_ne] == lons[i_se]
    assert lats[i_sw2] == lats[i_nw] and lons[i_sw2] == lons[i_nw]
    # SE/SW and NW/NE differ: there is a step
    assert not (lats[i_se2] == lats[i_sw] and lons[i_se2] == lons[i_sw])
    assert not (lats[i_nw2] == lats[i_ne2] and lons[i_nw2] == lons[i_ne2])


def test_quadrants_equal_radii_gives_circle_distances():
    for radius in [1.0, 1_000.0, 1_000_000.0, 10_000_000.0, np.pi * R - 1.0]:
        for lat, lon in CENTRES:
            lats_q, lons_q = circle_quadrants(lat, lon, [radius] * 4, n=100)
            lats_c, lons_c = circle(lat, lon, radius, n=100)
            c = _vec(lat, lon)
            dq = _angle(c, _vec(lats_q, lons_q)) * R
            dc = _angle(c, _vec(lats_c, lons_c)) * R
            err_msg = f"centre=({lat}, {lon}) radius={radius}"
            # all the points are at the radius, exactly like the points of the circle
            np.testing.assert_allclose(dq, dc, rtol=1e-9, atol=1e-6, err_msg=err_msg)
            np.testing.assert_allclose(dq, radius, rtol=1e-9, atol=1e-6, err_msg=err_msg)
            assert np.ptp(dq) <= 2 * (1e-9 * radius + 1e-6), err_msg


# ---------------------------------------------------------------------------
# longitude range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "radii",
    [(10_000.0, 500_000.0, 5_000_000.0, 1_000_000.0), (4_000_000.0, 100_000.0, 8_000_000.0, 50_000.0)],
)
@pytest.mark.parametrize("lon", [0.0, 1.0, 179.0, 180.0, 200.0, 359.0, 720.5])
@pytest.mark.parametrize("lat", [0.0, 60.0, -85.0, 90.0])
def test_quadrants_lon_range_non_negative_centre(lat, lon, radii):
    _, lons = circle_quadrants(lat, lon, radii)
    assert np.all(lons >= 0)
    assert np.all(lons < 360)


@pytest.mark.parametrize(
    "radii",
    [(10_000.0, 500_000.0, 5_000_000.0, 1_000_000.0), (4_000_000.0, 100_000.0, 8_000_000.0, 50_000.0)],
)
@pytest.mark.parametrize("lon", [-1e-9, -0.5, -179.0, -180.0, -181.0, -540.5])
@pytest.mark.parametrize("lat", [0.0, 60.0, -85.0, -90.0])
def test_quadrants_lon_range_negative_centre(lat, lon, radii):
    _, lons = circle_quadrants(lat, lon, radii)
    assert np.all(lons >= -180)
    assert np.all(lons < 180)


def test_quadrants_lon_range_dateline():
    # centre east of the dateline: points cross it and are given as 180..360
    _, lons = circle_quadrants(0, 179, [500_000, 500_000, 2_000_000, 2_000_000])
    assert np.all((lons >= 0) & (lons < 360))
    assert np.any(lons > 180)
    assert np.any(lons < 180)

    # centre west of the dateline: points cross it and are given as -180..180
    _, lons = circle_quadrants(0, -179, [500_000, 500_000, 2_000_000, 2_000_000])
    assert np.all((lons >= -180) & (lons < 180))
    assert np.any(lons > 0)
    assert np.any(lons < 0)


def test_quadrants_lon_range_negative_zero():
    _, lons = circle_quadrants(0.0, -0.0, [500_000, 500_000, 500_000, 500_000])
    assert np.all((lons >= 0) & (lons < 360))


# ---------------------------------------------------------------------------
# poles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pole", [90.0, -90.0])
@pytest.mark.parametrize("lon", [0.0, 123.0, -45.0])
def test_quadrants_pole_centre(pole, lon):
    radii = [1_000_000.0, 2_000_000.0, 1_500_000.0, 3_000_000.0]
    n = 13
    lats, lons = circle_quadrants(pole, lon, radii, n=n)
    d = [_delta_deg(r) for r in radii]
    expected = np.repeat([90 - x if pole > 0 else -90 + x for x in d], _split(n))
    np.testing.assert_allclose(lats[:-1], expected, rtol=0, atol=1e-9)
    assert lats[-1] == lats[0]
    assert np.all(np.isfinite(lons))
    # the bearings are those of the limit of approaching the pole along the meridian ``lon``: the
    # direction north, which at the north pole continues beyond it, is on the meridian lon + 180
    first = lon + 180.0 if pole > 0 else lon
    assert _angdiff(np.deg2rad(lons[0]), np.deg2rad(first)) < 1e-9
    minimum, maximum = (0.0, 360.0) if lon >= 0 else (-180.0, 180.0)
    assert np.all((lons >= minimum) & (lons < maximum))


def test_quadrants_containing_pole():
    radii = [1_500_000.0, 2_500_000.0, 1_000_000.0, 3_000_000.0]
    lats, lons = circle_quadrants(80, 0, radii, n=101)
    assert np.all(np.isfinite(lats))
    assert np.all(np.isfinite(lons))
    assert np.all(np.abs(lats) <= 90)
    deltas = np.repeat(np.asarray(radii) / R, _split(101))
    dist = _angle(_vec(80.0, 0.0), _vec(lats, lons)) * R
    np.testing.assert_allclose(dist, np.append(deltas, deltas[0]) * R, rtol=1e-9, atol=1e-6)
    # the NE quadrant is at 1500 km = 13.48 degrees: the first point (due north) goes beyond
    # the pole, to the other side (lon 180), at latitude 90 - (80 + d - 90)
    d = _delta_deg(radii[0])
    assert d > 10
    assert lats[0] == pytest.approx(90 - (80 + d - 90), abs=1e-9)
    assert lons[0] == pytest.approx(180.0, abs=1e-9)
    assert np.any(lats > 80)  # points on the other side of the pole
    # the southernmost point is due south, at the end of the SE arc (the 3000 km NW arc, which
    # runs between bearings 270 and 360, only reaches latitude 61.4 degrees at its west end)
    assert np.min(lats) == pytest.approx(80 - _delta_deg(radii[1]), abs=1e-9)
    # the outline passes on both sides of the pole in longitude
    assert np.any(np.cos(np.deg2rad(lons)) > 0.9)
    assert np.any(np.cos(np.deg2rad(lons)) < -0.9)


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------


def test_quadrants_radius_reaches_antipode_unit_sphere():
    # on the unit sphere delta = radius, so only radii below pi are below the antipode
    with pytest.raises(ValueError, match=r"^circle_quadrants: radii\[3\]=.* m reaches the antipode"):
        circle_quadrants(0, 0, [1e-6, 1e-6, 1e-6, np.pi + 1e-6], figure=UNIT_SPHERE)
    # but just under it works
    lats, lons = circle_quadrants(0, 0, [np.pi - 1e-6] * 4, figure=UNIT_SPHERE)
    assert np.all(np.isfinite(lats)) and np.all(np.isfinite(lons))


def test_quadrants_radius_above_antipode():
    # delta just above pi: use the unit sphere, where delta = radius
    with pytest.raises(ValueError, match=r"^circle_quadrants: "):
        circle_quadrants(0, 0, [np.pi * (1 + 1e-12)] * 4, figure=UNIT_SPHERE)


def test_quadrants_radius_equal_to_antipode():
    # the exact limit radius, pi * figure.radius metres, gives delta == pi exactly, which is
    # rejected too (delta >= pi raises)
    for figure in (IFS_SPHERE, UNIT_SPHERE):
        limit = np.pi * figure.radius
        with pytest.raises(ValueError, match=r"^circle_quadrants: radii\[0\]=.* m reaches the antipode"):
            circle_quadrants(0, 0, [limit, limit, limit, limit], figure=figure)


# ---------------------------------------------------------------------------
# validation: n
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [5.0, 5.5, np.float64(5), "5", None, True, False, np.bool_(True), [5], np.array([5]), 1j])
def test_quadrants_n_not_int(n):
    with pytest.raises(ValueError, match=r"^circle_quadrants: n must be an int"):
        circle_quadrants(0, 0, [100_000, 100_000, 100_000, 100_000], n=n)


@pytest.mark.parametrize("n", [8, 3, 2, 1, 0, -1, -100, np.int64(8), np.int32(0)])
def test_quadrants_n_too_small(n):
    with pytest.raises(ValueError, match=r"^circle_quadrants: n must be at least 9"):
        circle_quadrants(0, 0, [100_000, 100_000, 100_000, 100_000], n=n)


# ---------------------------------------------------------------------------
# validation: lat, lon
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("value", INVALID_SCALARS)
@pytest.mark.parametrize("arg", ["lat", "lon"])
def test_quadrants_invalid_scalar_type(arg, value):
    kwargs = {"lat": 10.0, "lon": 20.0, "radii": [100_000.0, 100_000.0, 100_000.0, 100_000.0]}
    kwargs[arg] = value
    with pytest.raises(ValueError, match=f"^circle_quadrants: {arg} must be a real-valued numeric scalar"):
        circle_quadrants(**kwargs)


@pytest.mark.parametrize("value", NOT_FINITE)
@pytest.mark.parametrize("arg", ["lat", "lon"])
def test_quadrants_not_finite(arg, value):
    kwargs = {"lat": 10.0, "lon": 20.0, "radii": [100_000.0, 100_000.0, 100_000.0, 100_000.0]}
    kwargs[arg] = value
    with pytest.raises(ValueError, match=f"^circle_quadrants: {arg} must be finite"):
        circle_quadrants(**kwargs)


@pytest.mark.parametrize("lat", [90.0001, -90.0001, 91, -100, 1e10])
def test_quadrants_lat_out_of_range(lat):
    with pytest.raises(ValueError, match=r"^circle_quadrants: lat must be in the range"):
        circle_quadrants(lat, 0, [10_000, 10_000, 10_000, 10_000])


# ---------------------------------------------------------------------------
# validation: radii
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "radii",
    [
        [],
        [1_000_000.0],
        [1_000_000.0, 2_000_000.0, 3_000_000.0],
        (1_000_000.0, 2_000_000.0, 3_000_000.0),
        [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0, 5_000_000.0],
        np.array([1_000_000.0, 2_000_000.0, 3_000_000.0]),
        np.ones(5),
        np.array([]),
    ],
)
def test_quadrants_radii_wrong_length(radii):
    with pytest.raises(ValueError, match=r"^circle_quadrants: radii must have exactly 4 elements"):
        circle_quadrants(0, 0, radii)


@pytest.mark.parametrize(
    "radii",
    [
        1000,
        1_000_000.0,
        np.float64(1_000_000.0),
        None,
        "1000000",
        b"1234",
        {"a": 1, "b": 2, "c": 3, "d": 4},
        {1_000_000.0},
        iter([1_000_000.0] * 4),
        range(4),
        np.array(1_000_000.0),
        np.full((2, 2), 1_000_000.0),
        np.array([[[1_000_000.0, 2_000_000.0], [3_000_000.0, 4_000_000.0]]]),
    ],
)
def test_quadrants_radii_not_a_sequence(radii):
    with pytest.raises(ValueError, match=r"^circle_quadrants: radii must be a sequence of 4"):
        circle_quadrants(0, 0, radii)


@pytest.mark.parametrize(
    "radii",
    [
        [True, 1_000_000.0, 1_000_000.0, 1_000_000.0],
        [1_000_000.0, np.bool_(False), 1_000_000.0, 1_000_000.0],
        [1_000_000.0, 1_000_000.0, 1_000_000.0, True],
        np.array([True, True, True, True]),
        ["1000000", 1_000_000.0, 1_000_000.0, 1_000_000.0],
        [1_000_000.0, b"1", 1_000_000.0, 1_000_000.0],
        np.array(["10", "20", "30", "40"]),
        [[1_000_000.0], 1_000_000.0, 1_000_000.0, 1_000_000.0],
        ((1_000_000.0, 2_000_000.0), 1_000_000.0, 1_000_000.0, 1_000_000.0),
        [np.array([1_000_000.0]), 1_000_000.0, 1_000_000.0, 1_000_000.0],
        [1_000_000.0 + 0j, 1_000_000.0, 1_000_000.0, 1_000_000.0],
        [1_000_000.0, decimal.Decimal("1000000"), 1_000_000.0, 1_000_000.0],
        [1_000_000.0, None, 1_000_000.0, 1_000_000.0],
        [1_000_000.0, object(), 1_000_000.0, 1_000_000.0],
    ],
)
def test_quadrants_radii_invalid_elements(radii):
    with pytest.raises(ValueError, match=r"^circle_quadrants: radii\[\d\] must be a real-valued numeric scalar"):
        circle_quadrants(0, 0, radii)


@pytest.mark.parametrize("value", NOT_FINITE)
@pytest.mark.parametrize("index", [0, 1, 2, 3])
def test_quadrants_radii_not_finite(index, value):
    radii = [1_000_000.0, 1_000_000.0, 1_000_000.0, 1_000_000.0]
    radii[index] = value
    with pytest.raises(ValueError, match=rf"^circle_quadrants: radii\[{index}\] must be finite"):
        circle_quadrants(0, 0, radii)


@pytest.mark.parametrize("value", [0, 0.0, -0.0, -1, -1e-9, np.float64(-5)])
@pytest.mark.parametrize("index", [0, 1, 2, 3])
def test_quadrants_radii_not_positive(index, value):
    radii = [1_000_000.0, 1_000_000.0, 1_000_000.0, 1_000_000.0]
    radii[index] = value
    with pytest.raises(ValueError, match=rf"^circle_quadrants: radii\[{index}\] must be positive"):
        circle_quadrants(0, 0, radii)


@pytest.mark.parametrize("value", [20_015_900, 20_016_000, 30_000_000, 1e12])
@pytest.mark.parametrize("index", [0, 1, 2, 3])
def test_quadrants_radii_reaches_antipode(index, value):
    radii = [1_000_000.0, 1_000_000.0, 1_000_000.0, 1_000_000.0]
    radii[index] = value
    # the message quotes the radius as a float and gives the exact limit, in metres
    limit = re.escape(f"{np.pi * R}")
    got = re.escape(f"{float(value)}")
    with pytest.raises(
        ValueError, match=rf"^circle_quadrants: radii\[{index}\]={got} m reaches the antipode.*at most {limit} m"
    ):
        circle_quadrants(0, 0, radii)


def test_quadrants_just_below_antipode_works():
    r = np.pi * R - 1.0
    lats, lons = circle_quadrants(0, 0, [r, 1.0, r / 2, 1_000.0], n=9)
    assert np.all(np.isfinite(lats)) and np.all(np.isfinite(lons))


# ---------------------------------------------------------------------------
# accepted input types
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lat, lon, radii, n",
    [
        (10, 20, [1_000_000, 2_000_000, 3_000_000, 4_000_000], 9),
        (np.int64(10), np.int32(20), (1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0), np.int64(9)),
        (
            np.float32(10),
            np.float64(20),
            np.array([1_000_000, 2_000_000, 3_000_000, 4_000_000], dtype=np.float32),
            np.int8(9),
        ),
        (np.array(10.0), np.array(20), np.array([1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0]), 9),
    ],
)
def test_quadrants_accepted_types(lat, lon, radii, n):
    lats, lons = circle_quadrants(lat, lon, radii, n=n)
    assert lats.shape == lons.shape == (n,)
    deltas = np.repeat(np.asarray(radii, dtype=float) / R, _split(n))
    dist = _angle(_vec(float(lat), float(lon)), _vec(lats[:-1], lons[:-1])) * R
    np.testing.assert_allclose(dist, deltas * R, rtol=1e-9, atol=1e-6)


def test_quadrants_does_not_change_inputs():
    radii = [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0]
    arr = np.array([1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0])
    circle_quadrants(10.0, 20.0, radii, n=20)
    circle_quadrants(10.0, 20.0, arr, n=20)
    assert radii == [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0]
    np.testing.assert_array_equal(arr, [1_000_000.0, 2_000_000.0, 3_000_000.0, 4_000_000.0])


# ---------------------------------------------------------------------------
# the error messages of circle() are unchanged (the functions share their helpers)
# ---------------------------------------------------------------------------


def test_circle_error_messages_unchanged():
    with pytest.raises(ValueError, match=r"^circle: lat must be a real-valued numeric scalar"):
        circle(True, 0, 10_000)
    with pytest.raises(ValueError, match=r"^circle: lon must be finite"):
        circle(0, np.nan, 10_000)
    with pytest.raises(ValueError, match=r"^circle: radius must be a real-valued numeric scalar"):
        circle(0, 0, "10")
    with pytest.raises(ValueError, match=r"^circle: lat must be in the range"):
        circle(91, 0, 10_000)
    with pytest.raises(ValueError, match=r"^circle: radius must be positive"):
        circle(0, 0, 0)
    with pytest.raises(ValueError, match=r"^circle: radius=.* m reaches the antipode"):
        circle(0, 0, 30_000_000)
    with pytest.raises(ValueError, match=r"^circle: n must be an int"):
        circle(0, 0, 10_000, n=5.0)
    with pytest.raises(ValueError, match=r"^circle: n must be at least 4"):
        circle(0, 0, 10_000, n=3)
