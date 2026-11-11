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
from earthkit.geo.polyline import _mercator_difference, polyline_normals, sample_polyline

METHODS = ["great_circle", "rhumb", "linear"]

SQRT_HALF = np.sqrt(0.5)


def _xyz(lat, lon):
    return np.stack(latlon_to_xyz(lat, lon), axis=-1)


def _basis(lat, lon):
    """East and north unit vectors (``(..., 3)``), written independently of the implementation."""
    phi = np.deg2rad(np.asarray(lat, dtype=float))
    lam = np.deg2rad(np.asarray(lon, dtype=float))
    east = np.stack([-np.sin(lam), np.cos(lam), 0 * lam], axis=-1)
    north = np.stack([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)], axis=-1)
    return east, north


def _normal_3d(lat, lon, east, north):
    e, u = _basis(lat, lon)
    return east[:, None] * e + north[:, None] * u


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


# ---------------------------------------------------------------------------
# basic structure and cardinal directions
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_structure(method):
    east, north = polyline_normals([0, 5, 10], [0, 10, 20], method=method)
    assert east.shape == (3,)
    assert north.shape == (3,)
    assert east.dtype == float
    assert north.dtype == float

    # array-likes are accepted
    east2, north2 = polyline_normals((0, 5, 10), np.array([0, 10, 20]), method=method)
    np.testing.assert_array_equal(east, east2)
    np.testing.assert_array_equal(north, north2)


def test_polyline_normals_default_method_is_great_circle():
    lats, lons = [10.0, 30.0, 20.0], [0.0, 40.0, 100.0]
    a = polyline_normals(lats, lons)
    b = polyline_normals(lats, lons, method="great_circle")
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons,expected",
    [
        pytest.param([0, 0], [0, 10], (0, 1), id="east"),
        pytest.param([0, 0], [10, 0], (0, -1), id="west"),
        pytest.param([0, 20], [0, 0], (-1, 0), id="north"),
        pytest.param([20, 0], [0, 0], (1, 0), id="south"),
        pytest.param([0, 0, 0], [0, 10, 20], (0, 1), id="east-3"),
        pytest.param([0, 10, 20], [30, 30, 30], (-1, 0), id="north-3"),
        pytest.param([-20, -10, 40], [-170, -170, -170], (-1, 0), id="north-southern-hemisphere"),
    ],
)
def test_polyline_normals_cardinal_directions(method, lats, lons, expected):
    east, north = polyline_normals(lats, lons, method=method)
    np.testing.assert_allclose(east, expected[0], atol=1e-12)
    np.testing.assert_allclose(north, expected[1], atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_unit_length_random(method):
    rng = np.random.default_rng(20261022)
    for _ in range(20):
        k = int(rng.integers(2, 12))
        lats = rng.uniform(-80, 80, k)
        lons = rng.uniform(-180, 180, k)
        east, north = polyline_normals(lats, lons, method=method)
        assert np.all(np.isfinite(east))
        assert np.all(np.isfinite(north))
        np.testing.assert_allclose(east**2 + north**2, 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# great circle
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "a,b",
    [
        ((0.0, 0.0), (0.0, 90.0)),
        ((51.45, -0.97), (44.49, 11.34)),
        ((-30.0, 20.0), (60.0, -100.0)),
        ((70.0, 170.0), (-40.0, -160.0)),
        ((10.0, 0.0), (-10.0, 175.0)),
    ],
)
def test_polyline_normals_great_circle_is_normal_to_the_plane(a, b):
    lat, lon = sample_polyline([a[0], b[0]], [a[1], b[1]], n=50, method="great_circle")
    east, north = polyline_normals(lat, lon, method="great_circle")

    expected = _unit(np.cross(_xyz(*a), _xyz(*b)))
    got = _normal_3d(lat, lon, east, north)
    np.testing.assert_allclose(got, np.broadcast_to(expected, got.shape), atol=1e-9)

    # the same for the two waypoints alone
    east, north = polyline_normals([a[0], b[0]], [a[1], b[1]], method="great_circle")
    got = _normal_3d([a[0], b[0]], [a[1], b[1]], east, north)
    np.testing.assert_allclose(got, np.broadcast_to(expected, got.shape), atol=1e-12)


def _initial_bearing(lat1, lon1, lat2, lon2):
    """Initial bearing (radians) of the great circle from point 1 to point 2."""
    p1, p2, dl = np.deg2rad(lat1), np.deg2rad(lat2), np.deg2rad(lon2 - lon1)
    return np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))


@pytest.mark.parametrize(
    "a,b",
    [
        ((51.45, -0.97), (44.49, 11.34)),
        ((-30.0, 20.0), (60.0, -100.0)),
        ((10.0, 170.0), (35.0, -150.0)),
    ],
)
def test_polyline_normals_great_circle_bearings(a, b):
    east, north = polyline_normals([a[0], b[0]], [a[1], b[1]], method="great_circle")

    # normal = (-cos(bearing), sin(bearing)) with the bearing of the direction of travel
    start = _initial_bearing(*a, *b)
    # the final bearing is the reverse of the initial bearing of the way back
    end = _initial_bearing(*b, *a) + np.pi

    np.testing.assert_allclose([east[0], north[0]], [-np.cos(start), np.sin(start)], atol=1e-12)
    np.testing.assert_allclose([east[1], north[1]], [-np.cos(end), np.sin(end)], atol=1e-12)


# ---------------------------------------------------------------------------
# rhumb and linear, independent values
# ---------------------------------------------------------------------------


def test_polyline_normals_rhumb_bearing():
    # from (0, 0) to (45, dlon): dpsi = asinh(tan(45°)) = asinh(1), so a longitude
    # difference of the same size in radians gives a bearing of exactly 45°
    dlon = np.rad2deg(np.arcsinh(1.0))
    east, north = polyline_normals([0, 45], [0, dlon], method="rhumb")
    np.testing.assert_allclose(east, [-SQRT_HALF, -SQRT_HALF], atol=1e-12)
    np.testing.assert_allclose(north, [SQRT_HALF, SQRT_HALF], atol=1e-12)


def test_polyline_normals_linear_tangent_is_evaluated_at_each_end():
    east, north = polyline_normals([0, 10], [0, 10], method="linear")

    c = np.cos(np.deg2rad(10))
    np.testing.assert_allclose([east[0], north[0]], [-SQRT_HALF, SQRT_HALF], atol=1e-12)
    np.testing.assert_allclose([east[1], north[1]], np.array([-1.0, c]) / np.hypot(1.0, c), atol=1e-12)


@pytest.mark.parametrize("method", ["rhumb", "linear"])
def test_polyline_normals_longitude_difference_as_given(method):
    # a full parallel heading east
    east, north = polyline_normals([0, 0], [0, 360], method=method)
    np.testing.assert_allclose(east, [0, 0], atol=1e-12)
    np.testing.assert_allclose(north, [1, 1], atol=1e-12)

    # the long way west: 340° of longitude, not the short way east
    east, north = polyline_normals([0, 0], [170, -170], method=method)
    np.testing.assert_allclose(east, [0, 0], atol=1e-12)
    np.testing.assert_allclose(north, [-1, -1], atol=1e-12)


@pytest.mark.parametrize(
    "lats,expected_east",
    [
        pytest.param([1e-15, 0.0], 1.0, id="south-1e-15"),
        pytest.param([0.0, 1e-15], -1.0, id="north-1e-15"),
        pytest.param([45.0, np.nextafter(45.0, -np.inf)], 1.0, id="south-one-ulp"),
        pytest.param([45.0, np.nextafter(45.0, np.inf)], -1.0, id="north-one-ulp"),
        pytest.param([-60.0, np.nextafter(-60.0, -np.inf)], 1.0, id="south-one-ulp-southern-hemisphere"),
        pytest.param([80.0, 80.0 + 1e-13], -1.0, id="north-1e-13-high-latitude"),
    ],
)
def test_polyline_normals_rhumb_tiny_latitude_difference(lats, expected_east):
    # the difference of the Mercator latitudes of nearby latitudes must not cancel to zero,
    # which would give a wrong (due north) bearing for a line going south
    east, north = polyline_normals(lats, [0.0, 0.0], method="rhumb")
    np.testing.assert_allclose(east, expected_east, atol=1e-12)
    np.testing.assert_allclose(north, 0.0, atol=1e-12)


def test_polyline_normals_rhumb_tiny_latitude_difference_with_longitude_difference():
    # the bearing is east-ish: a tiny drop in latitude over a normal longitude difference
    east, north = polyline_normals([1e-15, 0.0], [0.0, 10.0], method="rhumb")
    np.testing.assert_allclose(east, [0, 0], atol=1e-12)
    np.testing.assert_allclose(north, [1, 1], atol=1e-12)


def test_polyline_normals_rhumb_u_turn_with_tiny_latitude_differences():
    # one ulp down then one ulp up
    lats = [45.0, np.nextafter(45.0, -np.inf), 45.0]
    with pytest.raises(ValueError, match=r"polyline_normals: the line reverses its direction at vertex 1"):
        polyline_normals(lats, [0.0, 0.0, 0.0], method="rhumb")


def test_polyline_normals_rhumb_matches_the_mercator_formula():
    # away from tiny latitude differences, the bearing is the one of the textbook formula
    rng = np.random.default_rng(1234)
    count = 0
    while count < 200:
        lat = rng.uniform(-85, 85, 2)
        lon = rng.uniform(-180, 180, 2)
        if abs(lat[1] - lat[0]) < 0.5:
            continue
        count += 1

        psi = np.log(np.tan(np.pi / 4 + np.deg2rad(lat) / 2))
        beta = np.arctan2(np.deg2rad(lon[1] - lon[0]), psi[1] - psi[0])

        east, north = polyline_normals(lat, lon, method="rhumb")
        np.testing.assert_allclose(east, -np.cos(beta), atol=1e-12)
        np.testing.assert_allclose(north, np.sin(beta), atol=1e-12)


def test_polyline_normals_rhumb_bearing_is_continuous_for_small_latitude_differences():
    # the normal at a latitude difference going to zero tends to the one of a parallel
    for dlat in [1e-3, 1e-6, 1e-9, 1e-12, 1e-15]:
        east, north = polyline_normals([30.0, 30.0 + dlat], [0.0, 10.0], method="rhumb")
        beta = np.arctan2(np.deg2rad(10.0), np.deg2rad(dlat) / np.cos(np.deg2rad(30.0)))
        np.testing.assert_allclose(east, -np.cos(beta), atol=1e-6)
        np.testing.assert_allclose(north, np.sin(beta), atol=1e-6)


def _ref_psi(lat):
    """High precision isometric latitude of ``lat`` (degrees, exact float input) as a Decimal."""
    with decimal.localcontext() as ctx:
        ctx.prec = 80
        lat = decimal.Decimal(float(lat))
        half_colat = (90 - abs(lat)) * _DECIMAL_PI / 360  # (90° - |lat|) / 2 in radians
        sin, cos = _decimal_sin_cos(half_colat)
        value = -(sin / cos).ln()
        return value if lat > 0 else -value


def _decimal_sin_cos(x):
    """Taylor series of sine and cosine (the context must have a high precision)."""
    sin = term = x
    cos = decimal.Decimal(1)
    cos_term = decimal.Decimal(1)
    tiny = decimal.Decimal(10) ** -75
    for k in range(1, 200):
        cos_term = -cos_term * x * x / ((2 * k - 1) * (2 * k))
        term = -term * x * x / ((2 * k) * (2 * k + 1))
        cos += cos_term
        sin += term
        if abs(term) < tiny and abs(cos_term) < tiny:
            break
    return sin, cos


_DECIMAL_PI = decimal.Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")


def _ref_rhumb_normal(lats, lons):
    """Reference normal ``(east, north)`` of a two-vertex rhumb line at the first vertex, in high precision."""
    with decimal.localcontext() as ctx:
        ctx.prec = 80
        dpsi = _ref_psi(lats[1]) - _ref_psi(lats[0])
        dlam = (decimal.Decimal(float(lons[1])) - decimal.Decimal(float(lons[0]))) * _DECIMAL_PI / 180
    beta = np.arctan2(float(dlam), float(dpsi))
    return -np.cos(beta), np.sin(beta)


def test_polyline_normals_reference_psi_matches_mpmath():
    mpmath = pytest.importorskip("mpmath")
    mpmath.mp.dps = 50
    for lat in [-89.9999999, -45.0, -1e-9, 0.0, 1e-9, 30.0, 89.99, 89.9999999]:
        ref = mpmath.asinh(mpmath.tan(mpmath.radians(mpmath.mpf(lat))))
        assert abs(float(ref) - float(_ref_psi(lat))) <= 1e-13 * max(1.0, abs(float(ref)))


NEAR_POLE_LATS = [
    [0.0, 89.9999999],
    [-89.9999999, 89.9999999],
    [89.99, 89.9999999],
    [89.9999999, 0.0],
    [89.9999999, -89.9999999],
    [89.9999999, 89.99],
    # the southern hemisphere mirrors
    [0.0, -89.9999999],
    [89.9999999, -89.9999999],
    [-89.99, -89.9999999],
    [-89.9999999, 0.0],
    [-89.9999999, 89.9999999],
    [-89.9999999, -89.99],
    # the closest floats to a pole
    [0.0, np.nextafter(90.0, 0.0)],
    [np.nextafter(-90.0, 0.0), 45.0],
    [89.0, 89.99999999999999],
]


@pytest.mark.parametrize("dlon", [10.0, 90.0, -10.0, 0.0])
@pytest.mark.parametrize("lats", NEAR_POLE_LATS)
def test_polyline_normals_rhumb_near_the_poles_against_high_precision(lats, dlon):
    lons = [0.0, dlon]
    east, north = polyline_normals(lats, lons, method="rhumb")
    exp_east, exp_north = _ref_rhumb_normal(lats, lons)

    # a bearing within 1e-9 rad
    np.testing.assert_allclose(east, exp_east, atol=1e-9)
    np.testing.assert_allclose(north, exp_north, atol=1e-9)
    np.testing.assert_allclose(east**2 + north**2, 1.0, atol=1e-12)


def test_polyline_normals_rhumb_pole_reproduction():
    # a large difference of the Mercator latitudes: the bearing is not due north
    east, north = polyline_normals([0, 89.9999999], [0, 90], method="rhumb")
    np.testing.assert_allclose(north, 0.07509114, atol=1e-6)
    np.testing.assert_allclose(east, -np.sqrt(1 - 0.07509114**2), atol=1e-6)
    np.testing.assert_allclose(east**2 + north**2, 1.0, atol=1e-12)

    # towards the south pole and east: the bearing is nearly due south, so the normal (to its
    # left) points nearly east, and slightly north
    east, north = polyline_normals([0, -89.9999999], [0, 90], method="rhumb")
    np.testing.assert_allclose(north, 0.07509114, atol=1e-6)
    np.testing.assert_allclose(east, np.sqrt(1 - 0.07509114**2), atol=1e-6)


@pytest.mark.parametrize("sign", [1.0, -1.0])
@pytest.mark.parametrize("start", [0.0, 30.0, -50.0, 80.0])
@pytest.mark.parametrize("below_and_above", [(-1e-3, 1e-3), (-1e-6, 1e-6), (-1e-10, 1e-10)])
def test_polyline_normals_rhumb_is_continuous_at_the_switch_over(sign, start, below_and_above):
    # The accurate difference of the Mercator latitudes has two formulas, switching at
    # tanh(dpsi) = 0.5. Check the bearing just below and just above that, at both ends.
    threshold = np.arctanh(0.5)

    results = []
    for shift in below_and_above:
        psi_end = float(_ref_psi(start)) + sign * (threshold + shift)
        end = float(np.rad2deg(np.arcsin(np.tanh(psi_end))))
        lats = [start, end]

        # the argument of the atanh of the stable formula, written independently
        dphi = np.deg2rad(lats[1] - lats[0])
        mid = np.deg2rad((lats[0] + lats[1]) / 2)
        x = 2 * np.cos(mid) * np.sin(dphi / 2) / (np.sin(dphi / 2) ** 2 + np.cos(mid) ** 2)
        results.append(abs(x) < 0.5)

        for dlon in (10.0, 90.0):
            east, north = polyline_normals(lats, [0.0, dlon], method="rhumb")
            exp_east, exp_north = _ref_rhumb_normal(lats, [0.0, dlon])
            np.testing.assert_allclose(east, exp_east, atol=1e-9)
            np.testing.assert_allclose(north, exp_north, atol=1e-9)

    # the two sides are on different sides of the switch-over
    assert results == [True, False]


def _assert_rhumb_matches_reference(lats, dlons=(10.0, -10.0, 90.0)):
    for dlon in dlons:
        lons = [0.0, dlon]
        east, north = polyline_normals(lats, lons, method="rhumb")
        exp_east, exp_north = _ref_rhumb_normal(lats, lons)
        np.testing.assert_allclose(east, exp_east, atol=1e-9, err_msg=f"lats={lats!r}, dlon={dlon}")
        np.testing.assert_allclose(north, exp_north, atol=1e-9, err_msg=f"lats={lats!r}, dlon={dlon}")


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_polyline_normals_rhumb_both_ends_next_to_the_same_pole_reproduction(sign):
    lats = [sign * 89.99999999999997, sign * 89.99999999999996]

    # the difference of the Mercator latitudes: ln(tan(c1 / 2) / tan(c2 / 2)) = ln(2 / 3)
    dpsi = float(_ref_psi(lats[1]) - _ref_psi(lats[0]))
    np.testing.assert_allclose(abs(dpsi), 0.40546511, atol=1e-8)
    np.testing.assert_allclose(_mercator_difference(lats[0], lats[1]), dpsi, rtol=1e-9)
    np.testing.assert_allclose(_mercator_difference(lats[1], lats[0]), -dpsi, rtol=1e-9)

    _assert_rhumb_matches_reference(lats)
    _assert_rhumb_matches_reference(lats[::-1])


def _floats_down_from_the_pole(count):
    values = [90.0]
    for _ in range(count):
        values.append(float(np.nextafter(values[-1], 0.0)))
    return values[1:]  # not the pole itself


NEXT_TO_THE_POLE_PAIRS = (
    [
        (89.99999999999997, 89.99999999999996),
        (89.9999999, 89.99999999),
        (89.99999999, 89.9999999999),
        (89.999999999999, 89.9999999999999),
        (89.9999, 89.99999999999999),
        (89.99999999999999, 89.9),
        (89.99999999999997, 89.99),
        (89.99999999999997, 60.0),
        (89.99999999999997, 45.0),
    ]
    # all the pairs of the 8 floats closest to the pole
    + [(a, b) for a in _floats_down_from_the_pole(8) for b in _floats_down_from_the_pole(8) if a != b]
    # adjacent floats further from the pole
    + [(float(np.nextafter(v, 0.0)), v) for v in (89.99, 89.0, 80.0, 60.0, 50.0)]
)


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_polyline_normals_rhumb_next_to_the_poles_both_directions(sign):
    for a, b in NEXT_TO_THE_POLE_PAIRS:
        _assert_rhumb_matches_reference([sign * a, sign * b])
        _assert_rhumb_matches_reference([sign * b, sign * a])


def test_polyline_normals_rhumb_next_to_the_pole_is_not_a_u_turn_or_degenerate():
    # adjacent floats next to the pole: a short but real segment, not a zero-length one
    a = 89.99999999999999
    b = float(np.nextafter(a, 0.0))
    for sign in (1.0, -1.0):
        east, north = polyline_normals([sign * a, sign * b], [0.0, 0.0], method="rhumb")
        # due south in the north, due north in the south: the normal is east, resp. west
        np.testing.assert_allclose(east, sign, atol=1e-12)
        np.testing.assert_allclose(north, 0.0, atol=1e-12)


def _around(value):
    return [
        value,
        float(np.nextafter(value, 90.0)),
        float(np.nextafter(value, -90.0)),
        value + 1e-9,
        value - 1e-9,
        value + 1e-3,
        value - 1e-3,
    ]


def test_polyline_normals_rhumb_is_continuous_across_the_45_degree_switch():
    # The colatitude formula is used if both ends are at least 45° from the equator, in the same
    # hemisphere. Segments with an end at 45 ± tiny, just inside and just outside of that region.
    others = [0.0, 30.0, 44.9999999, 45.0, 45.0000001, 60.0, 89.9999999, 89.99999999999999, -10.0, -45.0, -60.0]
    for sign in (1.0, -1.0):
        for end in _around(45.0):
            for other in others:
                _assert_rhumb_matches_reference([sign * end, sign * other])
                _assert_rhumb_matches_reference([sign * other, sign * end])

    # both ends close to 45, on both sides, across the equator and across the switch
    for a in _around(45.0):
        for b in _around(45.0):
            if a != b:
                _assert_rhumb_matches_reference([a, b])
                _assert_rhumb_matches_reference([-a, -b])
                _assert_rhumb_matches_reference([a, -b])


# ---------------------------------------------------------------------------
# consistency with the sampled path
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons",
    [
        ([10.0, 40.0, 30.0, 60.0], [10.0, 30.0, 60.0, 80.0]),
        ([-40.0, -20.0, 10.0], [100.0, 140.0, 150.0]),
        ([20.0, 50.0], [5.0, 120.0]),
    ],
)
def test_polyline_normals_perpendicular_to_finite_difference_tangent(method, lats, lons):
    _check_perpendicular_to_finite_difference(method, lats, lons)


@pytest.mark.parametrize(
    "method,lats,lons",
    [
        # sampled longitudes jump from 180 to -180 (the interior points are in [-180, 180))
        ("great_circle", [10.0, 30.0], [150.0, -150.0]),
        ("great_circle", [10.0, 30.0, 20.0], [160.0, -170.0, -120.0]),
        ("rhumb", [0.0, 20.0, 30.0], [-10.0, 170.0, 190.0]),
        ("rhumb", [-20.0, 10.0], [-30.0, 200.0]),
        # no wrapping for the linear method, but the longitudes are beyond 180
        ("linear", [0.0, 20.0, 30.0], [-10.0, 170.0, 190.0]),
        ("linear", [10.0, 30.0], [170.0, 190.0]),
    ],
)
def test_polyline_normals_perpendicular_across_the_longitude_wrap(method, lats, lons):
    raw_lon = _check_perpendicular_to_finite_difference(method, lats, lons)
    if method != "linear":
        # the sampled longitudes do cross the wrap boundary
        assert np.max(np.abs(np.diff(raw_lon))) > 180


def _check_perpendicular_to_finite_difference(method, lats, lons):
    """Check the normals of a densely sampled line. Return the sampled longitudes as ``sample_polyline`` gives them."""
    n = 4000
    lat, lon = sample_polyline(lats, lons, n=n, method=method)
    raw_lon = lon
    if method == "rhumb":
        # the longitude difference is used as given: undo the wrapping of the sampled longitudes
        lon = np.unwrap(lon, period=360)
    east, north = polyline_normals(lat, lon, method=method)
    normal = _normal_3d(lat, lon, east, north)

    # central difference in 3-D, projected onto the local east/north basis
    xyz = _xyz(lat, lon)
    diff = xyz[2:] - xyz[:-2]
    e, u = _basis(lat[1:-1], lon[1:-1])
    t_east = np.sum(diff * e, axis=-1)
    t_north = np.sum(diff * u, axis=-1)
    norm = np.hypot(t_east, t_north)
    t_east, t_north = t_east / norm, t_north / norm

    # only segment interiors: leave out the waypoints and their neighbours
    waypoints = [i for i in range(n) if any(lat[i] == la and raw_lon[i] == lo for la, lo in zip(lats, lons))]
    assert len(waypoints) >= len(lats)
    keep = np.ones(n, dtype=bool)
    keep[0] = keep[-1] = False
    for i in waypoints:
        keep[max(i - 2, 0) : i + 3] = False
    inner = keep[1:-1]
    assert np.sum(inner) > n // 2

    dot = east[1:-1] * t_east + north[1:-1] * t_north
    np.testing.assert_allclose(dot[inner], 0.0, atol=1e-4)

    # and the normal is on the left: normal = (-t_north, t_east)
    np.testing.assert_allclose(east[1:-1][inner], -t_north[inner], atol=1e-4)
    np.testing.assert_allclose(north[1:-1][inner], t_east[inner], atol=1e-4)

    # the 3-D normal is tangent to the sphere
    np.testing.assert_allclose(np.sum(normal * xyz, axis=-1), 0.0, atol=1e-12)

    return raw_lon


# ---------------------------------------------------------------------------
# bisector, first and last vertices
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_bisector_of_a_right_angle_turn(method):
    # east along the equator then north along the meridian 0: a left turn at (0, 0)
    east, north = polyline_normals([0, 0, 10], [-10, 0, 0], method=method)
    np.testing.assert_allclose(east, [0, -SQRT_HALF, -1], atol=1e-12)
    np.testing.assert_allclose(north, [1, SQRT_HALF, 0], atol=1e-12)

    # north then east: a right turn, the normal is the bisector of (0, 1) and (-1, 0) pointing
    # to the left of (1, 1)/sqrt(2)
    east, north = polyline_normals([-10, 0, 0], [0, 0, 10], method=method)
    np.testing.assert_allclose(east, [-1, -SQRT_HALF, 0], atol=1e-12)
    np.testing.assert_allclose(north, [0, SQRT_HALF, 1], atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_bisector_is_the_normalised_sum_of_the_tangents(method):
    lats = [10.0, 35.0, 20.0]
    lons = [0.0, 25.0, 70.0]
    east, north = polyline_normals(lats, lons, method=method)

    e_in, n_in = polyline_normals(lats[:2], lons[:2], method=method)
    e_out, n_out = polyline_normals(lats[1:], lons[1:], method=method)

    # normals of the two segments at the middle vertex (the normal is the tangent rotated by 90°)
    expected = _unit([e_in[1] + e_out[0], n_in[1] + n_out[0]])
    np.testing.assert_allclose([east[1], north[1]], expected, atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_first_and_last_use_their_single_segment(method):
    lats = [10.0, 35.0, 20.0, -5.0]
    lons = [0.0, 25.0, 70.0, 60.0]
    east, north = polyline_normals(lats, lons, method=method)

    e_first, n_first = polyline_normals(lats[:2], lons[:2], method=method)
    e_last, n_last = polyline_normals(lats[-2:], lons[-2:], method=method)

    np.testing.assert_allclose([east[0], north[0]], [e_first[0], n_first[0]], atol=1e-12)
    np.testing.assert_allclose([east[-1], north[-1]], [e_last[1], n_last[1]], atol=1e-12)


# ---------------------------------------------------------------------------
# zero-length segments
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_leading_and_trailing_duplicates(method):
    east, north = polyline_normals([0, 0, 0, 0], [0, 0, 10, 10], method=method)
    np.testing.assert_allclose(east, 0.0, atol=1e-12)
    np.testing.assert_allclose(north, 1.0, atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_duplicates_get_the_same_normal(method):
    # east, then north, with the corner repeated, and the last vertex repeated
    east, north = polyline_normals([0, 0, 0, 10, 10], [-10, 0, 0, 0, 0], method=method)
    np.testing.assert_allclose(east, [0, -SQRT_HALF, -SQRT_HALF, -1, -1], atol=1e-12)
    np.testing.assert_allclose(north, [1, SQRT_HALF, SQRT_HALF, 0, 0], atol=1e-12)

    # same as without the duplicates
    e, n = polyline_normals([0, 0, 10], [-10, 0, 0], method=method)
    np.testing.assert_allclose(east[[0, 1, 3]], e, atol=1e-12)
    np.testing.assert_allclose(north[[0, 1, 3]], n, atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_duplicate_first_vertex_with_turn(method):
    east, north = polyline_normals([0, 0, 0, 10], [-10, -10, 0, 0], method=method)
    np.testing.assert_allclose(east, [0, 0, -SQRT_HALF, -1], atol=1e-12)
    np.testing.assert_allclose(north, [1, 1, SQRT_HALF, 0], atol=1e-12)


def test_polyline_normals_great_circle_longitude_aliases_are_duplicates():
    east, north = polyline_normals([0, 0, 0], [0, 360, 10], method="great_circle")
    np.testing.assert_allclose(east, 0.0, atol=1e-12)
    np.testing.assert_allclose(north, 1.0, atol=1e-12)

    east, north = polyline_normals([0, 0, 0], [-170, 190, 200], method="great_circle")
    np.testing.assert_allclose(east, 0.0, atol=1e-12)
    np.testing.assert_allclose(north, 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_error_not_1d(method):
    with pytest.raises(ValueError, match=r"polyline_normals: lats and lons must be 1D"):
        polyline_normals([[0, 1], [2, 3]], [[0, 1], [2, 3]], method=method)

    with pytest.raises(ValueError, match=r"polyline_normals: lats and lons must be 1D"):
        polyline_normals(0.0, 0.0, method=method)


def test_polyline_normals_error_length():
    with pytest.raises(ValueError, match=r"polyline_normals: lats and lons must have the same length! 3 != 2"):
        polyline_normals([0, 1, 2], [0, 1])


@pytest.mark.parametrize("lats,lons", [([], []), ([0.0], [0.0])])
def test_polyline_normals_error_fewer_than_two(lats, lons):
    with pytest.raises(ValueError, match=r"polyline_normals: at least 2 vertices are required"):
        polyline_normals(lats, lons)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_polyline_normals_error_not_finite(bad):
    with pytest.raises(ValueError, match=r"polyline_normals: lats and lons must be finite"):
        polyline_normals([0, bad], [0, 10])

    with pytest.raises(ValueError, match=r"polyline_normals: lats and lons must be finite"):
        polyline_normals([0, 0], [0, bad])


@pytest.mark.parametrize("lat", [90.0001, -90.0001, 100.0, -200.0])
def test_polyline_normals_error_latitude_range(lat):
    with pytest.raises(ValueError, match=r"polyline_normals: latitudes must be in the range"):
        polyline_normals([0, lat], [0, 10])


def test_polyline_normals_error_unknown_method():
    with pytest.raises(ValueError, match=r"polyline_normals: unknown method='bad'"):
        polyline_normals([0, 1], [0, 1], method="bad")


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons,index",
    [
        ([0, 90], [0, 0], 1),
        ([-90, 0], [0, 0], 0),
        ([0, 90, 0], [0, 0, 20], 1),
        ([0, 0, -90], [0, 10, 10], 2),
    ],
)
def test_polyline_normals_error_pole_vertex(method, lats, lons, index):
    with pytest.raises(ValueError, match=rf"polyline_normals: vertex {index} is at a pole"):
        polyline_normals(lats, lons, method=method)


@pytest.mark.parametrize(
    "lats,lons,index",
    [
        ([0, 0], [0, 180], 0),
        ([10, -10], [20, 200], 0),
        ([0, 0, 30, -30], [0, 10, 20, 200], 2),
    ],
)
def test_polyline_normals_error_antipodal_great_circle(lats, lons, index):
    with pytest.raises(ValueError, match=rf"polyline_normals: segment {index} .*antipodal"):
        polyline_normals(lats, lons, method="great_circle")


@pytest.mark.parametrize("method", ["rhumb", "linear"])
def test_polyline_normals_antipodal_is_fine_for_other_methods(method):
    east, north = polyline_normals([0, 0], [0, 180], method=method)
    np.testing.assert_allclose(east, 0.0, atol=1e-12)
    np.testing.assert_allclose(north, 1.0, atol=1e-12)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons",
    [
        ([0, 0, 0], [0, 10, 0]),
        ([0, 10, 0], [0, 0, 0]),
        ([0, 0, 0, 0], [0, 10, 10, 0]),
    ],
)
def test_polyline_normals_error_u_turn(method, lats, lons):
    with pytest.raises(ValueError, match=r"polyline_normals: the line reverses its direction at vertex (1|2)"):
        polyline_normals(lats, lons, method=method)


def test_polyline_normals_error_u_turn_names_the_vertex():
    with pytest.raises(ValueError, match=r"at vertex 2 "):
        polyline_normals([0, 0, 0, 0], [0, 5, 10, 5], method="great_circle")


@pytest.mark.parametrize("method", METHODS)
def test_polyline_normals_error_all_segments_degenerate(method):
    with pytest.raises(ValueError, match=r"polyline_normals: all the segments have zero length"):
        polyline_normals([10, 10, 10], [20, 20, 20], method=method)


def test_polyline_normals_error_all_segments_degenerate_great_circle_alias():
    with pytest.raises(ValueError, match=r"polyline_normals: all the segments have zero length"):
        polyline_normals([0, 0], [0, 360], method="great_circle")
