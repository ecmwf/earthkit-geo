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
from earthkit.geo.distance import haversine_distance
from earthkit.geo.figure import UNIT_SPHERE
from earthkit.geo.polyline import _distribute, _wrap_lon, sample_polyline

METHODS = ["great_circle", "rhumb", "linear"]


def _mercator_y(lat):
    return np.log(np.tan(np.pi / 4 + np.deg2rad(lat) / 2))


def _counts_between(lons, wp_lons):
    """Number of points strictly between consecutive waypoints (lons must be increasing)."""
    return [int(np.sum((lons > a) & (lons < b))) for a, b in zip(wp_lons[:-1], wp_lons[1:])]


def _xyz(lat, lon):
    return np.stack(latlon_to_xyz(lat, lon), axis=-1)


# ---------------------------------------------------------------------------
# basic structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_length_and_waypoint_indices(method):
    # all segments have the same length (20) for every method, so the extra
    # points are shared equally: 3 per segment
    lats = [-30.0, -10.0, 10.0, 30.0]
    lons = [5.0, 5.0, 5.0, 5.0]
    lat, lon = sample_polyline(lats, lons, n=13, method=method)

    assert lat.shape == (13,)
    assert lon.shape == (13,)
    assert lat.dtype == float
    assert lon.dtype == float

    for i, idx in enumerate([0, 4, 8, 12]):
        assert lat[idx] == lats[i]
        assert lon[idx] == lons[i]

    # strictly increasing latitude along the meridian
    assert np.all(np.diff(lat) > 0)
    np.testing.assert_allclose(lon, 5.0)
    np.testing.assert_allclose(np.diff(lat), 5.0)


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_waypoints_exact_unequal_segments(method):
    lats = [10.123456789, 20.456, -5.5, 33.3, 33.3, 41.7]
    lons = [170.1, 20.2, 33.3, 120.4, 120.4, 151.9]
    n = 60
    lat, lon = sample_polyline(lats, lons, n=n, method=method)

    assert lat.shape == (n,)
    assert lon.shape == (n,)
    assert np.all(np.isfinite(lat))
    assert np.all(np.isfinite(lon))

    idx = []
    for la, lo in zip(lats, lons):
        found = np.flatnonzero((lat == la) & (lon == lo))
        assert found.size >= 1
        idx.append(int(found[0]))

    assert idx[0] == 0
    assert idx[-1] == n - 1
    # waypoints 3 and 4 are the same place: consecutive points, no point in between
    assert np.flatnonzero((lat == 33.3) & (lon == 120.4)).tolist() == [idx[3], idx[3] + 1]
    assert idx[3] == idx[4]
    assert idx[0] < idx[1] < idx[2] < idx[3] < idx[5]


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_inputs_types_and_not_modified(method):
    lats = np.array([0, 10, 20])
    lons = np.array([0, 10, 20])
    lats_copy = lats.copy()
    lons_copy = lons.copy()
    lat, lon = sample_polyline(lats, lons, n=7, method=method)
    np.testing.assert_array_equal(lats, lats_copy)
    np.testing.assert_array_equal(lons, lons_copy)
    assert lat.ndim == 1 and lon.ndim == 1
    assert lat.dtype == float and lon.dtype == float

    # lists and tuples are accepted, numpy integers as n too
    lat2, lon2 = sample_polyline((0, 10, 20), (0, 10, 20), n=np.int64(7), method=method)
    np.testing.assert_array_equal(lat, lat2)
    np.testing.assert_array_equal(lon, lon2)


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_default_n_and_method(method):
    lat, lon = sample_polyline([0, 10], [0, 10], method=method)
    assert lat.shape == (100,)
    assert lon.shape == (100,)

    lat_d, lon_d = sample_polyline([0, 10], [0, 10])
    lat_g, lon_g = sample_polyline([0, 10], [0, 10], n=100, method="great_circle")
    np.testing.assert_array_equal(lat_d, lat_g)
    np.testing.assert_array_equal(lon_d, lon_g)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons",
    [
        ([0.0, 10.0], [0.0, 20.0]),
        ([10.5, 20.5, -30.25], [0.0, 400.0, -500.0]),
        ([0.0, 0.0, 0.0], [-180.0, 180.0, 540.0]),
    ],
)
def test_sample_polyline_n_equals_k(method, lats, lons):
    lat, lon = sample_polyline(lats, lons, n=len(lats), method=method)
    np.testing.assert_array_equal(lat, lats)
    np.testing.assert_array_equal(lon, lons)


# ---------------------------------------------------------------------------
# great circle
# ---------------------------------------------------------------------------


def test_sample_polyline_great_circle_equator():
    lat, lon = sample_polyline([0, 0], [0, 90], n=10)
    np.testing.assert_allclose(lon, np.arange(0, 100, 10), atol=1e-9)
    np.testing.assert_allclose(lat, 0, atol=1e-12)
    assert lat[0] == 0 and lon[0] == 0
    assert lat[-1] == 0 and lon[-1] == 90


def test_sample_polyline_great_circle_meridian():
    lat, lon = sample_polyline([-20, 70], [30, 30], n=10)
    np.testing.assert_allclose(lat, np.linspace(-20, 70, 10), atol=1e-9)
    np.testing.assert_allclose(lon, 30, atol=1e-9)


@pytest.mark.parametrize(
    "p1,p2",
    [
        ((51.45, -0.97), (44.49, 11.34)),
        ((10.0, 20.0), (-40.0, 150.0)),
        ((60.0, -120.0), (-10.0, 70.0)),
        ((0.0, 0.0), (89.0, 175.0)),
    ],
)
def test_sample_polyline_great_circle_equal_spacing_and_plane(p1, p2):
    n = 25
    lat, lon = sample_polyline([p1[0], p2[0]], [p1[1], p2[1]], n=n)

    # equally spaced along the great circle
    d = np.array([
        haversine_distance((lat[i], lon[i]), (lat[i + 1], lon[i + 1]), figure=UNIT_SPHERE) for i in range(n - 1)
    ])
    np.testing.assert_allclose(d, d[0], rtol=1e-9)

    # total distance is the great circle distance of the waypoints
    total = haversine_distance(p1, p2, figure=UNIT_SPHERE)
    np.testing.assert_allclose(d.sum(), total, rtol=1e-9)

    # all points are in the plane of the segment
    a = _xyz(*p1)
    b = _xyz(*p2)
    normal = np.cross(a, b)
    normal /= np.linalg.norm(normal)
    pts = _xyz(lat, lon)
    np.testing.assert_allclose(pts @ normal, 0, atol=1e-12)


def test_sample_polyline_great_circle_across_dateline():
    lat, lon = sample_polyline([0, 0], [170, -170], n=21)
    assert lon[0] == 170
    assert lon[-1] == -170
    np.testing.assert_allclose(lat, 0, atol=1e-12)
    # the short way, through +-180
    assert np.all(np.abs(lon) >= 170 - 1e-9)
    interior = lon[1:-1]
    assert np.all(interior >= -180)
    assert np.all(interior < 180)
    # equal 1 degree steps: 170, 171, ..., 179, -180, -179, ..., -170
    expected = np.concatenate([np.arange(171, 180), np.arange(-180, -170)]).astype(float)
    np.testing.assert_allclose(interior, expected, atol=1e-9)


def test_sample_polyline_great_circle_near_antipodal_meridian():
    # ends almost antipodal, but above the rejection threshold: the spacing must stay equal
    lat, lon = sample_polyline([-90, 89.9999999999], [0, 0], n=5)
    assert lat[0] == -90 and lat[-1] == 89.9999999999
    np.testing.assert_allclose(lat, np.linspace(-90, 89.9999999999, 5), atol=1e-6, rtol=0)
    np.testing.assert_allclose(lon, 0, atol=1e-9)

    lat, lon = sample_polyline([-90, 89.9999999999], [0, 0], n=181)
    np.testing.assert_allclose(lat, np.linspace(-90, 89.9999999999, 181), atol=1e-6, rtol=0)


@pytest.mark.parametrize("dlon", [179.999999999, 179.9999999999, 180.0000000001])
def test_sample_polyline_great_circle_near_antipodal_equator(dlon):
    n = 9
    lat, lon = sample_polyline([0, 0], [0, dlon], n=n)
    np.testing.assert_allclose(lat, 0, atol=1e-9)
    # shortest path: eastwards for dlon < 180, westwards (through 0) otherwise
    expected = np.linspace(0, dlon, n) if dlon < 180 else np.linspace(360, dlon, n) % 360
    expected[0], expected[-1] = 0, dlon
    np.testing.assert_allclose(lon, expected, atol=1e-6, rtol=0)


@pytest.mark.parametrize("lat2,lon2", [(-10.0, 180.0 - 1e-10), (20.0, -179.9999999999), (89.99999999, 180.0)])
def test_sample_polyline_great_circle_near_antipodal_equal_spacing(lat2, lon2):
    # p1 is 1e-9 degrees away from the antipode of (lat2, 180)
    p1 = (-lat2 + 1e-9, 0.0)
    n = 21
    lat, lon = sample_polyline([p1[0], lat2], [p1[1], lon2], n=n)
    d = np.array([
        haversine_distance((lat[i], lon[i]), (lat[i + 1], lon[i + 1]), figure=UNIT_SPHERE) for i in range(n - 1)
    ])
    np.testing.assert_allclose(d, d[0], rtol=1e-6)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))


def test_sample_polyline_great_circle_over_the_pole():
    lat, lon = sample_polyline([80, 80], [0, 180], n=11)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))
    # arcsin is ill-conditioned at the pole
    assert np.max(lat) == pytest.approx(90, abs=1e-5)
    # and the end at the pole itself
    lat, lon = sample_polyline([0, 90], [10, 0], n=10)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))
    np.testing.assert_allclose(lat, np.linspace(0, 90, 10), atol=1e-9)


@pytest.mark.parametrize("dlon", [1e-9, 1e-10, 1e-11])
def test_sample_polyline_great_circle_near_coincident(dlon):
    # central angle is around 1.7e-11 or smaller: close to the zero-length threshold
    lat, lon = sample_polyline([0, 0], [0, dlon], n=10)
    assert np.all(np.isfinite(lat))
    assert np.all(np.isfinite(lon))
    assert lon[0] == 0 and lon[-1] == dlon
    assert np.all(np.diff(lon) >= 0)
    assert np.all(lon >= 0) and np.all(lon <= dlon * (1 + 1e-9))
    np.testing.assert_allclose(lat, 0, atol=1e-12)


def test_sample_polyline_great_circle_coincident_is_zero_length():
    # the first segment is a zero-length one: gets no points
    lat, lon = sample_polyline([0, 0, 0], [0, 1e-13, 10], n=5)
    assert lon[0] == 0
    assert lon[1] == 1e-13
    np.testing.assert_allclose(lon[2:], [10 / 3, 20 / 3, 10], atol=1e-9)


def test_sample_polyline_great_circle_global_three_or_more_waypoints():
    lat, lon = sample_polyline([0, 0, 0, 0], [-180, -60, 60, 180], n=361)
    assert lat.shape == (361,)
    np.testing.assert_allclose(lat, 0, atol=1e-12)
    np.testing.assert_allclose(lon, np.arange(-180, 181), atol=1e-9)
    assert lon[0] == -180 and lon[-1] == 180
    assert np.all(lon[1:-1] >= -180) and np.all(lon[1:-1] < 180)


def test_sample_polyline_great_circle_global_zero_to_360():
    lat, lon = sample_polyline([0, 0, 0, 0], [0, 120, 240, 360], n=361)
    np.testing.assert_allclose(lat, 0, atol=1e-12)
    np.testing.assert_allclose(lon, np.arange(0, 361), atol=1e-9)
    assert lon[0] == 0 and lon[-1] == 360
    interior = lon[1:-1]
    assert np.all(interior >= 0) and np.all(interior < 360)


# ---------------------------------------------------------------------------
# rhumb
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "p1,p2",
    [
        ((10.0, -20.0), (50.0, 40.0)),
        ((-30.0, 100.0), (60.0, 120.0)),
        ((0.0, 0.0), (60.0, 60.0)),
        ((70.0, 10.0), (-20.0, -50.0)),
    ],
)
def test_sample_polyline_rhumb_straight_in_mercator(p1, p2):
    lat, lon = sample_polyline([p1[0], p2[0]], [p1[1], p2[1]], n=30, method="rhumb")
    x = np.deg2rad(lon)
    y = _mercator_y(lat)
    dx = x - x[0]
    dy = y - y[0]
    cross = dx * dy[-1] - dy * dx[-1]
    scale = max(abs(dx[-1]), abs(dy[-1]))
    np.testing.assert_allclose(cross / scale**2, 0, atol=1e-12)

    # constant bearing: dlon / dpsi is the same for every step
    ratio = np.diff(x) / np.diff(y)
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-8)


def test_sample_polyline_rhumb_equal_spacing():
    n = 13
    lat, lon = sample_polyline([0, 60], [0, 60], n=n, method="rhumb")
    # latitude is linear in the fraction
    np.testing.assert_allclose(np.diff(lat), 5.0, atol=1e-12)
    np.testing.assert_allclose(lat, np.linspace(0, 60, n), atol=1e-12)

    # rhumb distance (unit sphere) is the same for every step
    phi = np.deg2rad(lat)
    dphi = np.diff(phi)
    dpsi = np.diff(_mercator_y(lat))
    dlam = np.diff(np.deg2rad(lon))
    d = np.sqrt(dphi**2 + (dphi / dpsi * dlam) ** 2)
    np.testing.assert_allclose(d, d[0], rtol=1e-8)

    # endpoints exactly
    assert (lat[0], lon[0]) == (0, 0)
    assert (lat[-1], lon[-1]) == (60, 60)
    # longitude increases with latitude following the Mercator projection, not linearly
    assert np.all(np.diff(lon) > 0)
    assert not np.allclose(lon, np.linspace(0, 60, n))


def test_sample_polyline_rhumb_constant_latitude():
    lat, lon = sample_polyline([45, 45], [10, 50], n=5, method="rhumb")
    np.testing.assert_array_equal(lat, 45)
    np.testing.assert_allclose(lon, [10, 20, 30, 40, 50])


def test_sample_polyline_rhumb_meridian():
    lat, lon = sample_polyline([-10, 50], [20, 20], n=7, method="rhumb")
    np.testing.assert_allclose(lat, np.linspace(-10, 50, 7), atol=1e-12)
    np.testing.assert_array_equal(lon, 20)


def test_sample_polyline_rhumb_global_equator():
    lat, lon = sample_polyline([0, 0], [-180, 180], n=361, method="rhumb")
    np.testing.assert_array_equal(lat, 0)
    np.testing.assert_allclose(lon, np.arange(-180, 181), atol=1e-9)
    assert lon[0] == -180 and lon[-1] == 180
    assert np.all(lon[1:-1] >= -180) and np.all(lon[1:-1] < 180)


def test_sample_polyline_rhumb_global_parallel():
    lat, lon = sample_polyline([45, 45], [-180, 180], n=361, method="rhumb")
    np.testing.assert_array_equal(lat, 45)
    np.testing.assert_allclose(lon, np.arange(-180, 181), atol=1e-9)
    # all longitudes covered, one degree apart
    np.testing.assert_allclose(np.diff(lon), 1, atol=1e-9)


def test_sample_polyline_rhumb_goes_west_through_zero():
    # the longitude difference is taken as given: -340 degrees
    lat, lon = sample_polyline([0, 0], [170, -170], n=341, method="rhumb")
    np.testing.assert_array_equal(lat, 0)
    np.testing.assert_allclose(lon, np.arange(170, -171, -1), atol=1e-9)
    assert 0 in np.round(lon)


def test_sample_polyline_rhumb_crosses_dateline_eastwards():
    lat, lon = sample_polyline([0, 0], [170, 190], n=5, method="rhumb")
    assert lon[0] == 170 and lon[-1] == 190
    # interior is normalised to [0, 360) because all inputs are >= 0
    np.testing.assert_allclose(lon[1:-1], [175, 180, 185])


def test_sample_polyline_rhumb_global_zero_to_360():
    lat, lon = sample_polyline([0, 0], [0, 360], n=361, method="rhumb")
    np.testing.assert_array_equal(lat, 0)
    np.testing.assert_allclose(lon, np.arange(0, 361), atol=1e-9)
    assert lon[-1] == 360
    interior = lon[1:-1]
    assert np.all(interior >= 0) and np.all(interior < 360)


# ---------------------------------------------------------------------------
# linear
# ---------------------------------------------------------------------------


def test_sample_polyline_linear_exact():
    lat, lon = sample_polyline([0, 10], [0, 20], n=11, method="linear")
    np.testing.assert_allclose(lat, np.arange(11), rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(lon, 2 * np.arange(11), rtol=1e-14, atol=1e-14)
    assert (lat[0], lon[0]) == (0, 0)
    assert (lat[-1], lon[-1]) == (10, 20)


def test_sample_polyline_linear_global_zero_to_360():
    lat, lon = sample_polyline([0, 0], [0, 360], n=361, method="linear")
    np.testing.assert_array_equal(lat, 0)
    np.testing.assert_allclose(lon, np.arange(0, 361), atol=1e-9)
    assert lon[-1] == 360
    interior = lon[1:-1]
    assert np.all(interior >= 0) and np.all(interior < 360)


def test_sample_polyline_linear_global_minus180_to_180():
    lat, lon = sample_polyline([20, 20], [-180, 180], n=361, method="linear")
    np.testing.assert_array_equal(lat, 20)
    np.testing.assert_allclose(lon, np.arange(-180, 181), atol=1e-9)


def test_sample_polyline_linear_no_wrapping():
    # goes the long way round and is not normalised
    lat, lon = sample_polyline([0, 0], [170, -170], n=5, method="linear")
    np.testing.assert_allclose(lon, [170, 85, 0, -85, -170])

    lat, lon = sample_polyline([0, 0], [0, 720], n=5, method="linear")
    np.testing.assert_allclose(lon, [0, 180, 360, 540, 720])

    lat, lon = sample_polyline([0, 0], [-10, 400], n=3, method="linear")
    np.testing.assert_allclose(lon, [-10, 195, 400])


# ---------------------------------------------------------------------------
# sharing of the points
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_share_proportional_to_length_1_to_3(method):
    # along the equator the length ratio is 1:3 for all the methods; 8 extra points
    lat, lon = sample_polyline([0, 0, 0], [0, 10, 40], n=11, method=method)
    assert _counts_between(lon, [0, 10, 40]) == [2, 6]
    assert lon[3] == 10
    expected = [0, 10 / 3, 20 / 3, 10] + [10 + 30 * j / 7 for j in range(1, 7)] + [40]
    np.testing.assert_allclose(lon, expected, atol=1e-9)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lon_wp,n,expected",
    [
        # 3 equal segments, 2 extra: quotas 2/3 each -> the earlier ones win
        ([0, 1, 2, 3], 6, [1, 1, 0]),
        # tie between two equal segments: the earlier one wins
        ([0, 1, 2], 4, [1, 0]),
        ([0, 1, 2, 3, 4], 7, [1, 1, 0, 0]),
        # largest remainder in the later segment: quotas 0.4, 1.6, 2.0
        ([0, 1, 5, 10], 8, [0, 2, 2]),
        # quotas 3.5, 1, 0.5: tie between the first and the third
        ([0, 7, 9, 10], 9, [4, 1, 0]),
        # quotas 0.9, 0.9, 1.2: the two biggest remainders win
        ([0, 3, 6, 10], 7, [1, 1, 1]),
        # exact quotas
        ([0, 1, 3, 6], 10, [1, 2, 3]),
    ],
)
def test_sample_polyline_largest_remainder(method, lon_wp, n, expected):
    lats = [0] * len(lon_wp)
    lat, lon = sample_polyline(lats, lon_wp, n=n, method=method)
    assert lat.shape == (n,)
    assert _counts_between(lon, lon_wp) == expected
    assert sum(expected) == n - len(lon_wp)


@pytest.mark.parametrize(
    "lengths,extra,expected",
    [
        ([1, 1, 1], 2, [1, 1, 0]),
        ([1, 1], 1, [1, 0]),
        ([1, 4, 5], 4, [0, 2, 2]),
        ([7, 2, 1], 5, [4, 1, 0]),
        ([0, 1, 1], 1, [0, 1, 0]),
        ([1, 0, 1], 3, [2, 0, 1]),
        ([0, 0, 0], 4, [2, 1, 1]),
        ([0, 0, 0], 3, [1, 1, 1]),
        ([0, 0], 1, [1, 0]),
        ([1, 2, 3], 0, [0, 0, 0]),
        ([0, 0], 0, [0, 0]),
        # floating point noise must not break ties/exact values
        ([0.1 + 0.2, 0.3], 1, [1, 0]),
        ([0.1 * 3, 0.1 * 7], 10, [3, 7]),
        ([2.9999999999999996, 7.0], 10, [3, 7]),
        ([3.0000000000000004, 7.0], 10, [3, 7]),
        # real differences, even small ones, are not ties: the longer segment wins
        ([1, 1.000000001], 1, [0, 1]),
        ([1, 1 + 1e-10], 1, [0, 1]),
        ([1 + 1e-10, 1], 1, [1, 0]),
        ([1, 1.000000001, 1], 1, [0, 1, 0]),
        ([1, 1.000000002, 1.000000001], 2, [0, 1, 1]),
        ([1, 1.000000001], 101, [50, 51]),
    ],
)
def test_polyline_distribute(lengths, extra, expected):
    counts = _distribute(np.array(lengths, dtype=float), extra)
    assert counts.tolist() == expected
    assert counts.sum() == extra


def test_polyline_distribute_always_sums_to_extra():
    rng = np.random.default_rng(1)
    for _ in range(200):
        nseg = int(rng.integers(1, 12))
        lengths = rng.random(nseg) * rng.choice([0, 1, 1e-6, 1e6], size=nseg)
        extra = int(rng.integers(0, 200))
        counts = _distribute(lengths, extra)
        assert counts.sum() == extra
        assert np.all(counts >= 0)
        if lengths.sum() > 0:
            assert np.all(counts[lengths == 0] == 0)
            # no share is more than one away from its quota
            quota = extra * lengths / lengths.sum()
            assert np.all(np.abs(counts - quota) < 1 + 1e-6)


def test_polyline_distribute_many_quotas_just_below_integers():
    # Many quotas are just below an integer: together they must not give more than ``extra``.
    m = 10_000
    extra = 1_000_000_000
    lengths = np.full(m + 1, 100_000 - 2e-4)
    lengths[-1] = 2.0
    quota = extra * lengths / lengths.sum()

    counts = _distribute(lengths, extra)
    assert counts.sum() == extra
    assert np.all(counts >= 0)
    assert np.all(np.abs(counts - quota) < 1 + 1e-6)


def test_polyline_distribute_zero_length_tail_never_negative():
    lengths = np.r_[np.full(10_000, 100_000 - 2e-4), 2.0, 0.0, 0.0]
    extra = 1_000_000_000
    counts = _distribute(lengths, extra)
    assert counts.sum() == extra
    assert np.all(counts >= 0)
    assert counts[-2:].tolist() == [0, 0]
    quota = extra * lengths / lengths.sum()
    assert np.all(np.abs(counts - quota) < 1 + 1e-6)


def test_polyline_distribute_cumulative_upward_rounding_many_segments():
    m = 1_000_000
    lengths = np.full(m + 1, 1.999999)
    lengths[-1] = 1.0
    counts = _distribute(lengths, 2 * m)
    assert counts.sum() == 2 * m
    assert counts.min() >= 0
    assert counts.max() <= 2


def test_polyline_distribute_random_property():
    rng = np.random.default_rng(20241021)
    for i in range(300):
        nseg = int(rng.integers(1, 60))
        kind = i % 4
        if kind == 0:
            lengths = rng.random(nseg)
        elif kind == 1:
            # quotas that are integers up to floating point noise
            ints = rng.integers(0, 6, size=nseg).astype(float)
            lengths = ints * (1 + rng.choice([0.0, 1e-16, -1e-16, 1e-15], size=nseg))
        elif kind == 2:
            # quotas just below/above integers by tiny but real amounts
            ints = rng.integers(1, 6, size=nseg).astype(float)
            lengths = ints + rng.choice([0.0, 1e-13, -1e-13, 1e-9, -1e-9], size=nseg)
        else:
            lengths = rng.random(nseg) * rng.choice([0.0, 1.0, 1e-8, 1e8], size=nseg)

        if kind in (1, 2):
            extra = int(round(lengths.sum())) + int(rng.integers(0, 3))
        else:
            extra = int(rng.integers(0, 500))

        counts = _distribute(lengths, extra)
        assert counts.sum() == extra
        assert np.all(counts >= 0)
        total = lengths.sum()
        if total > 0:
            assert np.all(counts[lengths == 0] == 0)
            quota = extra * lengths / total
            assert np.all(np.abs(counts - quota) < 1 + 1e-6)
        else:
            assert counts.max() - counts.min() <= 1


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("longer", [1e-9, 1e-8, 1e-6])
def test_sample_polyline_near_but_unequal_quotas(method, longer):
    # the second segment is slightly longer, so it gets the only extra point
    lon_wp = [0, 1, 2 + longer]
    lat, lon = sample_polyline([0, 0, 0], lon_wp, n=4, method=method)
    assert _counts_between(lon, lon_wp) == [0, 1]
    assert lon[2] == pytest.approx(1 + (1 + longer) / 2, abs=1e-9)

    # and the first one when it is the longer one
    lon_wp = [0, 1 + longer, 2 + longer]
    lat, lon = sample_polyline([0, 0, 0], lon_wp, n=4, method=method)
    assert _counts_between(lon, lon_wp) == [1, 0]


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_repeated_waypoints(method):
    # the zero-length segment gets no points
    lat, lon = sample_polyline([0, 0, 0], [0, 0, 10], n=5, method=method)
    assert lat.shape == (5,)
    np.testing.assert_array_equal(lon[:2], [0, 0])
    np.testing.assert_allclose(lon[2:], [10 / 3, 20 / 3, 10], atol=1e-9)
    np.testing.assert_allclose(lat, 0, atol=1e-12)

    # repeated waypoint in the middle and at the end
    lat, lon = sample_polyline([0, 0, 0, 0], [0, 10, 10, 10], n=6, method=method)
    np.testing.assert_allclose(lon, [0, 10 / 3, 20 / 3, 10, 10, 10], atol=1e-9)


@pytest.mark.parametrize("method", METHODS)
def test_sample_polyline_all_zero_length(method):
    lat, lon = sample_polyline([12.5, 12.5, 12.5, 12.5], [30, 30, 30, 30], n=9, method=method)
    assert lat.shape == (9,) and lon.shape == (9,)
    np.testing.assert_allclose(lat, 12.5, atol=1e-12)
    np.testing.assert_allclose(lon, 30, atol=1e-9)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))
    assert lat[0] == 12.5 and lon[0] == 30
    assert lat[-1] == 12.5 and lon[-1] == 30


# ---------------------------------------------------------------------------
# longitude range of the interior points
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["great_circle", "rhumb"])
def test_sample_polyline_lon_range_0_360(method):
    # all non-negative inputs: the interior is in [0, 360)
    lat, lon = sample_polyline([0, 10], [350, 370], n=9, method=method)
    interior = lon[1:-1]
    assert np.all(interior >= 0) and np.all(interior < 360)
    # waypoints are copied, never normalised
    assert lon[0] == 350 and lon[-1] == 370
    # the points near the start are not wrapped to negative values
    np.testing.assert_allclose(interior[-1] % 360, interior[-1])


@pytest.mark.parametrize("method", ["great_circle", "rhumb"])
def test_sample_polyline_lon_range_mixed(method):
    lat, lon = sample_polyline([0, 10, 0], [-170, 175, 200], n=41, method=method)
    interior = np.delete(lon, [0, np.flatnonzero(lon == 175)[0], len(lon) - 1])
    assert np.all(interior >= -180) and np.all(interior < 180)
    assert lon[0] == -170 and lon[-1] == 200


@pytest.mark.parametrize("method", ["great_circle", "rhumb"])
def test_sample_polyline_lon_range_waypoints_not_normalised(method):
    lat, lon = sample_polyline([0, 0], [400, 500], n=3, method=method)
    assert lon[0] == 400 and lon[-1] == 500
    # 450 normalised to 90
    np.testing.assert_allclose(lon[1], 90, atol=1e-9)

    lat, lon = sample_polyline([0, 0], [-300, -100], n=3, method="rhumb")
    assert lon[0] == -300 and lon[-1] == -100
    # -200 normalised to 160
    np.testing.assert_allclose(lon[1], 160, atol=1e-9)


def test_sample_polyline_lon_range_excludes_180():
    # the interior point at 180 is returned as -180 (range [-180, 180)) for mixed-sign inputs
    lat, lon = sample_polyline([0, 0], [-10, 370], n=3, method="rhumb")
    assert lon[0] == -10 and lon[-1] == 370
    assert lon[1] == -180


def test_sample_polyline_lon_range_linear_not_normalised():
    lat, lon = sample_polyline([0, 0], [-100, 700], n=5, method="linear")
    np.testing.assert_allclose(lon, [-100, 100, 300, 500, 700])


@pytest.mark.parametrize(
    "lon,minimum,expected",
    [
        (0.0, 0.0, 0.0),
        (360.0, 0.0, 0.0),
        (-1e-20, 0.0, 0.0),
        (-1.0, 0.0, 359.0),
        (725.0, 0.0, 5.0),
        (180.0, -180.0, -180.0),
        (-180.0, -180.0, -180.0),
        (190.0, -180.0, -170.0),
        (-190.0, -180.0, 170.0),
        (-180.0 - 1e-20, -180.0, -180.0),
        (179.9999999999, -180.0, 179.9999999999),
    ],
)
def test_polyline_wrap_lon(lon, minimum, expected):
    out = _wrap_lon(np.array([lon]), minimum)
    assert out[0] == pytest.approx(expected, abs=1e-12)
    assert minimum <= out[0] < minimum + 360


# ---------------------------------------------------------------------------
# errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(
    "lats,lons,n,match",
    [
        ([[0, 1], [2, 3]], [0, 1], 5, "sample_polyline: .*1D"),
        ([0, 1], [[0, 1], [2, 3]], 5, "sample_polyline: .*1D"),
        (0.0, 1.0, 5, "sample_polyline: .*1D"),
        ([0, 1, 2], [0, 1], 5, "sample_polyline: .*same length"),
        ([0], [0], 5, "sample_polyline: .*at least 2 waypoints"),
        ([], [], 5, "sample_polyline: .*at least 2 waypoints"),
        ([0, 91], [0, 1], 5, "sample_polyline: .*latitudes"),
        ([-90.5, 0], [0, 1], 5, "sample_polyline: .*latitudes"),
        ([0, np.nan], [0, 1], 5, "sample_polyline: .*finite"),
        ([0, 1], [0, np.nan], 5, "sample_polyline: .*finite"),
        ([0, np.inf], [0, 1], 5, "sample_polyline: .*finite"),
        ([0, 1], [0, -np.inf], 5, "sample_polyline: .*finite"),
        ([0, 1, 2], [0, 1, 2], 2, "sample_polyline: n must be at least"),
        ([0, 1], [0, 1], 1, "sample_polyline: n must be at least"),
        ([0, 1], [0, 1], 0, "sample_polyline: n must be at least"),
        ([0, 1], [0, 1], -3, "sample_polyline: n must be at least"),
        ([0, 1], [0, 1], 5.0, "sample_polyline: n must be an int"),
        ([0, 1], [0, 1], "5", "sample_polyline: n must be an int"),
        ([0, 1], [0, 1], None, "sample_polyline: n must be an int"),
        ([0, 1], [0, 1], True, "sample_polyline: n must be an int"),
    ],
)
def test_sample_polyline_invalid_input(method, lats, lons, n, match):
    with pytest.raises(ValueError, match=match):
        sample_polyline(lats, lons, n=n, method=method)


@pytest.mark.parametrize("method", ["foo", "", "greatcircle", "Rhumb", None, 1])
def test_sample_polyline_unknown_method(method):
    with pytest.raises(ValueError, match="sample_polyline: unknown method"):
        sample_polyline([0, 1], [0, 1], n=5, method=method)


@pytest.mark.parametrize(
    "lats,lons",
    [
        ([0, 0], [0, 180]),
        ([0, 0], [-90, 90]),
        ([30, -30], [10, -170]),
        ([90, -90], [0, 0]),
    ],
)
def test_sample_polyline_great_circle_antipodal(lats, lons):
    with pytest.raises(ValueError, match="sample_polyline: .*antipodal"):
        sample_polyline(lats, lons, n=10, method="great_circle")


def test_sample_polyline_great_circle_antipodal_in_later_segment_and_without_extra_points():
    with pytest.raises(ValueError, match="sample_polyline: .*segment 1"):
        sample_polyline([0, 0, 0], [0, 10, 190], n=3, method="great_circle")


@pytest.mark.parametrize("method", ["rhumb", "linear"])
def test_sample_polyline_antipodal_ok_for_other_methods(method):
    lat, lon = sample_polyline([0, 0], [0, 180], n=5, method=method)
    np.testing.assert_allclose(lon, [0, 45, 90, 135, 180])


@pytest.mark.parametrize(
    "lats,lons",
    [
        ([90, 0], [0, 0]),
        ([0, 90], [0, 0]),
        ([-90, 0], [0, 10]),
        ([0, -90], [0, 10]),
        ([0, 90, 0], [0, 10, 20]),
    ],
)
def test_sample_polyline_rhumb_pole(lats, lons):
    with pytest.raises(ValueError, match="sample_polyline: .*poles"):
        sample_polyline(lats, lons, n=10, method="rhumb")


@pytest.mark.parametrize("method", ["great_circle", "linear"])
def test_sample_polyline_pole_ok_for_other_methods(method):
    lat, lon = sample_polyline([0, 90], [10, 0], n=4, method=method)
    assert lat.shape == (4,)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))
    assert (lat[-1], lon[-1]) == (90, 0)
