# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Tests for the extension surface of :mod:`earthkit.geo.grids.xarray`.

Covers :class:`RulesCoordinateGuesser`, the :class:`CoordinateGuesser`
subclass hooks (``_guess_non_geo``, ``_scalar``, ``_unsupported``,
``_grid_cache_key``), the datetime/timedelta handling in
:func:`extract_single_value` and the ``Coordinate`` kwargs pass-through.
"""

import datetime

import numpy as np
import pytest

from earthkit.geo.grids.xarray import (
    Coordinate,
    DefaultCoordinateGuesser,
    LatitudeCoordinate,
    LongitudeCoordinate,
    PointCoordinate,
    RulesCoordinateGuesser,
    ScalarCoordinate,
    UnsupportedCoordinate,
    XCoordinate,
    YCoordinate,
    extract_single_value,
)

xr = pytest.importorskip("xarray")


# --------------------------------------------------------------------------------------
# RulesCoordinateGuesser
# --------------------------------------------------------------------------------------


def test_rules_lon_lat():
    # coordinates with names the default heuristics do not know
    ds = xr.Dataset(
        {"a": (("lat_deg", "lon_deg"), np.zeros((3, 4)))},
        coords={
            "lat_deg": ("lat_deg", np.array([50.0, 40.0, 30.0]), {"standard_name": "grid_lat"}),
            "lon_deg": ("lon_deg", np.array([0.0, 10.0, 20.0, 30.0]), {"standard_name": "grid_lon"}),
        },
    )

    # the default guesser cannot classify them
    default = DefaultCoordinateGuesser(ds)
    assert isinstance(default.guess(ds["lat_deg"], "lat_deg"), UnsupportedCoordinate)
    assert isinstance(default.guess(ds["lon_deg"], "lon_deg"), UnsupportedCoordinate)

    guess = RulesCoordinateGuesser(
        ds,
        {
            "latitude": [{"standard_name": "grid_lat"}],
            "longitude": [{"standard_name": "grid_lon"}],
        },
    )
    lat = guess.guess(ds["lat_deg"], "lat_deg")
    lon = guess.guess(ds["lon_deg"], "lon_deg")
    assert isinstance(lat, LatitudeCoordinate)
    assert isinstance(lon, LongitudeCoordinate)

    # and the grid is built from them
    coordinates = [guess.guess(ds[c], c) for c in ds["a"].coords]
    grid = guess.grid(coordinates, ds["a"])
    lats, lons = grid.latlon_grid.latlons
    assert lats.shape == (12,)
    assert set(lats.tolist()) == {50.0, 40.0, 30.0}
    assert set(lons.tolist()) == {0.0, 10.0, 20.0, 30.0}


def test_rules_x_y_point():
    ds = xr.Dataset(
        {
            "f": (("ix", "iy"), np.zeros((3, 2))),
            "p": (("station",), np.zeros((2,))),
        },
        coords={
            "east": ("ix", np.array([0.0, 1.0, 2.0]), {"axis": "X", "long_name": "east offset"}),
            "north": ("iy", np.array([0.0, 1.0]), {"axis": "Y", "long_name": "north offset"}),
            "site": ("station", np.array([1, 2]), {"standard_name": "station_index"}),
        },
    )

    guess = RulesCoordinateGuesser(
        ds,
        {
            "x": [{"axis": "X", "long_name": "east offset"}],
            "y": [{"long_name": "north offset"}],
            "point": [{"standard_name": "station_index"}],
        },
    )
    assert isinstance(guess.guess(ds["east"], "east"), XCoordinate)
    assert isinstance(guess.guess(ds["north"], "north"), YCoordinate)
    assert isinstance(guess.guess(ds["site"], "site"), PointCoordinate)

    # the first matching rule wins; a single dict is accepted instead of a
    # list; non-string rule values are not attribute constraints
    guess2 = RulesCoordinateGuesser(
        ds,
        {
            "x": [{"long_name": "nope"}, {"axis": "X"}],
            "y": {"axis": "Y", "levtype": ["pl", "sfc"]},
        },
    )
    assert isinstance(guess2.guess(ds["east"], "east"), XCoordinate)
    assert isinstance(guess2.guess(ds["north"], "north"), YCoordinate)


def test_rules_ignores_non_geo_keys_and_misses():
    ds = xr.Dataset(
        {"a": (("lat", "lon"), np.zeros((3, 4)))},
        coords={
            "lat": ("lat", np.array([50.0, 40.0, 30.0]), {"standard_name": "latitude"}),
            "lon": ("lon", np.array([0.0, 10.0, 20.0, 30.0]), {"standard_name": "longitude"}),
            "time": ("time", np.array([0, 1, 2], dtype="datetime64[ns]")),
        },
    )

    # only the geography keys are consulted here: "time" is a non-geo key
    # (handled by a subclass on the earthkit-data side), and with no rules
    # for latitude/longitude the CF-convention heuristics are NOT applied
    guess = RulesCoordinateGuesser(ds, {"time": [{"standard_name": "time"}]})
    assert isinstance(guess.guess(ds["lat"], "lat"), UnsupportedCoordinate)
    assert isinstance(guess.guess(ds["lon"], "lon"), UnsupportedCoordinate)
    assert isinstance(guess.guess(ds["time"], "time"), UnsupportedCoordinate)

    # no rules at all: everything falls back
    empty = RulesCoordinateGuesser(ds)
    assert isinstance(empty.guess(ds["lat"], "lat"), UnsupportedCoordinate)


# --------------------------------------------------------------------------------------
# CoordinateGuesser subclass hooks
# --------------------------------------------------------------------------------------


class TimeLikeCoordinate(Coordinate):
    is_time = True
    mars_names = ("valid_datetime",)


class MyScalarCoordinate(ScalarCoordinate):
    pass


class MyUnsupportedCoordinate(UnsupportedCoordinate):
    pass


class HookGuesser(DefaultCoordinateGuesser):
    """Classifies "time" through the non-geo hook and uses custom fallbacks."""

    def _guess_non_geo(self, c, attributes):
        if attributes.name == "time":
            return TimeLikeCoordinate(c)
        return None

    def _scalar(self, c):
        return MyScalarCoordinate(c)

    def _unsupported(self, c):
        return MyUnsupportedCoordinate(c)


def _ds_with_time():
    return xr.Dataset(
        {"a": (("time", "lat", "lon"), np.zeros((2, 3, 4)))},
        coords={
            "time": ("time", np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")),
            "lat": ("lat", np.array([50.0, 40.0, 30.0]), {"standard_name": "latitude"}),
            "lon": ("lon", np.array([0.0, 10.0, 20.0, 30.0]), {"standard_name": "longitude"}),
            "reference": ((), np.array(2020)),
            "labels": ("time", np.array(["x", "y"])),
        },
    )


def test_guess_non_geo_hook_and_grid():
    ds = _ds_with_time()
    guess = HookGuesser(ds)

    time = guess.guess(ds["time"], "time")
    assert isinstance(time, TimeLikeCoordinate)
    assert time.is_time and time.is_dim
    # geography matchers keep priority over the hook
    assert isinstance(guess.guess(ds["lat"], "lat"), LatitudeCoordinate)
    # the fallback factories are used for the rest
    assert isinstance(guess.guess(ds["reference"], "reference"), MyScalarCoordinate)
    assert isinstance(guess.guess(ds["labels"], "labels"), MyUnsupportedCoordinate)

    # the grid of the variable still builds with the extra coordinate
    coordinates = [guess.guess(ds[c], c) for c in ds["a"].coords]
    grid = guess.grid(coordinates, ds["a"])
    lats, lons = grid.latlon_grid.latlons
    assert lats.shape == (12,)
    assert set(lats.tolist()) == {50.0, 40.0, 30.0}
    assert set(lons.tolist()) == {0.0, 10.0, 20.0, 30.0}

    # the default guesser does not classify "time" (and reports it unsupported),
    # while the base factories build the base classes
    default = DefaultCoordinateGuesser(ds)
    assert isinstance(default.guess(ds["time"], "time"), UnsupportedCoordinate)
    assert isinstance(default.guess(ds["reference"], "reference"), ScalarCoordinate)
    assert isinstance(default.guess(ds["labels"], "labels"), UnsupportedCoordinate)


def test_coordinate_flags_and_kwargs_pass_through():
    class LevelLikeCoordinate(Coordinate):
        def __init__(self, variable, levtype):
            super().__init__(variable, levtype=levtype)
            self.levtype = levtype

    ds = xr.Dataset(
        {"a": (("plev", "lat"), np.zeros((2, 3)))},
        coords={
            "plev": ("plev", np.array([700.0, 500.0])),
            "lat": ("lat", np.array([50.0, 40.0, 30.0])),
        },
    )

    c = LevelLikeCoordinate(ds["plev"], "pl")
    assert c.kwargs == {"levtype": "pl"}
    # base class flags, including the new ones
    assert (c.is_step, c.is_date, c.is_member) == (False, False, False)
    assert not c.is_time
    # reduced() carries the extra kwargs over
    r = c.reduced(1)
    assert isinstance(r, LevelLikeCoordinate)
    assert r.levtype == "pl"
    assert float(r.variable.values) == 500.0

    plain = Coordinate(ds["lat"])
    assert plain.kwargs == {}
    assert (plain.is_step, plain.is_date, plain.is_member) == (False, False, False)


def test_grid_cache_key_default_and_override():
    ds = xr.Dataset(
        {
            "a": (("lat", "lon"), np.zeros((3, 4)), {"grid_mapping": "rotated"}),
            "b": (("lat", "lon"), np.ones((3, 4)), {"grid_mapping": "rotated"}),
            "rotated": (
                (),
                np.array(0),
                {
                    "grid_mapping_name": "rotated_latitude_longitude",
                    "grid_north_pole_latitude": 40.0,
                    "grid_north_pole_longitude": -80.0,
                },
            ),
        },
        coords={
            "lat": ("lat", np.array([50.0, 40.0, 30.0]), {"standard_name": "latitude"}),
            "lon": ("lon", np.array([0.0, 10.0, 20.0, 30.0]), {"standard_name": "longitude"}),
        },
    )

    guess = DefaultCoordinateGuesser(ds)
    lat = LatitudeCoordinate(ds["lat"])
    lon = LongitudeCoordinate(ds["lon"])
    dim_vars = ("lat", "lon")

    # the default key is exactly what the cache used inline
    key = guess._grid_cache_key(lat, lon, dim_vars, ds["a"])
    assert key == (lat.name, lon.name, dim_vars, guess._mapping_key(ds["a"]))

    # with the default key, both variables share one grid (cache hit)
    g_a = guess.grid([guess.guess(ds["lat"], "lat"), guess.guess(ds["lon"], "lon")], ds["a"])
    g_b = guess.grid([guess.guess(ds["lat"], "lat"), guess.guess(ds["lon"], "lon")], ds["b"])
    assert g_a.latlon_grid is g_b.latlon_grid
    # ... and the built grid is cached under the default key
    assert key in guess._grid_cache
    assert guess._grid_cache[key] is g_a.latlon_grid

    # a subclass whose grids depend on anything else from the variable can
    # extend the key to stop the sharing
    class KeyedByVariableGuesser(DefaultCoordinateGuesser):
        def _grid_cache_key(self, first, second, dim_vars, variable):
            return super()._grid_cache_key(first, second, dim_vars, variable) + (variable.name,)

    keyed = KeyedByVariableGuesser(ds)
    g_a = keyed.grid([keyed.guess(ds["lat"], "lat"), keyed.guess(ds["lon"], "lon")], ds["a"])
    g_b = keyed.grid([keyed.guess(ds["lat"], "lat"), keyed.guess(ds["lon"], "lon")], ds["b"])
    assert g_a.latlon_grid is not g_b.latlon_grid


# --------------------------------------------------------------------------------------
# extract_single_value: datetime64 / timedelta64
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value,expected",
    [
        (np.datetime64("2021-01-01T12:00:00"), datetime.datetime(2021, 1, 1, 12, 0, 0)),
        (np.datetime64("2020-03-04T05:06:07"), datetime.datetime(2020, 3, 4, 5, 6, 7)),
        (np.datetime64("2020-03-04"), datetime.datetime(2020, 3, 4, 0, 0, 0)),
    ],
)
def test_extract_single_value_datetime_0d(value, expected):
    da = xr.DataArray(value)
    got = extract_single_value(da)
    assert isinstance(got, datetime.datetime)
    assert got == expected


def test_extract_single_value_datetime_len1_and_nat():
    da = xr.DataArray(np.array(["2020-03-04T05:06:07"], dtype="datetime64[s]"), dims="time")
    got = extract_single_value(da)
    assert got == datetime.datetime(2020, 3, 4, 5, 6, 7)

    nat = extract_single_value(xr.DataArray(np.datetime64("NaT", "s")))
    assert nat is None


def test_extract_single_value_timedelta():
    da = xr.DataArray(np.timedelta64(6, "h"))
    got = extract_single_value(da)
    assert isinstance(got, datetime.timedelta)
    assert got == datetime.timedelta(hours=6)

    da = xr.DataArray(np.timedelta64(90, "m"))
    assert extract_single_value(da) == datetime.timedelta(minutes=90)

    da = xr.DataArray(np.array([np.timedelta64(3, "h")], dtype="timedelta64[h]"), dims="step")
    assert extract_single_value(da) == datetime.timedelta(hours=3)

    assert extract_single_value(xr.DataArray(np.timedelta64("NaT", "s"))) is None


def test_extract_single_value_plain():
    # unchanged behaviour for other dtypes
    assert extract_single_value(xr.DataArray(np.array([3.5]))) == 3.5
    assert extract_single_value(xr.DataArray(np.array(7))) == 7
