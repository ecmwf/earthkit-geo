# (C) Copyright 2024 Anemoi contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import logging
import warnings

import numpy as np
import pytest

from earthkit.geo.grids._regrid.data.xarray.coordinates import Coordinate
from earthkit.geo.grids._regrid.data.xarray.grid import (
    MeshedGrid,
    MeshProjectionGrid,
    UnstructuredGrid,
    UnstructuredProjectionGrid,
)
from earthkit.geo.grids._regrid.data.xarray.grid_spec import (
    healpix_grid_spec_from_points,
    reduced_gg_grid_spec_from_points,
)

xr = pytest.importorskip("xarray")
eckit_geo = pytest.importorskip("eckit.geo")

INCREMENTS = [0.05, 0.08333333333333333, 0.1, 0.125, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 2.5, 3.0, 5.0, 10.0]


def _coordinate(name, values):
    return Coordinate(xr.DataArray(np.asarray(values), dims=[name], name=name))


def _meshed_grid(lat_values, lon_values):
    lat = _coordinate("latitude", lat_values)
    lon = _coordinate("longitude", lon_values)
    return MeshedGrid(lat, lon, variable_dims=("latitude", "longitude"))


def _unstructured_grid(lat_values, lon_values, dims):
    lat = Coordinate(xr.DataArray(np.asarray(lat_values), dims=dims, name="latitude"))
    lon = Coordinate(xr.DataArray(np.asarray(lon_values), dims=dims, name="longitude"))
    return UnstructuredGrid(lat, lon, variable_dims=dims)


@pytest.mark.parametrize("inc", INCREMENTS)
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_subarea(inc, ascending):
    npts = 4
    north = 10.0
    west = 20.0

    lat_values = north - inc * np.arange(npts)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = west + inc * np.arange(npts)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])

    if ascending:
        assert spec.get("order") == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (npts, npts)
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


@pytest.mark.parametrize("inc", INCREMENTS)
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_global(inc, ascending):
    reference = eckit_geo.Grid({"grid": [inc, inc]})
    lat_values = np.asarray(reference.distinct_latitudes())
    lon_values = np.asarray(reference.distinct_longitudes())

    if ascending:
        lat_values = lat_values[::-1]

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])

    if ascending:
        assert spec.get("order") == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    expected_keys = {"grid", "order"} if ascending else {"grid"}

    # a grid spanning the full globe is normalised back to the simple global form,
    # i.e. the "area" key is dropped
    assert set(eckit_grid.spec.keys()) == expected_keys
    assert eckit_grid.spec["grid"] == pytest.approx([inc, inc])
    if ascending:
        assert eckit_grid.spec["order"] == "i+j+"
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


@pytest.mark.parametrize("inc", INCREMENTS)
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_subarea_shifted(inc, ascending):
    npts = 4
    lon_offset = 0.3 * inc
    lat_offset = 0.6 * inc

    # a corner aligned on the default (0, 0) reference, plus a deliberate sub-increment shift
    north = round(50.0 / inc) * inc + lat_offset
    west = round(-10.0 / inc) * inc + lon_offset

    lat_values = north - inc * np.arange(npts)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = west + inc * np.arange(npts)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])
    assert spec["reference"] == pytest.approx([lon_offset, lat_offset])

    if ascending:
        assert spec.get("order") == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert set(eckit_grid.spec.keys()) == set(spec.keys())
    assert eckit_grid.spec["grid"] == pytest.approx(spec["grid"])
    assert eckit_grid.spec["area"] == pytest.approx(spec["area"])
    assert eckit_grid.spec["reference"] == pytest.approx(spec["reference"])
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (npts, npts)
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


@pytest.mark.parametrize("inc", INCREMENTS)
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_global_shifted(inc, ascending):
    lon_offset = 0.3 * inc
    lat_offset = 0.6 * inc

    reference = eckit_geo.Grid({"grid": [inc, inc], "reference": [lon_offset, lat_offset]})
    lat_values = np.asarray(reference.distinct_latitudes())
    lon_values = np.asarray(reference.distinct_longitudes())

    if ascending:
        lat_values = lat_values[::-1]

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["reference"] == pytest.approx([lon_offset, lat_offset])

    if ascending:
        assert spec.get("order") == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.spec["grid"] == pytest.approx([inc, inc])
    assert eckit_grid.spec["reference"] == pytest.approx([lon_offset, lat_offset])
    assert eckit_grid.spec.get("order") == ("i+j+" if ascending else None)

    # re-normalising an already-normalised spec is a no-op
    assert eckit_geo.Grid(eckit_grid.spec).spec == eckit_grid.spec

    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_fine_subarea(ascending):
    # a 0.2 degree grid with 899 latitudes (-89.8 to 89.8) and 1800 longitudes
    # (-180 to 179.8); it excludes the poles, so it is a subarea, stored either
    # with latitude ascending (grid A, "i+j+") or descending (grid B, "i+j-")
    inc = 0.2
    north = 89.8
    west = -180.0
    nlat = 899
    nlon = 1800

    lat_values = north - inc * np.arange(nlat)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = west + inc * np.arange(nlon)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])
    # the area corners lie on the default (0, 0) reference
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)

    # the points are numbered in the order of the data: the first row is
    # lat_values[0], running west to east
    lat_points, lon_points = eckit_grid.to_latlons()
    assert lat_points[0] == pytest.approx(lat_values[0])
    assert lon_points[0] == pytest.approx(lon_values[0])
    assert lat_points[-1] == pytest.approx(lat_values[-1])
    assert lon_points[-1] == pytest.approx(lon_values[-1])


@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_fine_subarea_float32(ascending):
    # the same 0.2 degree subarea as test_meshed_grid_spec_fine_subarea, but
    # stored in float32, as in a real netCDF file: the coordinates are then
    # quantised to an ulp of the coordinate magnitude (about 1e-5 degrees),
    # so the steps are not uniform to within np.allclose's default tolerance
    inc = 0.2
    north = 89.8
    west = -180.0
    nlat = 899
    nlon = 1800

    lat_values = (north - inc * np.arange(nlat)).astype(np.float32)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (west + inc * np.arange(nlon)).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])
    # the area corners lie on the default (0, 0) reference; the float32
    # quantisation must not produce a spurious "reference"
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)


@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_float32_off_lattice(ascending):
    # a small float32 subarea whose corners are deliberately off the increment
    # lattice (10.03125 is not a multiple of 0.25). Its values are exactly
    # representable, so the 6 digit "cleaning" a float32 axis goes through must
    # refuse them and keep the raw corners: had it not, the rounded north and south
    # would no longer span exactly (nlat - 1) increments and the grid would lose a
    # row of points
    lat_values = (10.03125 - 0.25 * np.arange(4)).astype(np.float32)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (99.96875 + 0.25 * np.arange(4)).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([0.25, 0.25])
    assert spec["area"] == pytest.approx([10.03125, 99.96875, 9.28125, 100.71875])

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.shape == (4, 4)

    # the points are the ones of the axis, to within a float32 ulp
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending), rtol=0, atol=1e-5)
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values), rtol=0, atol=1e-5)


def test_meshed_grid_spec_float32_irregular_drift():
    # a float32 axis every step of which is within a few ulps of the average step,
    # but whose points wander off the regular lattice by up to 0.01 degrees: only
    # quantisation noise may be forgiven, a genuinely irregular axis is not a grid
    lat_values = np.array([50.0, 49.0, 48.0], dtype=np.float32)
    t = np.linspace(0, np.pi, 1001)
    lon_values = (np.linspace(100, 110, 1001) + 0.01 * np.sin(t)).astype(np.float32)

    # every step is within 4 float32 ulps of the average one, so a per-step
    # tolerance alone cannot reject this axis
    steps = np.diff(lon_values)
    assert np.max(np.abs(steps - np.mean(steps))) <= 4 * np.finfo(np.float32).eps * 110

    grid = _meshed_grid(lat_values, lon_values)
    assert grid.grid_spec is None


@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_float32_twelfths(ascending):
    # the grid of latlon_1_per_12_ocean_nopole_jplus.nc: a 1/12 degree float32 grid
    # from -80 to 90 (2041 latitudes) and -180 to 179.9167 (4320 longitudes). Its
    # increment has more than 6 significant digits, so the 6 digit cleaning cannot
    # recover it; the axis is instead snapped onto the lattice of the fraction 1/12,
    # and the gridspec carries the increment at the 12 decimals eckit uses
    inc = round(1 / 12, 12)
    nlat, nlon = 2041, 4320

    lat_values = (90.0 - np.arange(nlat) / 12).astype(np.float32)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (-180.0 + np.arange(nlon) / 12).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == [inc, inc]
    assert spec["area"] == pytest.approx([90.0, -180.0, -80.0, round(-180.0 + (nlon - 1) / 12, 12)], abs=1e-9)
    # the corners lie on the default (0, 0) reference: no spurious "reference"
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)

    # the grid's coordinates are the file's, to within a float32 ulp
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending), rtol=0, atol=1e-5)
    assert np.allclose(eckit_grid.distinct_longitudes(), lon_values, rtol=0, atol=1e-5)


def test_meshed_grid_spec_float32_simple_fraction():
    # the snapping is not tied to 1/12: a 1/6 degree grid, whose increment rounds
    # up at the 12th decimal, is detected the same way and keeps its point count
    inc = round(1 / 6, 12)
    nlat, nlon = 1021, 2160

    lat_values = (-80.0 + np.arange(nlat) / 6).astype(np.float32)
    lon_values = (-180.0 + np.arange(nlon) / 6).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == [inc, inc]
    assert eckit_geo.Grid(spec).shape == (nlat, nlon)


def test_meshed_grid_spec_float32_wandering_twelfths():
    # an axis that merely averages out near 1/12, its points wandering +/-0.0005
    # degrees along it, is not snapped onto the 1/12 lattice: the snapping may
    # forgive float32 quantisation (a few 1e-5 degrees), not genuine wander
    n = 4320
    wander = 5e-4 * np.sin(2 * np.pi * np.arange(n) / 100)
    lat_values = (90.0 - np.arange(2041) / 12).astype(np.float32)
    lon_values = (-180.0 + np.arange(n) / 12 + wander).astype(np.float32)

    # the latitude axis is a clean 1/12 lattice: it is the wandering longitude
    # axis that must be rejected, and with it the whole grid
    grid = _meshed_grid(lat_values, lon_values)
    assert grid.grid_spec is None


def test_meshed_grid_spec_float32_increment_not_a_fraction():
    # a float32 grid written on a 0.08334 degree increment, which sits within 1e-4
    # relative of 1/12 but is not one, must not be pulled onto the 1/12 lattice:
    # it keeps its own increment (its steps are uniform enough that the raw fit
    # accepts it, but the snapping must stay out of it)
    n = 2041
    lat_values = (90.0 - np.arange(n) / 12).astype(np.float32)
    lon_values = (-180.0 + 0.08334 * np.arange(n)).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([0.08334, round(1 / 12, 12)], abs=1e-9)


@pytest.mark.parametrize("past", ["north", "south"])
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_float32_lat_past_pole(past, ascending):
    # a float32 latitude axis reaching a whole degree past a pole is not a grid:
    # the gridspec fitting it would be clamped by eckit into fewer rows than the
    # axis has. The guard is lenient with a few ulps of float32 noise (a degree
    # near the poles is quantised to ~7.6e-6), but one degree is ~130000 of them
    if past == "north":
        lat_values = (91.0 - 0.5 * np.arange(364)).astype(np.float32)  # 91 down to -90.5
    else:
        lat_values = (-91.0 + 0.5 * np.arange(364)).astype(np.float32)  # -91 up to 90.5
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (0.5 * np.arange(720)).astype(np.float32)

    # the axis is exactly regular in float32 (all values are multiples of 0.5), so
    # only the latitude-bounds guard can reject it
    grid = _meshed_grid(lat_values, lon_values)
    assert grid.grid_spec is None


@pytest.mark.parametrize("inc", [0.5, 0.25])
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_float32_global_at_poles(inc, ascending):
    # a global float32 axis ends exactly at +-90, which is exact in float32, so the
    # bounds guard accepts it: the axis is 90 - inc * arange(180/inc + 1), reaching
    # from pole to pole, stored descending ("i+j-") or ascending ("i+j+")
    nlat = int(round(180 / inc)) + 1
    nlon = int(round(360 / inc))

    lat_values = (90.0 - inc * np.arange(nlat)).astype(np.float32)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (-180.0 + inc * np.arange(nlon)).astype(np.float32)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([90.0, -180.0, -90.0, 180.0 - inc])
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    # the grid built from the spec has as many points as the axes
    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)


@pytest.mark.parametrize("end", ["north", "south"])
@pytest.mark.parametrize("ascending", [False, True])
def test_meshed_grid_spec_float32_global_poles_within_noise(end, ascending):
    # a global float32 axis whose far end sits a couple of ulps past a pole, as
    # axis quantisation noise can leave it, is still a grid: the guard forgives a
    # few ulps, and the 6 digit cleaning brings the spec's ends back to exactly +-90
    inc = 0.5
    nlat = 361
    nlon = 720
    ulp = np.spacing(np.float32(90.0))  # the ulp of a float32 near 90

    # the whole axis nudged 2 ulps past one pole: the other end stays within bounds
    shift = 2 * ulp if end == "north" else -2 * ulp
    lat_values = (90.0 + shift - inc * np.arange(nlat)).astype(np.float32)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = (-180.0 + inc * np.arange(nlon)).astype(np.float32)

    # the far end really does reach past its pole, by more than a single ulp
    if end == "north":
        assert float(np.max(lat_values)) > 90.0 + ulp
    else:
        assert float(np.min(lat_values)) < -90.0 - ulp

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([90.0, -180.0, -90.0, 180.0 - inc])
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)


def test_meshed_grid_spec_specific_area():
    north, west, south, east = 73, -27, 33, 45
    inc = 4

    lat_values = np.arange(north, south - inc / 2, -inc)
    lon_values = np.arange(west, east + inc / 2, inc)

    grid = _meshed_grid(lat_values, lon_values)
    spec = grid.grid_spec

    assert spec == {"grid": [inc, inc], "area": [north, west, south, east], "reference": [1.0, 1.0]}

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.spec == spec
    assert "order" not in eckit_grid.spec
    assert eckit_grid.order == "i+j-"
    assert eckit_grid.shape == (len(lat_values), len(lon_values))
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=True))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


def test_meshed_grid_spec_irregular():
    grid = _meshed_grid([10.0, 9.0, 7.0, 6.0], [0.0, 1.0, 2.0, 3.0])
    assert grid.grid_spec is None


def test_meshed_grid_spec_too_few_points():
    grid = _meshed_grid([10.0], [0.0])
    assert grid.grid_spec is None


@pytest.mark.parametrize("inc", INCREMENTS)
@pytest.mark.parametrize("ascending", [False, True])
def test_unstructured_grid_spec_1d(inc, ascending):
    npts = 4
    north = 10.0
    west = 20.0

    lat_values = north - inc * np.arange(npts)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = west + inc * np.arange(npts)

    # row-major: lat slowest-varying, lon fastest-varying, as flattened "values" points
    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")
    lat_points = lat_points.flatten()
    lon_points = lon_points.flatten()

    grid = _unstructured_grid(lat_points, lon_points, dims=("values",))
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])
    assert spec.get("order") == ("i+j+" if ascending else None)

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (npts, npts)
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


def test_unstructured_grid_spec_1d_lon_major():
    # column-major storage: lon slowest-varying, lat fastest-varying
    lat_values = np.array([30.0, 20.0, 10.0])
    lon_values = np.array([0.0, 10.0, 20.0, 30.0])

    lon_points, lat_points = np.meshgrid(lon_values, lat_values, indexing="ij")
    lat_points = lat_points.flatten()
    lon_points = lon_points.flatten()

    grid = _unstructured_grid(lat_points, lon_points, dims=("values",))
    spec = grid.grid_spec

    assert spec == {"grid": [10.0, 10.0], "area": [30.0, 0.0, 10.0, 30.0]}


@pytest.mark.parametrize("ascending", [False, True])
def test_unstructured_grid_spec_2d(ascending):
    inc = 2.5
    lat_values = np.array([15.0, 12.5, 10.0, 7.5])
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = np.array([-5.0, -2.5, 0.0, 2.5])

    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")

    grid = _unstructured_grid(lat_points, lon_points, dims=("y", "x"))
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec.get("order") == ("i+j+" if ascending else None)

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert np.allclose(eckit_grid.distinct_latitudes(), sorted(lat_values, reverse=not ascending))
    assert np.allclose(eckit_grid.distinct_longitudes(), sorted(lon_values))


@pytest.mark.parametrize("ascending", [False, True])
def test_unstructured_grid_spec_2d_fine_subarea(ascending):
    # the same 0.2 degree subarea as test_meshed_grid_spec_fine_subarea, but
    # given as 2D per-point latitudes/longitudes
    inc = 0.2
    north = 89.8
    west = -180.0
    nlat = 899
    nlon = 1800

    lat_values = north - inc * np.arange(nlat)
    if ascending:
        lat_values = lat_values[::-1]
    lon_values = west + inc * np.arange(nlon)

    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")

    grid = _unstructured_grid(lat_points, lon_points, dims=("y", "x"))
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([inc, inc])
    assert spec["area"] == pytest.approx([
        float(np.max(lat_values)),
        float(np.min(lon_values)),
        float(np.min(lat_values)),
        float(np.max(lon_values)),
    ])
    assert "reference" not in spec

    if ascending:
        assert spec["order"] == "i+j+"
    else:
        assert "order" not in spec

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.type == "regular_ll"
    assert eckit_grid.order == ("i+j+" if ascending else "i+j-")
    assert eckit_grid.shape == (nlat, nlon)

    # the points are numbered in the order of the (row-major) 2D data
    ek_lat, ek_lon = eckit_grid.to_latlons()
    assert ek_lat[0] == pytest.approx(lat_points[0, 0])
    assert ek_lon[0] == pytest.approx(lon_points[0, 0])
    assert ek_lat[-1] == pytest.approx(lat_points[-1, -1])
    assert ek_lon[-1] == pytest.approx(lon_points[-1, -1])


def test_unstructured_grid_spec_global():
    inc = 10.0
    reference = eckit_geo.Grid({"grid": [inc, inc]})
    lat_values = np.asarray(reference.distinct_latitudes())
    lon_values = np.asarray(reference.distinct_longitudes())

    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")

    grid = _unstructured_grid(lat_points.flatten(), lon_points.flatten(), dims=("values",))
    spec = grid.grid_spec

    assert spec["grid"] == pytest.approx([inc, inc])
    assert "order" not in spec

    # a grid spanning the full globe is normalised back to the simple global form,
    # i.e. the "area" key is dropped
    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.spec == {"grid": [inc, inc]}


def test_unstructured_grid_spec_shifted():
    inc = 2.0
    lon_offset = 0.5
    lat_offset = 1.0

    north = 10.0 + lat_offset
    west = -10.0 + lon_offset
    lat_values = north - inc * np.arange(4)
    lon_values = west + inc * np.arange(4)

    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")

    grid = _unstructured_grid(lat_points.flatten(), lon_points.flatten(), dims=("values",))
    spec = grid.grid_spec

    assert spec == {
        "grid": [inc, inc],
        "area": [
            float(np.max(lat_values)),
            float(np.min(lon_values)),
            float(np.min(lat_values)),
            float(np.max(lon_values)),
        ],
        "reference": [lon_offset, lat_offset],
    }

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.spec == spec


def test_unstructured_grid_spec_irregular_spacing():
    lat_values = np.array([10.0, 9.0, 7.0, 6.0])
    lon_values = np.array([0.0, 1.0, 2.0, 3.0])
    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")

    grid = _unstructured_grid(lat_points.flatten(), lon_points.flatten(), dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_missing_point():
    lat_values = np.array([10.0, 9.0, 8.0])
    lon_values = np.array([0.0, 1.0, 2.0])
    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")
    lat_points = lat_points.flatten()
    lon_points = lon_points.flatten()

    # drop one point: the remaining points no longer form a complete regular grid
    lat_points = lat_points[1:]
    lon_points = lon_points[1:]

    grid = _unstructured_grid(lat_points, lon_points, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_duplicate_point():
    lat_values = np.array([10.0, 9.0, 8.0])
    lon_values = np.array([0.0, 1.0, 2.0])
    lat_points, lon_points = np.meshgrid(lat_values, lon_values, indexing="ij")
    lat_points = lat_points.flatten()
    lon_points = lon_points.flatten()

    # duplicate the first point in place of the last one
    lat_points[-1] = lat_points[0]
    lon_points[-1] = lon_points[0]

    grid = _unstructured_grid(lat_points, lon_points, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_not_a_grid():
    # purely scattered points, not a grid at all
    lat_values = [34.05, 36.12, 40.71]
    lon_values = [-118.24, -115.17, -74.01]

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_too_few_points():
    grid = _unstructured_grid([10.0, 9.0], [0.0, 1.0], dims=("values",))
    assert grid.grid_spec is None


HEALPIX_NSIDES = [2, 4, 8, 16, 32, 64]


@pytest.mark.parametrize("nside", HEALPIX_NSIDES)
@pytest.mark.parametrize("order", ["ring", "nested"])
def test_unstructured_grid_spec_healpix_1d(nside, order):
    reference = eckit_geo.Grid({"grid": f"H{nside}", "order": order})
    lat_values, lon_values = reference.to_latlons()

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    spec = grid.grid_spec

    expected = {"grid": f"H{nside}"}
    if order == "nested":
        expected["order"] = "nested"
    assert spec == expected

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.shape == (nside * nside * 12,)
    eckit_lat, eckit_lon = eckit_grid.to_latlons()
    assert np.allclose(eckit_lat, lat_values)
    assert np.allclose(eckit_lon, lon_values)


def test_healpix_grid_spec_from_points_rejects_2d():
    # a HEALPix grid is always stored as a 1D array of points, unlike a
    # regular lat/lon grid
    reference = eckit_geo.Grid({"grid": "H8"})
    lat_values, lon_values = reference.to_latlons()

    lat_values = np.asarray(lat_values).reshape(-1, 1)
    lon_values = np.asarray(lon_values).reshape(-1, 1)

    assert healpix_grid_spec_from_points(lat_values, lon_values) is None


def test_unstructured_grid_spec_healpix_wrong_size():
    nside = 4
    reference = eckit_geo.Grid({"grid": f"H{nside}", "order": "ring"})
    lat_values, lon_values = reference.to_latlons()

    # drop one point: no longer a valid HEALPix point count
    lat_values = np.asarray(lat_values)[:-1]
    lon_values = np.asarray(lon_values)[:-1]

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_healpix_scattered_points():
    # the right point count (H4 == 192 points) but not actually a HEALPix grid
    rng = np.random.default_rng(0)
    lat_values = rng.uniform(-90, 90, size=192)
    lon_values = rng.uniform(0, 360, size=192)

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_healpix_wrong_order():
    # points from a "ring"-ordered grid should not be mistaken for "nested"-only
    # matches, i.e. the detection must not match a wrong order
    reference_ring = eckit_geo.Grid({"grid": "H4", "order": "ring"})
    lat_values, lon_values = reference_ring.to_latlons()

    reference_nested = eckit_geo.Grid({"grid": "H4", "order": "nested"})
    nested_lat, nested_lon = reference_nested.to_latlons()

    assert not np.allclose(lat_values, nested_lat) or not np.allclose(lon_values, nested_lon)

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec == {"grid": "H4"}


REDUCED_GG_NUMBERS = [16, 24, 32, 48, 64, 80, 96, 128]


@pytest.mark.parametrize("n", REDUCED_GG_NUMBERS)
@pytest.mark.parametrize("kind", ["O", "N"])
def test_unstructured_grid_spec_reduced_gg_1d(n, kind):
    grid_name = f"{kind}{n}"
    reference = eckit_geo.Grid({"grid": grid_name})
    lat_values, lon_values = reference.to_latlons()

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    spec = grid.grid_spec

    assert spec == {"grid": grid_name}

    eckit_grid = eckit_geo.Grid(spec)
    assert eckit_grid.shape == (len(lat_values),)
    eckit_lat, eckit_lon = eckit_grid.to_latlons()
    assert np.allclose(eckit_lat, lat_values)
    assert np.allclose(eckit_lon, lon_values)


def test_reduced_gg_grid_spec_from_points_rejects_2d():
    # a reduced Gaussian grid has a varying number of longitudes per latitude row,
    # so it cannot be represented as a 2D array, unlike a regular lat/lon or HEALPix grid
    reference = eckit_geo.Grid({"grid": "O32"})
    lat_values, lon_values = reference.to_latlons()

    lat_values = np.asarray(lat_values).reshape(-1, 1)
    lon_values = np.asarray(lon_values).reshape(-1, 1)

    assert reduced_gg_grid_spec_from_points(lat_values, lon_values) is None


def test_unstructured_grid_spec_reduced_gg_wrong_size():
    reference = eckit_geo.Grid({"grid": "O32"})
    lat_values, lon_values = reference.to_latlons()

    # drop one point: no longer a known reduced Gaussian grid point count
    lat_values = np.asarray(lat_values)[:-1]
    lon_values = np.asarray(lon_values)[:-1]

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_reduced_gg_scattered_points():
    # the right point count (O32 == 5248 points) but not actually a reduced Gaussian grid
    rng = np.random.default_rng(0)
    lat_values = rng.uniform(-90, 90, size=5248)
    lon_values = rng.uniform(0, 360, size=5248)

    grid = _unstructured_grid(lat_values, lon_values, dims=("values",))
    assert grid.grid_spec is None


def test_unstructured_grid_spec_reduced_gg_octahedral_vs_classical():
    # O32 and N32 have different point counts, so they must not be confused
    reference_o = eckit_geo.Grid({"grid": "O32"})
    reference_n = eckit_geo.Grid({"grid": "N32"})
    assert reference_o.shape != reference_n.shape

    lat_o, lon_o = reference_o.to_latlons()
    grid_o = _unstructured_grid(lat_o, lon_o, dims=("values",))
    assert grid_o.grid_spec == {"grid": "O32"}

    lat_n, lon_n = reference_n.to_latlons()
    grid_n = _unstructured_grid(lat_n, lon_n, dims=("values",))
    assert grid_n.grid_spec == {"grid": "N32"}


# the CF grid mapping attributes of a rotated latitude/longitude grid whose
# north pole sits at (39.25, -162), as in the cordex sample
ROTATED_MAPPING = {
    "grid_mapping_name": "rotated_latitude_longitude",
    "grid_north_pole_latitude": 39.25,
    "grid_north_pole_longitude": -162.0,
    "north_pole_grid_longitude": 0.0,
}


def _rotated_grid(rlat_values, rlon_values, lat2d, lon2d, mapping=ROTATED_MAPPING, dims=("rlat", "rlon")):
    # a 2D lat/lon coordinate pair over the two 1D frame axes, as a CF
    # rotated latitude/longitude grid stores them
    rlat = xr.DataArray(
        np.asarray(rlat_values, dtype=np.float64),
        dims=["rlat"],
        name="rlat",
        attrs={"standard_name": "grid_latitude", "axis": "Y", "units": "degrees"},
    )
    rlon = xr.DataArray(
        np.asarray(rlon_values, dtype=np.float64),
        dims=["rlon"],
        name="rlon",
        attrs={"standard_name": "grid_longitude", "axis": "X", "units": "degrees"},
    )
    lat = xr.DataArray(np.asarray(lat2d), dims=list(dims), name="lat", coords={"rlat": rlat, "rlon": rlon})
    lon = xr.DataArray(np.asarray(lon2d), dims=list(dims), name="lon", coords={"rlat": rlat, "rlon": rlon})
    return UnstructuredGrid(Coordinate(lat), Coordinate(lon), variable_dims=dims, grid_mapping=mapping)


def _rotated_points(rlat_values, rlon_values, projection):
    # the true lat/lon of a regular grid in the rotated frame, as eckit
    # projects them: the reference the synthetic dataset is built on
    rlat_values = np.asarray(rlat_values, dtype=np.float64)
    rlon_values = np.asarray(rlon_values, dtype=np.float64)
    reference = eckit_geo.Grid({
        "grid": [2.0, 2.0],
        "area": [
            float(rlat_values.max()),
            float(rlon_values.min()),
            float(rlat_values.min()),
            float(rlon_values.max()),
        ],
        "order": "i+j+",
        "projection": projection,
    })
    lat, lon = reference.to_latlons()
    shape = (len(rlat_values), len(rlon_values))
    return np.asarray(lat).reshape(shape), np.asarray(lon).reshape(shape)


def test_unstructured_grid_spec_rotated_coord_detection():
    # the frame axes recognised as plain (non-dimension) coordinates of the
    # lat/lon variables, by their axis attribute rather than their
    # standard_name
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    rlat = xr.DataArray(rlat_values, dims=["y"], name="y", attrs={"axis": "Y", "units": "degrees"})
    rlon = xr.DataArray(rlon_values, dims=["x"], name="x", attrs={"axis": "X", "units": "degrees"})
    lat = xr.DataArray(lat2d, dims=["y", "x"], name="lat", coords={"y": rlat, "x": rlon})
    lon = xr.DataArray(lon2d, dims=["y", "x"], name="lon", coords={"y": rlat, "x": rlon})
    grid = UnstructuredGrid(Coordinate(lat), Coordinate(lon), variable_dims=("y", "x"), grid_mapping=ROTATED_MAPPING)

    spec = grid.grid_spec
    assert isinstance(spec, dict)
    assert spec["projection"] == {"type": "rotation", "south_pole": [18.0, -39.25]}


def test_unstructured_grid_spec_rotated_unrecognised_axes():
    # the frame axes cannot be recognised (no CF hints on them): no spec
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    rlat = xr.DataArray(rlat_values, dims=["rlat"], name="rlat")
    rlon = xr.DataArray(rlon_values, dims=["rlon"], name="rlon")
    lat = xr.DataArray(lat2d, dims=["rlat", "rlon"], name="lat", coords={"rlat": rlat, "rlon": rlon})
    lon = xr.DataArray(lon2d, dims=["rlat", "rlon"], name="lon", coords={"rlat": rlat, "rlon": rlon})
    grid = UnstructuredGrid(
        Coordinate(lat), Coordinate(lon), variable_dims=("rlat", "rlon"), grid_mapping=ROTATED_MAPPING
    )

    assert grid.grid_spec is None


def test_rotation_projection():
    from earthkit.geo.grids._regrid.data.xarray.grid_spec import rotation_projection

    projection = rotation_projection(ROTATED_MAPPING)
    # eckit describes the rotation by its south pole, [longitude, latitude]
    assert projection == {"type": "rotation", "south_pole": [18.0, -39.25]}

    # absent means the CF default of 0
    mapping = dict(ROTATED_MAPPING)
    del mapping["north_pole_grid_longitude"]
    assert rotation_projection(mapping) == projection

    assert rotation_projection({}) is None
    assert rotation_projection(None) is None
    assert rotation_projection({"grid_mapping_name": "lambert_azimuthal_equal_area"}) is None

    incomplete = dict(ROTATED_MAPPING)
    del incomplete["grid_north_pole_latitude"]
    assert rotation_projection(incomplete) is None

    # an angle of the grid about its own pole cannot be carried by eckit
    angled = dict(ROTATED_MAPPING, north_pole_grid_longitude=10.0)
    assert rotation_projection(angled) is None


def test_unstructured_grid_spec_rotated():
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["grid"] == pytest.approx([2.0, 2.0])
    assert spec["area"] == pytest.approx([4.0, -4.0, -4.0, 8.0])
    assert spec["order"] == "i+j+"
    assert spec["projection"] == {"type": "rotation", "south_pole": [18.0, -39.25]}

    # the spec builds a grid of the data's shape reproducing its coordinates
    ek_grid = eckit_geo.Grid(spec)
    assert ek_grid.type == "regular_ll"
    assert ek_grid.shape == (5, 7)

    ek_lat, ek_lon = ek_grid.to_latlons()
    data_lat, data_lon = grid.latlons
    assert np.max(np.abs(np.asarray(ek_lat) - data_lat)) <= 1e-3
    dlon = (np.asarray(ek_lon) - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


def test_unstructured_grid_spec_rotated_descending_rlat():
    # a rotated grid with its rlat axis descending: the default, north to
    # south latitude order carries it, no "order" in the spec
    rlat_values = np.arange(4.0, -4.5, -2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    projection = {"type": "rotation", "south_pole": [18.0, -39.25]}

    # the file's points, in its own (descending rlat) row-major order
    lat2d_ascending, lon2d_ascending = _rotated_points(rlat_values[::-1], rlon_values, projection)
    lat2d = lat2d_ascending[::-1]
    lon2d = lon2d_ascending[::-1]

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert "order" not in spec
    assert spec["projection"] == projection

    ek_grid = eckit_geo.Grid(spec)
    assert ek_grid.shape == (5, 7)
    ek_lat, ek_lon = ek_grid.to_latlons()
    data_lat, data_lon = grid.latlons
    assert np.max(np.abs(np.asarray(ek_lat) - data_lat)) <= 1e-3
    dlon = (np.asarray(ek_lon) - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


def test_unstructured_grid_spec_rotated_irregular_rlat():
    # a non-regular rlat axis: no rotated gridspec, so the grid falls back
    # to being handled as unstructured
    rlat_values = np.array([-4.0, -2.0, 0.0, 2.5, 4.0])
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(
        np.arange(-4.0, 4.5, 2.0), rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]}
    )

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_rotated_metadata_mismatch():
    # the grid mapping disagrees with the dataset's own 2D lat/lon (the pole
    # is moved): the spec is rejected and the grid falls back to unstructured
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    wrong_pole = dict(ROTATED_MAPPING, grid_north_pole_latitude=55.0)
    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d, mapping=wrong_pole)
    assert grid.grid_spec is None


@pytest.mark.parametrize("coordinate", ["lat", "lon"])
def test_unstructured_grid_spec_rotated_nan_coordinate(coordinate):
    # one of the dataset's own 2D coordinates is a NaN: it cannot be compared to
    # the points of the spec, so the spec is rejected and the grid falls back to
    # unstructured, like every dataset whose coordinates are not all finite
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    if coordinate == "lat":
        lat2d[1, 2] = np.nan
    else:
        lon2d[1, 2] = np.nan

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_rotated_nan_hides_a_shift():
    # the same rotated grid with its latitudes moved 5 degrees away from the pole
    # the grid mapping states, and one of them replaced by a NaN: the shift alone
    # is rejected by the point check, and the NaN must not hide it -- the largest
    # difference of an array holding a NaN is a NaN, which is not greater than the
    # tolerance, so the comparison has to be an affirmative one over finite points
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    lat2d = lat2d + 5.0
    lat2d[1, 2] = np.nan

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_rotated_nan_free_shift_is_rejected():
    # the check without the NaN: the same shift alone is caught, which is what the
    # test above needs to hold the fix to
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    assert _rotated_grid(rlat_values, rlon_values, lat2d + 5.0, lon2d).grid_spec is None
    assert _rotated_grid(rlat_values, rlon_values, lat2d, lon2d).grid_spec is not None


def test_unstructured_grid_spec_rotated_without_mapping():
    # the same points, but no grid mapping to tell us they are rotated: the
    # behaviour is unchanged and the grid has no gridspec
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, {"type": "rotation", "south_pole": [18.0, -39.25]})

    grid = _rotated_grid(rlat_values, rlon_values, lat2d, lon2d, mapping=None)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_rotated_lon_major():
    # a rotated grid whose coordinates vary the rotated latitude fastest: the
    # frame axes keep their eckit slots (the rotation projects the latitude
    # slot as the rotated latitude), only the iteration order is swapped so
    # that the spec's points come out in the same order as the latlons
    rlat_values = np.arange(-4.0, 4.5, 2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    projection = {"type": "rotation", "south_pole": [18.0, -39.25]}
    lat2d, lon2d = _rotated_points(rlat_values, rlon_values, projection)

    rlat = xr.DataArray(
        rlat_values,
        dims=["rlat"],
        name="rlat",
        attrs={"standard_name": "grid_latitude", "axis": "Y", "units": "degrees"},
    )
    rlon = xr.DataArray(
        rlon_values,
        dims=["rlon"],
        name="rlon",
        attrs={"standard_name": "grid_longitude", "axis": "X", "units": "degrees"},
    )
    lat = xr.DataArray(lat2d.T, dims=["rlon", "rlat"], name="lat", coords={"rlat": rlat, "rlon": rlon})
    lon = xr.DataArray(lon2d.T, dims=["rlon", "rlat"], name="lon", coords={"rlat": rlat, "rlon": rlon})
    grid = UnstructuredGrid(
        Coordinate(lat), Coordinate(lon), variable_dims=("rlon", "rlat"), grid_mapping=ROTATED_MAPPING
    )
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["order"] == "j+i+"
    assert spec["projection"] == projection

    # the eckit grid is shaped by its axes, not by the data indexing
    ek_grid = eckit_geo.Grid(spec)
    assert ek_grid.shape == (5, 7)

    ek_lat, ek_lon = ek_grid.to_latlons()
    data_lat, data_lon = grid.latlons
    assert np.max(np.abs(np.asarray(ek_lat) - data_lat)) <= 1e-3
    dlon = (np.asarray(ek_lon) - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


def test_unstructured_grid_spec_rotated_lon_major_descending_rlat():
    # the same, with a descending rlat axis: the swapped order carries the
    # north to south iteration of the default order
    rlat_values = np.arange(4.0, -4.5, -2.0)
    rlon_values = np.arange(-4.0, 9.0, 2.0)
    projection = {"type": "rotation", "south_pole": [18.0, -39.25]}

    lat2d_ascending, lon2d_ascending = _rotated_points(rlat_values[::-1], rlon_values, projection)
    lat2d, lon2d = lat2d_ascending[::-1], lon2d_ascending[::-1]

    rlat = xr.DataArray(
        rlat_values,
        dims=["rlat"],
        name="rlat",
        attrs={"standard_name": "grid_latitude", "axis": "Y", "units": "degrees"},
    )
    rlon = xr.DataArray(
        rlon_values,
        dims=["rlon"],
        name="rlon",
        attrs={"standard_name": "grid_longitude", "axis": "X", "units": "degrees"},
    )
    lat = xr.DataArray(lat2d.T, dims=["rlon", "rlat"], name="lat", coords={"rlat": rlat, "rlon": rlon})
    lon = xr.DataArray(lon2d.T, dims=["rlon", "rlat"], name="lon", coords={"rlat": rlat, "rlon": rlon})
    grid = UnstructuredGrid(
        Coordinate(lat), Coordinate(lon), variable_dims=("rlon", "rlat"), grid_mapping=ROTATED_MAPPING
    )
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["order"] == "j-i+"
    assert spec["projection"] == projection

    ek_grid = eckit_geo.Grid(spec)
    assert ek_grid.shape == (5, 7)

    ek_lat, ek_lon = ek_grid.to_latlons()
    data_lat, data_lon = grid.latlons
    assert np.max(np.abs(np.asarray(ek_lat) - data_lat)) <= 1e-3
    dlon = (np.asarray(ek_lon) - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


# a Lambert azimuthal equal-area projection on GRS80, with the false easting and
# northing of the efas sample: a CF grid mapping describes it by its PROJ string
LAEA_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

PROJECTED_MAPPING = {
    "grid_mapping_name": "lambert_azimuthal_equal_area",
    "proj4_params": LAEA_PROJ4,
}


def _projected_points(x_values, y_values, proj4_params=LAEA_PROJ4):
    # the true lat/lon of the points of the plane, as PROJ computes them: the
    # reference a synthetic projected dataset is built on
    from pyproj import CRS, Transformer

    transformer = Transformer.from_crs(CRS.from_proj4(proj4_params), "EPSG:4326", always_xy=True)
    xv, yv = np.meshgrid(np.asarray(x_values, dtype=np.float64), np.asarray(y_values, dtype=np.float64))
    lon, lat = transformer.transform(xv, yv)
    return lat, lon


def _projected_grid(x_values, y_values, lat2d, lon2d, mapping=PROJECTED_MAPPING, dims=("y", "x")):
    # a 2D lat/lon coordinate pair over the two 1D axes of the projection plane,
    # as a CF projected grid stores them
    x = xr.DataArray(
        np.asarray(x_values),
        dims=["x"],
        name="x",
        attrs={"standard_name": "projection_x_coordinate", "axis": "X", "units": "m"},
    )
    y = xr.DataArray(
        np.asarray(y_values),
        dims=["y"],
        name="y",
        attrs={"standard_name": "projection_y_coordinate", "axis": "Y", "units": "m"},
    )
    lat = xr.DataArray(np.asarray(lat2d), dims=list(dims), name="latitude", coords={"x": x, "y": y})
    lon = xr.DataArray(np.asarray(lon2d), dims=list(dims), name="longitude", coords={"x": x, "y": y})
    return UnstructuredGrid(Coordinate(lat), Coordinate(lon), variable_dims=dims, grid_mapping=mapping)


# the axes of the synthetic projected grid: 5 km, y descending as in the efas sample
X_VALUES = 4321000.0 + 5000.0 * np.arange(6)
Y_VALUES = 3210000.0 - 5000.0 * np.arange(4)


@pytest.mark.parametrize("ascending_y", [False, True])
def test_unstructured_grid_spec_projected(ascending_y):
    y_values = Y_VALUES[::-1] if ascending_y else Y_VALUES
    lat2d, lon2d = _projected_points(X_VALUES, y_values)

    grid = _projected_grid(X_VALUES, y_values, lat2d, lon2d)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["type"] == "regular_xy"
    # the ranges are [start, stop, step], in the direction each axis runs
    assert spec["x"] == pytest.approx([4321000.0, 4346000.0, 5000.0])
    assert spec["y"] == pytest.approx([float(y_values[0]), float(y_values[-1]), -5000.0 if not ascending_y else 5000.0])
    assert spec["projection"] == {"type": "proj", "source": LAEA_PROJ4}

    # the spec builds a grid of the data's shape reproducing its coordinates
    ek_grid = eckit_geo.Grid(spec)
    assert ek_grid.type == "regular_xy"
    assert ek_grid.shape == (4, 6)

    ek_lat, ek_lon = ek_grid.to_latlons()
    data_lat, data_lon = grid.latlons
    assert np.max(np.abs(np.asarray(ek_lat) - data_lat)) <= 1e-3
    dlon = (np.asarray(ek_lon) - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


def test_unstructured_grid_spec_projected_float32_axes():
    # the axes are the ones of the efas sample, integers below 2**23, so float32
    # carries them exactly: the grid is still detected, with the same spec
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    grid = _projected_grid(X_VALUES.astype(np.float32), Y_VALUES.astype(np.float32), lat2d, lon2d)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["x"] == pytest.approx([4321000.0, 4346000.0, 5000.0])
    assert spec["y"] == pytest.approx([3210000.0, 3195000.0, -5000.0])
    assert eckit_geo.Grid(spec).shape == (4, 6)


def test_unstructured_grid_spec_projected_irregular_x():
    # the x axis is not regular: no projected gridspec, so the grid falls back to
    # being handled as unstructured. Its coordinates are the ones of the axis as
    # the file holds it, so it is the axis alone that is rejected here
    x_values = X_VALUES.copy()
    x_values[3] += 1000.0
    lat2d, lon2d = _projected_points(x_values, Y_VALUES)

    grid = _projected_grid(x_values, Y_VALUES, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_projected_metadata_mismatch():
    # the axes are regular, but the projection the grid mapping describes is not
    # the one the dataset's own 2D lat/lon were computed with: the spec is
    # rejected by the point check and the grid falls back to unstructured
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    wrong = dict(PROJECTED_MAPPING, proj4_params=LAEA_PROJ4.replace("+lat_0=52", "+lat_0=40"))
    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d, mapping=wrong)
    assert grid.grid_spec is None


@pytest.mark.parametrize("coordinate", ["lat", "lon"])
def test_unstructured_grid_spec_projected_nan_coordinate(coordinate):
    # one of the dataset's own 2D coordinates is a NaN: it cannot be compared to
    # the points of the spec, so the spec is rejected and the grid falls back to
    # unstructured, like every dataset whose coordinates are not all finite
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    if coordinate == "lat":
        lat2d[1, 2] = np.nan
    else:
        lon2d[1, 2] = np.nan

    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_projected_nan_hides_a_shift():
    # the same projected grid with its latitudes moved 5 degrees away from the
    # projection the grid mapping states, and one of them replaced by a NaN: the
    # shift alone is rejected by the point check, and the NaN must not hide it --
    # the largest difference of an array holding a NaN is a NaN, which is not
    # greater than the tolerance, so the comparison has to be an affirmative one
    # over finite points
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    lat2d = lat2d + 5.0
    lat2d[1, 2] = np.nan

    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_projected_nan_free_shift_is_rejected():
    # the check without the NaN: the same shift alone is caught, which is what the
    # test above needs to hold the fix to
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    assert _projected_grid(X_VALUES, Y_VALUES, lat2d + 5.0, lon2d).grid_spec is None
    assert _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d).grid_spec is not None


def test_grid_spec_matches_rejects_non_finite_points():
    # the point check itself, on a plain regular lat/lon grid: a non-finite
    # coordinate of either side fails it, it never passes on a NaN
    from earthkit.geo.grids._regrid.data.xarray.grid import _grid_spec_matches

    spec = {"grid": [10.0, 10.0], "area": [50.0, 0.0, -50.0, 50.0]}
    grid = eckit_geo.Grid(spec)
    shape = tuple(grid.shape)
    lat, lon = grid.to_latlons()
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)

    assert _grid_spec_matches(spec, shape, lat, lon)

    shifted = lat + 5.0
    assert not _grid_spec_matches(spec, shape, shifted, lon)

    with_nan = shifted.copy()
    with_nan[3] = np.nan
    assert not _grid_spec_matches(spec, shape, with_nan, lon)

    with_nan_lon = lon.copy()
    with_nan_lon[0] = np.nan
    assert not _grid_spec_matches(spec, shape, lat, with_nan_lon)

    # a grid of a different shape than the data is rejected, finite or not
    assert not _grid_spec_matches(spec, (shape[0] + 1, shape[1]), lat, lon)


def test_unstructured_grid_spec_projected_crs_wkt():
    # a mapping that describes its projection by CF crs_wkt rather than by
    # proj4_params: the PROJ definition is translated with pyproj and then
    # checked against the dataset's coordinates, as any other
    from pyproj import CRS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wkt = CRS.from_proj4(LAEA_PROJ4).to_wkt()

    mapping = {"grid_mapping_name": "lambert_azimuthal_equal_area", "crs_wkt": wkt}
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d, mapping=mapping)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["projection"]["type"] == "proj"
    assert eckit_geo.Grid(spec).shape == (4, 6)


def test_unstructured_grid_spec_projected_unrecognised_axes():
    # the plane axes cannot be recognised (no CF hints on them): no spec
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    x = xr.DataArray(X_VALUES, dims=["x"], name="x")
    y = xr.DataArray(Y_VALUES, dims=["y"], name="y")
    lat = xr.DataArray(lat2d, dims=["y", "x"], name="latitude", coords={"x": x, "y": y})
    lon = xr.DataArray(lon2d, dims=["y", "x"], name="longitude", coords={"x": x, "y": y})
    grid = UnstructuredGrid(Coordinate(lat), Coordinate(lon), variable_dims=("y", "x"), grid_mapping=PROJECTED_MAPPING)

    assert grid.grid_spec is None


def test_unstructured_grid_spec_projected_without_mapping():
    # the same points, but no grid mapping to say what projection they are in:
    # the behaviour is unchanged and the grid has no gridspec
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d, lon2d, mapping=None)
    assert grid.grid_spec is None


def test_unstructured_grid_spec_projected_y_major():
    # the data varies fastest along y: eckit scans a regular_xy grid with its x
    # axis fastest and cannot be told otherwise, so the point check rejects the
    # spec and the grid falls back to unstructured
    lat2d, lon2d = _projected_points(X_VALUES, Y_VALUES)

    grid = _projected_grid(X_VALUES, Y_VALUES, lat2d.T, lon2d.T, dims=("x", "y"))
    assert grid.grid_spec is None


# A projected grid held as the two 1D axes of its plane, with no 2D latitude and
# longitude to check a gridspec against, as in the polar stereographic AMSR2
# sample. pyproj reads the projection of a CF grid mapping from its crs_wkt or
# from the individual parameters, and never from proj4_params, so the mapping
# here carries both: the parameters for the grid's own points, the PROJ string
# for the gridspec. WGS84 rather than the GRS80 of LAEA_PROJ4, so that the two
# descriptions of it are the same projection.
LAEA_WGS84 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +a=6378137 +rf=298.257223563 +units=m +no_defs"

LAEA_CF = {
    "grid_mapping_name": "lambert_azimuthal_equal_area",
    "latitude_of_projection_origin": 52.0,
    "longitude_of_projection_origin": 10.0,
    "false_easting": 4321000.0,
    "false_northing": 3210000.0,
    "semi_major_axis": 6378137.0,
    "inverse_flattening": 298.257223563,
}

MESH_PROJECTED_MAPPING = dict(LAEA_CF, proj4_params=LAEA_WGS84)


def _mesh_projected_grid(
    x_values, y_values, mapping=MESH_PROJECTED_MAPPING, variable_dims=("y", "x"), dim_names=("x", "y")
):
    # a projected grid over its two 1D axes only, as a file which stores neither
    # a 2D lat/lon field nor a rotated frame holds it; dim_names lets a test hold
    # the axes on dimensions named differently from the coordinates carrying them
    x = xr.DataArray(
        np.asarray(x_values),
        dims=[dim_names[0]],
        name="x",
        attrs={"standard_name": "projection_x_coordinate", "axis": "X", "units": "m"},
    )
    y = xr.DataArray(
        np.asarray(y_values),
        dims=[dim_names[1]],
        name="y",
        attrs={"standard_name": "projection_y_coordinate", "axis": "Y", "units": "m"},
    )
    return MeshProjectionGrid(Coordinate(x), Coordinate(y), mapping, variable_dims=variable_dims)


def _assert_mesh_points_match(grid, spec):
    # the eckit grid of the spec holds the same points as the grid's own, in the
    # same order, the check a detected spec has to pass
    ek_grid = eckit_geo.Grid(spec)
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in ek_grid.to_latlons())
    data_lat, data_lon = grid.latlons

    assert np.max(np.abs(ek_lat - data_lat)) <= 1e-3
    dlon = (ek_lon - data_lon + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(dlon)) <= 1e-3


@pytest.mark.parametrize("ascending_y", [False, True])
def test_mesh_projection_grid_spec_projected(ascending_y):
    # the axes of the projection plane and a PROJ description of the projection
    # are enough: the grid is a regular_xy one in that plane, whichever way its
    # y axis runs, because its ranges are taken in the direction of the axis
    y_values = Y_VALUES[::-1] if ascending_y else Y_VALUES

    grid = _mesh_projected_grid(X_VALUES, y_values)
    spec = grid.grid_spec

    assert grid.variable_dims == ("y", "x")
    assert grid.grid_type == "regular_xy"

    assert isinstance(spec, dict)
    assert spec["type"] == "regular_xy"
    # the ranges are [start, stop, step], in the direction each axis runs
    assert spec["x"] == pytest.approx([4321000.0, 4346000.0, 5000.0])
    assert spec["y"] == pytest.approx([float(y_values[0]), float(y_values[-1]), 5000.0 if ascending_y else -5000.0])
    assert spec["projection"] == {"type": "proj", "source": LAEA_WGS84}

    assert eckit_geo.Grid(spec).shape == (4, 6)
    _assert_mesh_points_match(grid, spec)

    # the points are the mesh of the axes, x fastest, as the (y, x) data holds them
    lat, lon = grid.latlons
    assert lat.shape == (24,)
    ref_lat, ref_lon = _projected_points(X_VALUES, y_values, proj4_params=LAEA_WGS84)
    assert np.max(np.abs(lat - ref_lat.reshape(-1))) <= 1e-3
    assert np.max(np.abs(lon - ref_lon.reshape(-1))) <= 1e-3


def test_mesh_projection_grid_spec_projected_crs_wkt():
    # the polar stereographic sample's case: a mapping that carries its PROJ
    # description as crs_wkt, and no proj4_params at all
    from pyproj import CRS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wkt = CRS.from_proj4(LAEA_WGS84).to_wkt()

    mapping = dict(LAEA_CF, crs_wkt=wkt)
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, mapping=mapping)
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert spec["type"] == "regular_xy"
    assert spec["projection"]["type"] == "proj"
    assert eckit_geo.Grid(spec).shape == (4, 6)
    _assert_mesh_points_match(grid, spec)


def test_mesh_projection_grid_spec_projected_irregular_x():
    # the axes are not regular: no gridspec, but the grid still describes its
    # points, which are then taken as unstructured ones
    x_values = X_VALUES.copy()
    x_values[3] += 1000.0

    grid = _mesh_projected_grid(x_values, Y_VALUES)

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"
    assert grid.latlons[0].shape == (24,)


def test_mesh_projection_grid_spec_projected_too_few_points():
    # an axis with a single point is not a grid either
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES[:1])

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"
    assert grid.latlons[0].shape == (6,)


def test_mesh_projection_grid_spec_projected_metadata_mismatch():
    # the axes are regular, but the PROJ string and the CF parameters the points
    # are computed from describe different projections: the spec is rejected by
    # the point check, and the grid keeps its unstructured treatment
    mapping = dict(MESH_PROJECTED_MAPPING, proj4_params=LAEA_WGS84.replace("+lat_0=52", "+lat_0=40"))
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, mapping=mapping)

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"
    assert grid.latlons[0].shape == (24,)


def test_mesh_projection_grid_spec_projected_without_proj_string():
    # a mapping holding its projection in the CF parameters only, with neither
    # proj4_params nor crs_wkt: the points are computable, but there is nothing
    # a gridspec can be built from, so the grid falls back to unstructured
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, mapping=LAEA_CF)

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"

    lat, lon = grid.latlons
    assert np.isfinite(lat).all() and np.isfinite(lon).all()
    ref_lat, ref_lon = _projected_points(X_VALUES, Y_VALUES, proj4_params=LAEA_WGS84)
    assert np.max(np.abs(lat - ref_lat.reshape(-1))) <= 1e-3


def test_mesh_projection_grid_spec_projected_y_major():
    # the data varies fastest along y: the points follow the data, eckit scans a
    # regular_xy grid with its x axis fastest and cannot be told otherwise, so the
    # point check rejects the spec and the grid falls back to unstructured
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, variable_dims=("x", "y"))

    assert grid.variable_dims == ("x", "y")
    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"

    # the points are the mesh of the axes with y fastest, in the order the data
    # is flattened into them, not the x-fastest one
    lat, lon = grid.latlons
    ref_lat, ref_lon = (
        a.reshape(4, 6).T.reshape(-1) for a in _projected_points(X_VALUES, Y_VALUES, proj4_params=LAEA_WGS84)
    )
    assert np.max(np.abs(lat - ref_lat)) <= 1e-3
    assert np.max(np.abs(lon - ref_lon)) <= 1e-3

    x, y = grid.xys
    assert np.max(np.abs(x - np.repeat(X_VALUES, 4))) == 0.0
    assert np.max(np.abs(y - np.tile(Y_VALUES, 6))) == 0.0


def test_mesh_projection_grid_spec_projected_renamed_dims():
    # the axes sit on dimensions whose names differ from the coordinates
    # carrying them, as in a dataset renamed with xarray.Dataset.rename_dims:
    # the point order is read off the axes' dimensions, never their names
    # a variable indexed (j, i) -- (y, x) by another name -- holds the mesh x
    # fastest and keeps its regular_xy spec
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, variable_dims=("j", "i"), dim_names=("i", "j"))
    spec = grid.grid_spec

    assert isinstance(spec, dict)
    assert grid.grid_type == "regular_xy"
    _assert_mesh_points_match(grid, spec)

    x, y = grid.xys
    assert np.max(np.abs(x - np.tile(X_VALUES, 4))) == 0.0
    assert np.max(np.abs(y - np.repeat(Y_VALUES, 6))) == 0.0

    # a variable indexed (i, j) varies fastest along y: the spec is not claimed
    # and the fallback points follow the data's own order
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, variable_dims=("i", "j"), dim_names=("i", "j"))

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"

    x, y = grid.xys
    assert np.max(np.abs(x - np.repeat(X_VALUES, 4))) == 0.0
    assert np.max(np.abs(y - np.tile(Y_VALUES, 6))) == 0.0


def test_mesh_projection_grid_spec_projected_unknown_dim_order():
    # an indexing which matches neither order of the axes' dimensions: no
    # regular_xy spec is claimed on an assumed order, and the points keep the
    # (y, x) layout a projected grid is usually stored in
    grid = _mesh_projected_grid(X_VALUES, Y_VALUES, variable_dims=("p", "q"), dim_names=("i", "j"))

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"

    lat, lon = grid.latlons
    assert lat.shape == (24,)
    assert np.isfinite(lat).all() and np.isfinite(lon).all()


# the points of a projected grid held as one dimension of pairs: x and y are two
# coordinates over that same dimension, so the grid is a scatter over the plane
# and not the mesh of two axes. The pairs are 6, where the mesh of the axes would
# be 36.
SCATTER_X = 4321000.0 + 5000.0 * np.array([0.0, 1.0, 2.0, 0.0, 4.0, 6.0])
SCATTER_Y = 3210000.0 - 5000.0 * np.array([0.0, 1.0, 2.0, 6.0, 3.0, 5.0])


def _unstructured_projected_grid(
    x_values, y_values, mapping=MESH_PROJECTED_MAPPING, variable_dims=None, dim_names=("values",)
):
    # a projected grid over the dimension(s) its coordinates share: one point per
    # index, the case the guesser hands to UnstructuredProjectionGrid. dim_names
    # gives the dimensions the coordinates are defined on, which the variable may
    # index in another order.
    dims = list(dim_names)
    x = xr.DataArray(
        np.asarray(x_values),
        dims=dims,
        name="x",
        attrs={"standard_name": "projection_x_coordinate", "axis": "X", "units": "m"},
    )
    y = xr.DataArray(
        np.asarray(y_values),
        dims=dims,
        name="y",
        attrs={"standard_name": "projection_y_coordinate", "axis": "Y", "units": "m"},
    )
    return UnstructuredProjectionGrid(Coordinate(x), Coordinate(y), mapping, variable_dims=variable_dims)


def _pair_latlons(x_values, y_values, proj4_params=LAEA_WGS84):
    # the (lat, lon) pyproj computes for a flat list of (x, y) pairs
    from pyproj import CRS, Transformer

    transformer = Transformer.from_crs(CRS.from_proj4(proj4_params), "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(np.asarray(x_values, dtype=np.float64), np.asarray(y_values, dtype=np.float64))
    return lat, lon


def test_unstructured_projection_grid_points():
    # the x/y coordinates are paired along one dimension: the grid is the scatter
    # of those pairs, with no gridspec to claim for it
    grid = _unstructured_projected_grid(SCATTER_X, SCATTER_Y, variable_dims=("values",))

    assert grid.variable_dims == ("values",)
    assert grid.grid_type == "unstructured"
    assert grid.grid_spec is None

    # the points are the pairs, not the mesh of the two coordinates
    x, y = grid.xys
    assert x.shape == (6,)
    assert np.max(np.abs(x - SCATTER_X)) == 0.0
    assert np.max(np.abs(y - SCATTER_Y)) == 0.0

    lat, lon = grid.latlons
    ref_lat, ref_lon = _pair_latlons(SCATTER_X, SCATTER_Y)
    assert lat.shape == (6,)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-6
    assert np.max(np.abs((lon - ref_lon + 180.0) % 360.0 - 180.0)) <= 1e-6


def test_unstructured_projection_grid_points_default_dims():
    # a grid built without the variable's indexing takes the order its own
    # coordinates are defined on
    grid = _unstructured_projected_grid(SCATTER_X, SCATTER_Y)

    assert grid.variable_dims == ("values",)

    lat, lon = grid.latlons
    ref_lat, ref_lon = _pair_latlons(SCATTER_X, SCATTER_Y)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-6


def test_unstructured_projection_grid_points_2d(caplog):
    # the coordinates are two dimensional over the two dimensions the variable
    # shares with them: the points are still the pairs of the coordinates, in the
    # order the data is flattened into them
    x_values = np.arange(6, dtype=np.float64).reshape(2, 3) * 5000.0 + 4321000.0
    y_values = 3210000.0 - np.arange(6, dtype=np.float64).reshape(2, 3) * 5000.0

    grid = _unstructured_projected_grid(x_values, y_values, variable_dims=("i", "j"), dim_names=("i", "j"))
    assert grid.variable_dims == ("i", "j")
    assert grid.grid_spec is None

    x, y = grid.xys
    assert x.shape == (6,)
    assert np.max(np.abs(x - x_values.flatten())) == 0.0
    assert np.max(np.abs(y - y_values.flatten())) == 0.0

    # the same indexing spelled as a list of dimensions is still the order the
    # coordinates are defined on, not a mis-reading of it
    listed = _unstructured_projected_grid(x_values, y_values, variable_dims=["i", "j"], dim_names=("i", "j"))

    x, y = listed.xys
    assert np.max(np.abs(x - x_values.flatten())) == 0.0

    # the same variable indexed the other way round: the points follow the data,
    # which the grid says out loud rather than handing out a mis-ordered scatter
    transposed = _unstructured_projected_grid(x_values, y_values, variable_dims=("j", "i"), dim_names=("i", "j"))

    with caplog.at_level(logging.WARNING, logger="earthkit.geo.grids._regrid.data.xarray.grid"):
        x, y = transposed.xys

    assert "does not match grid indexing" in caplog.text
    assert x.shape == (6,)
    assert np.max(np.abs(x - x_values.transpose().flatten())) == 0.0
    assert np.max(np.abs(y - y_values.transpose().flatten())) == 0.0

    lat, lon = transposed.latlons
    ref_lat, ref_lon = _pair_latlons(x_values.transpose().flatten(), y_values.transpose().flatten())
    assert np.max(np.abs(lat - ref_lat)) <= 1e-6
    assert np.max(np.abs(lon - ref_lon)) <= 1e-6


def test_unstructured_projection_grid_points_no_proj_string():
    # the points come from the CF parameters of the mapping even when it holds no
    # PROJ description; the grid keeps its unstructured treatment either way
    grid = _unstructured_projected_grid(SCATTER_X, SCATTER_Y, mapping=LAEA_CF, variable_dims=("values",))

    assert grid.grid_spec is None
    assert grid.grid_type == "unstructured"

    lat, lon = grid.latlons
    ref_lat, ref_lon = _pair_latlons(SCATTER_X, SCATTER_Y)
    assert lat.shape == (6,)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-6
    assert np.max(np.abs(lon - ref_lon)) <= 1e-6


def test_unstructured_projection_grid_points_lonlat_projection():
    # a projection given as the string "epsg:4326", which the guesser never builds
    # (it always passes the attributes of a grid mapping variable, a dict): the
    # axes are already degrees, and xys has always returned them the other way
    # round for that projection. latlons still pairs them the way the plane does.
    grid = _unstructured_projected_grid(SCATTER_X, SCATTER_Y, mapping="epsg:4326", variable_dims=("values",))

    x, y = grid.xys
    assert np.max(np.abs(x - SCATTER_Y)) == 0.0
    assert np.max(np.abs(y - SCATTER_X)) == 0.0

    lat, lon = grid.latlons
    assert np.max(np.abs(lat - SCATTER_Y)) <= 1e-9
    assert np.max(np.abs(lon - SCATTER_X)) <= 1e-9


def test_unstructured_projection_grid_rejects_unpaired_coordinates():
    # coordinates which are not paired index by index are not this grid's case
    with pytest.raises(AssertionError):
        _unstructured_projected_grid(SCATTER_X, SCATTER_Y[:4], variable_dims=("values",))

    with pytest.raises(AssertionError):
        # a variable indexed by one dimension cannot hold a grid defined on two
        _unstructured_projected_grid(
            np.arange(6.0).reshape(2, 3),
            np.arange(6.0).reshape(2, 3),
            variable_dims=("values",),
            dim_names=("i", "j"),
        )


def test_proj_projection():
    from earthkit.geo.grids._regrid.data.xarray.grid_spec import proj_projection

    assert proj_projection(PROJECTED_MAPPING) == {"type": "proj", "source": LAEA_PROJ4}

    # a netCDF attribute may come as bytes
    assert proj_projection(dict(PROJECTED_MAPPING, proj4_params=LAEA_PROJ4.encode())) == proj_projection(
        PROJECTED_MAPPING
    )

    assert proj_projection({}) is None
    assert proj_projection(None) is None
    assert proj_projection("not a mapping") is None
    assert proj_projection({"grid_mapping_name": "lambert_azimuthal_equal_area"}) is None
    assert proj_projection(dict(PROJECTED_MAPPING, proj4_params="   ")) is None

    # a rotated grid mapping is left to rotation_projection, even with a PROJ string
    rotated = dict(ROTATED_MAPPING, proj4_params=LAEA_PROJ4)
    assert proj_projection(rotated) is None


def test_regular_xy_grid_spec():
    from earthkit.geo.grids._regrid.data.xarray.grid_spec import regular_xy_grid_spec

    projection = {"type": "proj", "source": LAEA_PROJ4}

    spec = regular_xy_grid_spec(X_VALUES, Y_VALUES, projection)
    assert spec == {
        "type": "regular_xy",
        "x": [4321000.0, 4346000.0, 5000.0],
        "y": [3210000.0, 3195000.0, -5000.0],
        "projection": projection,
    }

    # an irregular axis, or an axis with a single point, is not a grid
    assert regular_xy_grid_spec(np.array([0.0, 5000.0, 12000.0]), Y_VALUES, projection) is None
    assert regular_xy_grid_spec(X_VALUES, np.array([3210000.0]), projection) is None
