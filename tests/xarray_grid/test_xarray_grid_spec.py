# (C) Copyright 2024 Anemoi contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest

from earthkit.geo.grids._regrid.data.xarray.coordinates import Coordinate
from earthkit.geo.grids._regrid.data.xarray.grid import MeshedGrid, UnstructuredGrid
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
