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
