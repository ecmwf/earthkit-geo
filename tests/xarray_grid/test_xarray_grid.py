# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest

from earthkit.geo.grids._regrid.data.xarray.loader import variables as get_variables
from earthkit.geo.utils.testing import NO_COVJSON, covjson_to_xarray, get_test_data

xr = pytest.importorskip("xarray")

GRID_SPEC = {"grid": [10, 10], "area": [50, 0, 30, 20]}


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_xarray_grid_2d_1(lat_name, lon_name):
    # Dimensions:  (level: 2, lat: 3, lon: 3)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #   * lat      (lat) int64 24B 50 40 30
    #   * lon      (lon) int64 24B 0 10 20
    # Data variables:
    #     a        (level, lat, lon) int64 144B 11 12 13 21 22 23 ... 25 26 34 35 36

    dims = {"level": 2, lat_name: 3, lon_name: 3}
    coords = {
        "level": np.array([700, 500]),
        lat_name: np.array([50, 40, 30]),
        lon_name: np.array([0, 10, 20]),
    }

    lat_ref = np.array([[50, 50, 50], [40, 40, 40], [30, 30, 30]])
    lon_ref = np.array([[0, 10, 20], [0, 10, 20], [0, 10, 20]])

    data = np.array([
        [[11, 12, 13], [21, 22, 23], [31, 32, 33]],
        [[14, 15, 16], [24, 25, 26], [34, 35, 36]],
    ])

    a = xr.Variable(
        dims,
        data,
    )
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    grid_spec = {"grid": [10, 10], "area": [50, 0, 30, 20]}
    variable_grid_dims = (lat_name, lon_name)

    variables = get_variables(ds_in, user_ek_grid=grid_spec)

    v = variables[0]

    assert v.name == "a"
    assert v.ek_grid.spec == grid_spec

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert np.allclose(lat, lat_ref.flatten())
    assert np.allclose(lon, lon_ref.flatten())


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_xarray_grid_2d_2(lat_name, lon_name):
    # Dimensions:  (level: 2, y: 3, x: 2)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #     lat      (y, x) int64 48B 50 50 40 40 30 30
    #     lon      (y, x) int64 48B 0 10 0 10 0 10
    # Dimensions without coordinates: y, x
    # Data variables:
    #     a        (level, y, x) int64 96B 11 12 21 22 31 32 14 15 24 25 34 35

    dims = {"level": 2, "y": 3, "x": 2}
    coords = {
        "level": np.array([700, 500]),
        lat_name: (["y", "x"], np.array([[50, 50], [40, 40], [30, 30]])),
        lon_name: (["y", "x"], np.array([[0, 10], [0, 10], [0, 10]])),
    }

    lat_ref = np.array([[50, 50], [40, 40], [30, 30]])
    lon_ref = np.array([[0, 10], [0, 10], [0, 10]])

    data = np.array([
        [[11, 12], [21, 22], [31, 32]],
        [[14, 15], [24, 25], [34, 35]],
    ])

    a = xr.Variable(dims, data)
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    variable_grid_dims = ("y", "x")

    variables = get_variables(ds_in, user_ek_grid=GRID_SPEC)

    v = variables[0]

    assert v.name == "a"
    assert v.ek_grid.spec == GRID_SPEC

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (6,)
    assert lon.shape == (6,)
    assert np.allclose(lat, lat_ref.flatten())
    assert np.allclose(lon, lon_ref.flatten())


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_xarray_grid_1d_1(lat_name, lon_name):
    # Dimensions:  (level: 2, values: 9)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #     lat      (values) int64 72B 50 50 50 40 40 40 30 30 30
    #     lon      (values) int64 72B 0 10 20 0 10 20 0 10 20
    # Dimensions without coordinates: values
    # Data variables:
    #     a        (level, values) int64 144B 11 12 13 21 22 23 ... 24 25 26 34 35 36

    dims = {"level": 2, "values": 9}
    coords = {
        "level": np.array([700, 500]),
        lat_name: ("values", np.array([50, 50, 50, 40, 40, 40, 30, 30, 30])),
        lon_name: ("values", np.array([0, 10, 20, 0, 10, 20, 0, 10, 20])),
    }

    data = np.array([
        [11, 12, 13, 21, 22, 23, 31, 32, 33],
        [14, 15, 16, 24, 25, 26, 34, 35, 36],
    ])

    a = xr.Variable(dims, data)
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    variable_grid_dims = ("values",)

    variables = get_variables(ds_in, user_ek_grid=GRID_SPEC)

    v = variables[0]

    assert v.name == "a"
    assert v.ek_grid.spec == GRID_SPEC

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (9,)
    assert lon.shape == (9,)
    assert np.allclose(lat, coords[lat_name][1])
    assert np.allclose(lon, coords[lon_name][1])


def test_xarray_grid_netcdf_ll_1():
    path = get_test_data("test_single.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    variable_grid_dims = ("latitude", "longitude")
    points_num = 7 * 12

    # auto grid-spec detection
    variables = get_variables(ds)

    assert len(variables) == 1
    v = variables[0]
    assert v.name == "t2m"

    ek_grid = v.ek_grid
    assert ek_grid.type == "unstructured_ll"
    assert ek_grid.shape == (points_num,)

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)


def test_xarray_grid_netcdf_ll_2():
    path = get_test_data("test_single.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    variable_grid_dims = ("latitude", "longitude")
    points_num = 7 * 12

    # user defined grid-spec
    grid_spec = {"grid": [30.0, 30.0]}
    variables = get_variables(ds, user_ek_grid=grid_spec)

    assert len(variables) == 1
    v = variables[0]
    assert v.name == "t2m"

    ek_grid = v.ek_grid
    assert ek_grid.spec == grid_spec
    assert ek_grid.shape == (7, 12)

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_xarray_grid_laea():
    path = get_test_data("efas.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    variable_grid_dims = ("y", "x")
    points_num = 950 * 1000

    variables = get_variables(ds)

    assert len(variables) == 3
    v = variables[0]
    assert v.name == "dis06"

    ek_grid = v.ek_grid
    assert ek_grid.type == "unstructured_ll"
    assert ek_grid.shape == (points_num,)

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)
    assert np.allclose(lat, ds["latitude"].values.flatten())
    assert np.allclose(lon, ds["longitude"].values.flatten())


@pytest.mark.download
def test_xarray_grid_cordex():
    path = get_test_data("cordex.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    variable_grid_dims = ("rlat", "rlon")
    points_num = 19 * 15

    variables = get_variables(ds)

    assert len(variables) == 1

    v = variables[0]
    assert v.name == "sfcWind"

    ek_grid = v.ek_grid
    assert ek_grid.type == "unstructured_ll"
    assert ek_grid.shape == (points_num,)

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)
    assert np.allclose(lat, ds["lat"].values.flatten())
    assert np.allclose(lon, ds["lon"].values.flatten())


@pytest.mark.download
@pytest.mark.skipif(NO_COVJSON, reason="No covjsonkit available")
def test_xarray_grid_covjson_points():
    path = get_test_data("points.covjson", subfolder="xr")
    ds = covjson_to_xarray(path)

    variable_grid_dims = ("points",)
    points_num = 2770

    variables = get_variables(ds)

    assert len(variables) == 2

    v = variables[0]
    assert v.name == "10u"

    ek_grid = v.ek_grid
    assert ek_grid.type == "unstructured_ll"
    assert ek_grid.shape == (points_num,)

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)
    assert np.allclose(lat, ds["latitude"].values.flatten())
    assert np.allclose(lon, ds["longitude"].values.flatten())
