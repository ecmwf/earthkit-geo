# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest

from earthkit.geo import regrid
from earthkit.geo.utils.testing import (
    NO_COVJSON,  # noqa: E402
    NO_EKD,  # noqa: E402
    NO_MIR,  # noqa: E402
    compare_dims,
    covjson_to_xarray,
    get_test_data,
)

if not NO_EKD:
    from earthkit.data import from_source  # noqa


xr = pytest.importorskip("xarray")

REFS = [
    ({"grid": [10, 10]}, {"grid": [10, 10]}, {"step": 2, "latitude": 19, "longitude": 36}),
    ({"grid": "N32"}, {"grid": "N32"}, {"step": 2, "values": 6114}),
    ({"grid": "H4", "order": "nested"}, {"grid": "H4", "order": "nested"}, {"step": 2, "values": 192}),
    ({"grid": "H4", "order": "ring"}, {"grid": "H4"}, {"step": 2, "values": 192}),
]

REFS_SUBAREA = [
    (
        {"grid": [10, 10], "area": [80, -20, -10, 60]},
        {"grid": [10, 10], "area": [80, -20, -10, 60]},
        {"step": 2, "latitude": 10, "longitude": 9},
        [80, -20, -10, 60],
    ),
    (
        {"grid": [10, 10], "area": [85, -25, -10, 60]},
        {"grid": [10, 10], "area": [80, -20, -10, 60]},
        {"step": 2, "latitude": 10, "longitude": 9},
        [80, -20, -10, 60],
    ),
]


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.parametrize("out_grid,out_grid_ref,dims", REFS)
def test_regrid_xarray_from_ogg(out_grid, out_grid_ref, dims):

    ds_in = from_source("sample", "O32_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec == out_grid_ref


@pytest.mark.skip(reason="This test is currently failing")
@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.parametrize("out_grid,out_grid_ref,dims,area_ref", REFS_SUBAREA)
def test_regrid_xarray_from_ogg_to_subarea(out_grid, out_grid_ref, dims, area_ref):

    ds_in = from_source("sample", "O32_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec == out_grid_ref

    lat = r["latitude"].values
    lon = r["longitude"].values
    north = lat.max()
    south = lat.min()
    east = lon.max()
    west = lon.min()

    assert np.isclose(north, area_ref[0])
    assert np.isclose(south, area_ref[2])
    assert np.isclose(east, area_ref[3])
    assert np.isclose(west, area_ref[1])


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.parametrize(
    "out_grid,out_grid_ref,dims",
    REFS,
)
def test_regrid_xarray_from_h_nested(out_grid, out_grid_ref, dims):

    ds_in = from_source("sample", "H8_nested_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec == out_grid_ref


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.parametrize(
    "sample,out_grid,out_grid_ref,dims",
    [
        ("H8_nested_t2.grib2", {"grid": "N32"}, {"grid": "N32"}, {"step": 2, "values": 6114}),
        (
            "test.grib",
            {"grid": [10, 10]},
            {"grid": [10, 10], "area": [70, -20, 40, 40]},
            {"latitude": 4, "longitude": 7},
        ),
    ],
)
def test_regrid_xarray_transposed_dims(sample, out_grid, out_grid_ref, dims):
    ds_in = from_source("sample", sample).to_fieldlist()

    ds = ds_in.to_xarray()
    da = ds["2t"].transpose("longitude", "latitude", "values", "step", missing_dims="ignore")

    r = regrid(da, out_grid=out_grid, interpolation="linear")  # transposed dims (non-contiguous array)
    r_ref = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")  # original dim order (reference)

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec == out_grid_ref
    np.testing.assert_allclose(r.transpose(*r_ref.dims).values, r_ref.values)


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.parametrize(
    "out_grid,out_grid_ref,dims",
    REFS,
)
def test_regrid_xarray_dataset_from_h_nested(out_grid, out_grid_ref, dims):

    ds_in = from_source("sample", "H8_nested_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds, out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec == out_grid_ref


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_regrid_xarray_2d_1(lat_name, lon_name):
    # Dimensions:  (level: 2, lat: 3, lon: 3)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #   * lat      (lat) int64 24B 50 40 30
    #   * lon      (lon) int64 24B 0 10 20
    # Data variables:
    #     a        (level, lat, lon) int64 144B 11 12 13 21 22 23 ... 25 26 34 35 36

    import xarray as xr

    dims = {"level": 2, lat_name: 3, lon_name: 3}
    coords = {
        "level": np.array([700, 500]),
        lat_name: np.array([50, 40, 30]),
        lon_name: np.array([0, 10, 20]),
    }

    data = np.array(
        [
            [[11, 12, 13], [21, 22, 23], [31, 32, 33]],
            [[14, 15, 16], [24, 25, 26], [34, 35, 36]],
        ],
        dtype=np.float64,
    )

    a = xr.Variable(dims, data)
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    in_grid = {"grid": [10, 10], "area": [50, 0, 30, 20]}
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    out_dims = {"level": 2, "latitude": 5, "longitude": 5}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        [
            [11.0, 11.60916513, 12.0, 12.60709173, 13.0],
            [16.0, 16.64829292, 17.0, 17.56959565, 18.0],
            [21.0, 21.60735172, 22.0, 22.60876203, 23.0],
            [26.0, 26.57403598, 27.0, 27.63057889, 28.0],
            [31.0, 31.5, 32.0, 32.5, 33.0],
        ],
        [
            [14.0, 14.60916513, 15.0, 15.60709173, 16.0],
            [19.0, 19.64829292, 20.0, 20.56959565, 21.0],
            [24.0, 24.60735172, 25.0, 25.60876203, 26.0],
            [29.0, 29.57403598, 30.0, 30.63057889, 31.0],
            [34.0, 34.5, 35.0, 35.5, 36.0],
        ],
    ])

    ref_lat = np.array([50.0, 45.0, 40.0, 35.0, 30.0])
    ref_lon = np.array([0.0, 5.0, 10.0, 15.0, 20.0])

    assert np.allclose(r.values, ref_data)
    assert np.allclose(r.latitude.values, ref_lat)
    assert np.allclose(r.longitude.values, ref_lon)


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_regrid_xarray_2d_2(lat_name, lon_name):
    # Dimensions:  (level: 2, y: 3, x: 2)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #     lat      (y, x) int64 48B 50 50 40 40 30 30
    #     lon      (y, x) int64 48B 0 10 0 10 0 10
    # Dimensions without coordinates: y, x
    # Data variables:
    #     a        (level, y, x) int64 96B 11 12 21 22 31 32 14 15 24 25 34 35

    import xarray as xr

    dims = {"level": 2, "y": 3, "x": 2}
    coords = {
        "level": np.array([700, 500]),
        lat_name: (["y", "x"], np.array([[50, 50], [40, 40], [30, 30]])),
        lon_name: (["y", "x"], np.array([[0, 10], [0, 10], [0, 10]])),
    }

    data = np.array(
        [
            [[11, 12], [21, 22], [31, 32]],
            [[14, 15], [24, 25], [34, 35]],
        ],
        dtype=np.float64,
    )

    a = xr.Variable(dims, data)
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    in_grid = {
        "type": "unstructured_ll",
        "latitudes": [50.0, 50.0, 40.0, 40.0, 30.0, 30.0],
        "longitudes": [0.0, 10.0, 0.0, 10.0, 0.0, 10.0],
    }
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    out_dims = {"level": 2, "latitude": 37, "longitude": 72}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        [[11.0, 12.0], [21.0, 22.0], [31.0, 32.0]],
        [[14.0, 15.0], [24.0, 25.0], [34.0, 35.0]],
    ])

    ref_lat = np.array([50.0, 40.0, 30.0])
    ref_lon = np.array([0.0, 10.0])

    r_sub = r.sel(latitude=ref_lat, longitude=ref_lon)

    assert np.allclose(r_sub.values, ref_data)
    assert np.allclose(r_sub.latitude.values, ref_lat)
    assert np.allclose(r_sub.longitude.values, ref_lon)


@pytest.mark.parametrize("lat_name,lon_name", [("lat", "lon"), ("latitude", "longitude")])
def test_regrid_xarray_1d_1(lat_name, lon_name):
    # Dimensions:  (level: 2, values: 9)
    # Coordinates:
    #   * level    (level) int64 16B 700 500
    #     lat      (values) int64 72B 50 50 50 40 40 40 30 30 30
    #     lon      (values) int64 72B 0 10 20 0 10 20 0 10 20
    # Dimensions without coordinates: values
    # Data variables:
    #     a        (level, values) int64 144B 11 12 13 21 22 23 ... 24 25 26 34 35 36

    import xarray as xr

    dims = {"level": 2, "values": 9}
    coords = {
        "level": np.array([700, 500]),
        lat_name: ("values", np.array([50, 50, 50, 40, 40, 40, 30, 30, 30])),
        lon_name: ("values", np.array([0, 10, 20, 0, 10, 20, 0, 10, 20])),
    }

    data = np.array(
        [
            [11, 12, 13, 21, 22, 23, 31, 32, 33],
            [14, 15, 16, 24, 25, 26, 34, 35, 36],
        ],
        dtype=np.float64,
    )

    a = xr.Variable(dims, data)
    v = {"a": a}
    ds_in = xr.Dataset(v, coords=coords)

    in_grid = {
        "type": "unstructured_ll",
        "latitudes": [50.0, 50.0, 50.0, 40.0, 40.0, 40.0, 30.0, 30.0, 30.0],
        "longitudes": [0.0, 10.0, 20.0, 0.0, 10.0, 20.0, 0.0, 10.0, 20.0],
    }
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    out_dims = {"level": 2, "latitude": 37, "longitude": 72}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        [[11.0, 12.0, 13.0], [21.0, 22.30354886, 23.0], [31.0, 32.0, 33.0]],
        [[14.0, 15.0, 16.0], [24.0, 25.30354886, 26.0], [34.0, 35.0, 36.0]],
    ])

    ref_lat = np.array([50.0, 40.0, 30.0])
    ref_lon = np.array([0.0, 10.0, 20.0])

    r_sub = r.sel(latitude=ref_lat, longitude=ref_lon)

    assert np.allclose(r_sub.values, ref_data)
    assert np.allclose(r_sub.latitude.values, ref_lat)
    assert np.allclose(r_sub.longitude.values, ref_lon)


@pytest.mark.parametrize("in_grid", [None, {"grid": [30.0, 30.0]}])
def test_regrid_xarray_from_netcdf_ll_to_ll_1(in_grid):
    path = get_test_data("test_single.nc", subfolder="xr")
    ds = xr.open_dataset(path)
    da = ds["t2m"]

    out_grid = {"grid": [10, 10]}
    r = regrid(da, in_grid=in_grid, out_grid=out_grid, interpolation="nn")

    out_dims = {"latitude": 19, "longitude": 36}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        280.8106,
        280.8106,
        277.0606,
        277.0606,
        277.0606,
        284.4356,
        284.4356,
        284.4356,
        292.3106,
        292.3106,
        292.3106,
        274.8106,
        274.8106,
        274.8106,
        272.1856,
        272.1856,
        272.1856,
        273.9356,
        273.9356,
        273.9356,
        270.3106,
        270.3106,
        270.3106,
        272.8106,
        272.8106,
        272.8106,
        261.1856,
        261.1856,
        261.1856,
        264.3106,
        264.3106,
        264.3106,
        275.8106,
        275.8106,
        275.8106,
        280.8106,
    ])
    ref_lat = np.linspace(90.0, -90.0, 19)
    ref_lon = np.linspace(0.0, 350.0, 36)

    assert np.allclose(r.to_numpy()[2], ref_data)
    assert np.allclose(r.latitude.values, ref_lat)
    assert np.allclose(r.longitude.values, ref_lon)


def test_regrid_xarray_cordex_rotated_ll_to_ll():
    path = get_test_data("cordex.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    out_grid = {"grid": [10, 10]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    out_dims = {"time": 2, "latitude": 19, "longitude": 36}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        np.nan,
        4.9234657,
        np.nan,
        np.nan,
        3.0825672,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ])

    ref_lat = np.linspace(90.0, -90.0, 19)
    ref_lon = np.linspace(0.0, 350.0, 36)

    assert np.allclose(r["sfcWind"].to_numpy()[0][5], ref_data, equal_nan=True)
    assert np.allclose(r.latitude.values, ref_lat)
    assert np.allclose(r.longitude.values, ref_lon)


@pytest.mark.skipif(NO_COVJSON, reason="No covjsonkit available")
def test_regrid_xarray_covjson_unstructured_to_ll():
    path = get_test_data("points.covjson", subfolder="xr")
    ds = covjson_to_xarray(path)

    out_grid = {"grid": [10, 10]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    out_dims = {"datetimes": 1, "number": 1, "steps": 1, "latitude": 19, "longitude": 36}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([
        7.43803024,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        6.08695602,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
        7.43803024,
    ])

    ref_lat = np.linspace(90.0, -90.0, 19)
    ref_lon = np.linspace(0.0, 350.0, 36)

    assert np.allclose(r["10u"].to_numpy()[0][0][0][2], ref_data, equal_nan=True)
    assert np.allclose(r.latitude.values, ref_lat)
    assert np.allclose(r.longitude.values, ref_lon)
