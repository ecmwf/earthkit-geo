# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Regrid tests for xarray objects regridded by MIR through ``earthkit.geo.regrid``.

Three kinds of input are regridded: the ``xr.DataArray`` / ``xr.Dataset`` made from an
earthkit-data "sample" GRIB (O32_t2.grib2, H8_nested_t2.grib2, test.grib), the NetCDF and
covjson files of the "xr" test data subfolder (test_single.nc, efas.nc, cordex.nc,
points.covjson), and small synthetic ``xr.Dataset`` objects built by the tests themselves.

The shared parts come from ``regrid_common``: ``requires_mir`` and ``requires_ekd`` are
the skip marks of the sample tests, ``IN_GITHUB`` is re-exported there, and
``data_path("legacy", ..., subfolder="xr")`` is the single data access (it delegates to
``earthkit.geo.utils.testing.get_test_data``). ``NO_COVJSON``, ``compare_dims`` and
``covjson_to_xarray`` have no counterpart in ``regrid_common`` and are imported from
``earthkit.geo.utils.testing``.
"""

import numpy as np
import pytest
from regrid_common import IN_GITHUB, NO_EKD, data_path, requires_ekd, requires_mir

from earthkit.geo import regrid
from earthkit.geo.utils.testing import NO_COVJSON, compare_dims, covjson_to_xarray

if not NO_EKD:
    from earthkit.data import from_source


xr = pytest.importorskip("xarray")

# Being unable to run in the GitHub CI is a property of the tests of this file, and
# regrid_common.py only provides the MIR and earthkit-data marks.
skip_in_github = pytest.mark.skipif(IN_GITHUB, reason="Skipping test in GitHub CI")

# The output grids of test_regrid_xarray_from_sample: the requested out_grid, the
# gridspec the result is expected to report (its items must include these) and the
# expected sizes of the "step" and coordinate / "values" dimensions.
REFS = [
    ({"grid": [10, 10]}, {"grid": [10, 10]}, {"step": 2, "latitude": 19, "longitude": 36}),
    ({"grid": "N32"}, {"grid": "N32"}, {"step": 2, "values": 6114}),
    ({"grid": "H4", "order": "nested"}, {"grid": "H4", "order": "nested"}, {"step": 2, "values": 192}),
    ({"grid": "H4", "order": "ring"}, {"grid": "H4"}, {"step": 2, "values": 192}),
]
REFS_IDS = ["ll-10x10", "n32", "h4-nested", "h4-ring"]

# The same, for the regional output grids of test_regrid_xarray_from_ogg_to_subarea.
# The fourth element is the area (north, west, south, east) the result coordinates must
# span; the second case asks for an area larger than the source grid can cover.
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

# The inputs of test_regrid_xarray_from_sample: the earthkit-data sample GRIB, holding
# two fields, and the object handed to regrid() -- the name of its "2t" data array, or
# None for the whole Dataset.
SAMPLE_CASES = [
    pytest.param("O32_t2.grib2", "2t", id="ogg-da"),
    pytest.param("H8_nested_t2.grib2", "2t", id="h-nested-da"),
    pytest.param("H8_nested_t2.grib2", None, id="h-nested-ds"),
]

# A global {"grid": [10, 10]} result spans the poles: MIR returns 19 latitudes and 36
# longitudes for it.
OUT_DIMS_GLOBAL = {"latitude": 19, "longitude": 36}
REF_LAT_GLOBAL = np.linspace(90.0, -90.0, 19)
REF_LON_GLOBAL = np.linspace(0.0, 350.0, 36)

# The {"grid": [5, 5]} output over the area of the synthetic 3x3 field of the 2D / 1D
# tests: the linear regridding of its two levels (11..36 and 14..36, see
# _synthetic_dataset) onto the area [50, 0, 30, 20]. The 1D test's auto-detected input
# grid is the same grid, and the curved 2D test sees the first three longitudes of it.
REF_3X3_TO_5X5 = np.array([
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
REF_LAT_5 = np.array([50.0, 45.0, 40.0, 35.0, 30.0])
REF_LON_5 = np.array([0.0, 5.0, 10.0, 15.0, 20.0])

# The coordinates of the coarser grids the curved-coordinate tests select from their
# global result.
REF_LAT_3 = np.array([50.0, 40.0, 30.0])

# The parametrisations of the synthetic tests: both spellings of the coordinate names,
# the input grid detected by MIR or stated by the test, and the input grids of the
# NetCDF tests (None lets MIR detect the grid of the file).
LATLON_NAMES = [("lat", "lon"), ("latitude", "longitude")]
AUTO_IN_GRIDS = [False, True]
IN_GRID_CASES = [None, {"grid": [30.0, 30.0]}]


def _sample_fields(sample):
    """The fieldlist of an earthkit-data "sample" GRIB source."""
    return from_source("sample", sample).to_fieldlist()


def _open_xr(name):
    """The xr.Dataset of a file of the "xr" subfolder of the test data."""
    return xr.open_dataset(data_path("legacy", name, subfolder="xr"))


def _synthetic_dataset(dims, coords, data):
    """The two-level xr.Dataset with the single data variable "a" the tests regrid."""
    return xr.Dataset({"a": xr.Variable(dims, data)}, coords=coords)


def _assert_latlons(r, ref_lat, ref_lon):
    """Assert that a regridded object spans the reference latitude / longitude grids."""
    assert np.allclose(r.latitude.values, ref_lat)
    assert np.allclose(r.longitude.values, ref_lon)


@requires_mir
@requires_ekd
@pytest.mark.parametrize("out_grid,out_grid_ref,dims", REFS, ids=REFS_IDS)
@pytest.mark.parametrize("sample,var", SAMPLE_CASES)
def test_regrid_xarray_from_sample(out_grid, out_grid_ref, dims, sample, var):
    """Regrid a DataArray, or a whole Dataset, of a sample fieldlist onto every output grid."""
    ds_in = _sample_fields(sample)
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    obj = ds if var is None else ds[var]
    r = regrid(obj, out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec.items() >= out_grid_ref.items()


@pytest.mark.skip(reason="This test is currently failing")
@requires_mir
@requires_ekd
@pytest.mark.parametrize("out_grid,out_grid_ref,dims,area_ref", REFS_SUBAREA)
def test_regrid_xarray_from_ogg_to_subarea(out_grid, out_grid_ref, dims, area_ref):
    """Regrid the OGG sample onto a regional output grid and check the area it covers."""
    ds_in = _sample_fields("O32_t2.grib2")
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec.items() >= out_grid_ref.items()

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


@requires_mir
@requires_ekd
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
    """A transposed (non-contiguous) data array regrids like the original dim order."""
    ds = _sample_fields(sample).to_xarray()
    da = ds["2t"].transpose("longitude", "latitude", "values", "step", missing_dims="ignore")

    r = regrid(da, out_grid=out_grid, interpolation="linear")  # transposed dims (non-contiguous array)
    r_ref = regrid(ds["2t"], out_grid=out_grid, interpolation="linear")  # original dim order (reference)

    compare_dims(r, dims, sizes=True)

    assert r.earthkit.grid_spec.items() >= out_grid_ref.items()
    np.testing.assert_allclose(r.transpose(*r_ref.dims).values, r_ref.values)


@skip_in_github
@pytest.mark.parametrize("lat_name,lon_name", LATLON_NAMES)
@pytest.mark.parametrize("auto_in_grid", AUTO_IN_GRIDS)
def test_regrid_xarray_2d_coords_to_ll(lat_name, lon_name, auto_in_grid):
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

    data = np.array(
        [
            [[11, 12, 13], [21, 22, 23], [31, 32, 33]],
            [[14, 15, 16], [24, 25, 26], [34, 35, 36]],
        ],
        dtype=np.float64,
    )

    ds_in = _synthetic_dataset(dims, coords, data)

    in_grid = None if auto_in_grid else {"grid": [10, 10], "area": [50, 0, 30, 20]}
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    out_dims = {"level": 2, "latitude": 5, "longitude": 5}
    compare_dims(r, out_dims, sizes=True)

    assert np.allclose(r.values, REF_3X3_TO_5X5)
    _assert_latlons(r, REF_LAT_5, REF_LON_5)


@skip_in_github
@pytest.mark.parametrize("lat_name,lon_name", LATLON_NAMES)
@pytest.mark.parametrize("auto_in_grid", AUTO_IN_GRIDS)
def test_regrid_xarray_2d_curved_to_ll(lat_name, lon_name, auto_in_grid):
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

    data = np.array(
        [
            [[11, 12], [21, 22], [31, 32]],
            [[14, 15], [24, 25], [34, 35]],
        ],
        dtype=np.float64,
    )

    ds_in = _synthetic_dataset(dims, coords, data)

    in_grid = (
        None
        if auto_in_grid
        else {
            "latitudes": [50.0, 50.0, 40.0, 40.0, 30.0, 30.0],
            "longitudes": [0.0, 10.0, 0.0, 10.0, 0.0, 10.0],
        }
    )
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    if auto_in_grid:
        # the 6 points form a complete regular lat/lon grid (area [50, 0, 30, 10]), so the
        # auto-detected input grid is confined to that area, unlike the explicit unstructured
        # "latitudes"/"longitudes" in_grid below, which MIR interpolates over the whole globe
        out_dims = {"level": 2, "latitude": 5, "longitude": 3}
        compare_dims(r, out_dims, sizes=True)

        assert np.allclose(r.values, REF_3X3_TO_5X5[:, :, :3])
        _assert_latlons(r, REF_LAT_5, REF_LON_5[:3])
    else:
        out_dims = {"level": 2, "latitude": 37, "longitude": 72}
        compare_dims(r, out_dims, sizes=True)

        ref_data = np.array([
            [[11.0, 12.0], [21.0, 22.0], [31.0, 32.0]],
            [[14.0, 15.0], [24.0, 25.0], [34.0, 35.0]],
        ])

        ref_lon = np.array([0.0, 10.0])

        r_sub = r.sel(latitude=REF_LAT_3, longitude=ref_lon)

        assert np.allclose(r_sub.values, ref_data)
        _assert_latlons(r_sub, REF_LAT_3, ref_lon)


@skip_in_github
@pytest.mark.parametrize("lat_name,lon_name", LATLON_NAMES)
@pytest.mark.parametrize("auto_in_grid", AUTO_IN_GRIDS)
def test_regrid_xarray_1d_unstructured_to_ll(lat_name, lon_name, auto_in_grid):
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

    data = np.array(
        [
            [11, 12, 13, 21, 22, 23, 31, 32, 33],
            [14, 15, 16, 24, 25, 26, 34, 35, 36],
        ],
        dtype=np.float64,
    )

    ds_in = _synthetic_dataset(dims, coords, data)

    in_grid = (
        None
        if auto_in_grid
        else {
            "latitudes": [50.0, 50.0, 50.0, 40.0, 40.0, 40.0, 30.0, 30.0, 30.0],
            "longitudes": [0.0, 10.0, 20.0, 0.0, 10.0, 20.0, 0.0, 10.0, 20.0],
        }
    )
    out_grid = {"grid": [5, 5]}

    r = regrid(ds_in["a"], in_grid=in_grid, out_grid=out_grid, interpolation="linear")

    if auto_in_grid:
        # the 9 points form a complete regular lat/lon grid (area [50, 0, 30, 20]), identical
        # to the one used in test_regrid_xarray_2d_coords_to_ll, so the auto-detected input
        # grid is confined to that area, unlike the explicit unstructured
        # "latitudes"/"longitudes" in_grid below, which MIR interpolates over the whole globe
        out_dims = {"level": 2, "latitude": 5, "longitude": 5}
        compare_dims(r, out_dims, sizes=True)

        assert np.allclose(r.values, REF_3X3_TO_5X5)
        _assert_latlons(r, REF_LAT_5, REF_LON_5)
    else:
        out_dims = {"level": 2, "latitude": 37, "longitude": 72}
        compare_dims(r, out_dims, sizes=True)

        ref_data = np.array([
            [[11.0, 12.0, 13.0], [21.0, 22.30354886, 23.0], [31.0, 32.0, 33.0]],
            [[14.0, 15.0, 16.0], [24.0, 25.30354886, 26.0], [34.0, 35.0, 36.0]],
        ])

        ref_lon = np.array([0.0, 10.0, 20.0])

        r_sub = r.sel(latitude=REF_LAT_3, longitude=ref_lon)

        assert np.allclose(r_sub.values, ref_data)
        _assert_latlons(r_sub, REF_LAT_3, ref_lon)


@skip_in_github
@pytest.mark.parametrize("in_grid", IN_GRID_CASES)
def test_regrid_xarray_netcdf_ll_to_ll(in_grid):
    """Regrid the NetCDF global lat/lon field onto a coarser global lat/lon grid."""
    da = _open_xr("test_single.nc")["t2m"]

    out_grid = {"grid": [10, 10]}
    r = regrid(da, in_grid=in_grid, out_grid=out_grid, interpolation="nn")

    compare_dims(r, OUT_DIMS_GLOBAL, sizes=True)

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

    assert np.allclose(r.to_numpy()[2], ref_data)
    _assert_latlons(r, REF_LAT_GLOBAL, REF_LON_GLOBAL)


@skip_in_github
@pytest.mark.parametrize("in_grid", IN_GRID_CASES)
def test_regrid_xarray_netcdf_ll_to_points(in_grid):
    """Regrid the NetCDF global lat/lon field onto two lat/lon points."""
    da = _open_xr("test_single.nc")["t2m"]

    lats = [40.0, 50.0]
    lons = [10.0, 20.0]
    out_grid = {"latitudes": lats, "longitudes": lons}
    r = regrid(da, in_grid=in_grid, out_grid=out_grid, interpolation="nn")

    out_dims = {"values": len(lats)}
    compare_dims(r, out_dims, sizes=True)

    ref_data = np.array([307.1856, 277.0606])
    assert np.allclose(r.to_numpy(), ref_data)
    assert np.allclose(r.latitude.values, np.array(lats))
    assert np.allclose(r.longitude.values, np.array(lons))


@skip_in_github
# @pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_regrid_xarray_laea_to_ll():
    """Regrid the LAEA projected EFAS field onto a global lat/lon grid."""
    ds = _open_xr("efas.nc")

    da = ds["dis06"]

    out_grid = {"grid": [10, 10]}
    r = regrid(da, out_grid=out_grid, interpolation="nn")

    assert isinstance(r, xr.DataArray)

    compare_dims(r, OUT_DIMS_GLOBAL, sizes=True)

    ref_data = np.array([
        np.nan,
        2.35173828e02,
        1.15234375e-01,
        4.24804688e-01,
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
        np.nan,
    ])

    assert np.allclose(r.to_numpy()[3], ref_data, equal_nan=True)
    _assert_latlons(r, REF_LAT_GLOBAL, REF_LON_GLOBAL)


@skip_in_github
def test_regrid_xarray_cordex_rotated_ll_to_ll():
    """Regrid the CORDEX Dataset on a rotated lat/lon grid onto a global lat/lon grid."""
    ds = _open_xr("cordex.nc")

    out_grid = {"grid": [10, 10]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    out_dims = {"time": 2, **OUT_DIMS_GLOBAL}
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

    assert np.allclose(r["sfcWind"].to_numpy()[0][5], ref_data, equal_nan=True)
    _assert_latlons(r, REF_LAT_GLOBAL, REF_LON_GLOBAL)


@skip_in_github
@pytest.mark.skipif(NO_COVJSON, reason="No covjsonkit available")
def test_regrid_xarray_covjson_unstructured_to_ll():
    """Regrid the unstructured covjson Dataset onto a global lat/lon grid."""
    ds = covjson_to_xarray(data_path("legacy", "points.covjson", subfolder="xr"))

    out_grid = {"grid": [10, 10]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    out_dims = {"datetimes": 1, "number": 1, "steps": 1, **OUT_DIMS_GLOBAL}
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

    assert np.allclose(r["10u"].to_numpy()[0][0][0][2], ref_data, equal_nan=True)
    _assert_latlons(r, REF_LAT_GLOBAL, REF_LON_GLOBAL)
