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

The tests near the end of the file regrid a chunked, dask-backed object. Some of them use a
stub backend (:func:`_recording_backend`) instead of MIR, to watch the calls the handler
makes and the gridspec it is answered with; they need dask and the synthetic input, but
neither MIR nor earthkit-data. The last tests regrid onto a rotated output grid without an
out_dims argument, including the case where only the grid the backend returns dictates the
result's dimensions; the geography-builder and stub-backend checks among them need no MIR.

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

# Dask, like xarray, is an optional dependency. It has no NO_DASK flag in
# earthkit.geo.utils.testing to take one from, so the test file looks for it itself.
try:
    import dask.array

    NO_DASK = False
except ImportError:
    NO_DASK = True

# Being unable to run in the GitHub CI is a property of the tests of this file, and
# regrid_common.py only provides the MIR and earthkit-data marks.
skip_in_github = pytest.mark.skipif(IN_GITHUB, reason="Skipping test in GitHub CI")

# The mark of the tests that regrid a chunked (dask-backed) xarray object.
requires_dask = pytest.mark.skipif(NO_DASK, reason="No dask available")

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
@pytest.mark.download
@pytest.mark.timeout(180)
def test_regrid_xarray_laea_dataset_to_ll():
    """Regrid the whole LAEA projected EFAS Dataset onto a lat/lon grid over its area.

    A Dataset keeps the scalar grid mapping variable of the file, so its variables are
    detected on a ``regular_xy`` grid in the plane of the PROJ projection of that mapping.
    MIR cannot read that gridspec as an input grid, so the xarray regrid handler hands the
    backend the same points as an unstructured lat/lon grid (via
    :func:`~earthkit.geo.grids._regrid.data.utils.scatter_projected_grid`). This checks
    that the scattered points are the grid of the data: every finite output point carries
    the value of its nearest input point, with the input points read from the file's own 2D
    coordinates, like the CORDEX test below. The output grid covers the area of the sample
    instead of the whole globe, for the reason given in that test.
    """
    ds = _open_xr("efas.nc")[["dis06", "upArea", "lambert_azimuthal_equal_area"]]

    from earthkit.geo.grids.xarray.loader import variables

    # the fallback is only exercised when the input grid is the detected projected one
    assert [v.ek_grid.type for v in variables(ds)] == ["regular_xy", "regular_xy"]

    out_grid = {"grid": [1, 1], "area": [60, -10, 35, 40]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    assert isinstance(r, xr.Dataset)

    out_dims = {"latitude": 26, "longitude": 51}
    compare_dims(r, out_dims, sizes=True)
    _assert_latlons(r, np.linspace(60.0, 35.0, 26), np.linspace(-10.0, 40.0, 51))

    from scipy.spatial import cKDTree

    def _xyz(lats, lons):
        lats = np.radians(lats)
        lons = np.radians(lons)
        return np.stack([np.cos(lats) * np.cos(lons), np.cos(lats) * np.sin(lons), np.sin(lats)], axis=-1)

    # the input points and values, in the (y, x) row-major order of the projected grid
    in_xyz = _xyz(ds["latitude"].values.reshape(-1), ds["longitude"].values.reshape(-1))
    grid_lat, grid_lon = np.meshgrid(r["latitude"].values, r["longitude"].values, indexing="ij")
    out_xyz = _xyz(grid_lat.reshape(-1), grid_lon.reshape(-1))
    _, nearest = cKDTree(in_xyz).query(out_xyz, k=1)

    for name in ("dis06", "upArea"):
        values = r[name].to_numpy().reshape(-1)
        finite = np.isfinite(values)
        assert int(finite.sum()) == 781
        expected = ds[name].to_numpy().reshape(-1)[nearest[finite]]
        assert np.max(np.abs(values[finite] - expected)) == 0.0


@skip_in_github
@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(180)
def test_regrid_xarray_polar_stereo_north_to_ll():
    """Regrid the polar stereographic AMSR2 sample onto a regional lat/lon grid.

    The sample keeps its grid as two 1D projected axes and a CF grid mapping described by
    its ``crs_wkt`` alone, so its variable is detected on the ``regular_xy`` grid of that
    PROJ projection. Like the efas Dataset above, MIR cannot read that gridspec as an input
    grid, so the xarray regrid handler hands the backend the same points as an
    unstructured lat/lon grid (via
    :func:`~earthkit.geo.grids._regrid.data.utils.scatter_projected_grid`). A
    ``xr.DataArray`` only keeps the grid when the mapping variable travels with it as a
    scalar coordinate: without it nothing interprets the projected axes any more and the
    detection raises ``NotImplementedError``, as it does for the x/y files with no grid
    mapping. The output grid covers part of the data area instead of the whole globe, for
    the reason given in the CORDEX test below. The sample is a sea ice field with missing
    sea, so the finite pattern is checked too: an output point carries a value exactly when
    the projection covers the axes at it and the nearest input point found below has one.
    """
    ds = _open_xr("polar_stereo_north_amsr2.nc")

    from earthkit.geo.grids.xarray.loader import variables

    # the fallback is only exercised when the input grid is the detected projected one
    assert [v.ek_grid.type for v in variables(ds)] == ["regular_xy"]

    out_grid = {"grid": [2, 2], "area": [80, -170, 60, -100]}
    r = regrid(ds, out_grid=out_grid, interpolation="nearest-neighbour")

    assert isinstance(r, xr.Dataset)

    out_dims = {"latitude": 11, "longitude": 36}
    compare_dims(r, out_dims, sizes=True)
    _assert_latlons(r, np.linspace(80.0, 60.0, 11), np.linspace(-170.0, -100.0, 36))

    # the DataArray with its grid mapping kept as a scalar coordinate regrids identically
    da = ds["z"].assign_coords(polar_stereographic=ds["polar_stereographic"])
    r_da = regrid(da, out_grid=out_grid, interpolation="nearest-neighbour")

    assert isinstance(r_da, xr.DataArray)
    compare_dims(r_da, out_dims, sizes=True)
    assert np.array_equal(r_da.to_numpy(), r["z"].to_numpy(), equal_nan=True)

    values = r["z"].to_numpy()

    from pyproj import CRS, Transformer

    # the output points in the plane of the projection of the file itself
    transformer = Transformer.from_crs(
        "EPSG:4326", CRS.from_wkt(str(ds["polar_stereographic"].attrs["crs_wkt"])), always_xy=True
    )
    x = ds["x"].values
    y = ds["y"].values
    grid_lon, grid_lat = np.meshgrid(r["longitude"].values, r["latitude"].values)
    px, py = transformer.transform(grid_lon, grid_lat)

    # the nearest cell of the file's regular axes, and the offset of the point from it:
    # only points far from a cell boundary have an unambiguous nearest input point
    fx = (px - x[0]) / (x[1] - x[0])
    fy = (py - y[0]) / (y[1] - y[0])
    ix = np.rint(fx).astype(int)
    iy = np.rint(fy).astype(int)
    inside = (ix >= 0) & (ix < x.size) & (iy >= 0) & (iy < y.size)
    nearest = ds["z"].to_numpy()[np.clip(iy, 0, y.size - 1), np.clip(ix, 0, x.size - 1)]

    finite = np.isfinite(values)
    assert np.array_equal(finite, inside & np.isfinite(nearest))
    assert float(values[finite].min()) >= 0.0
    assert float(values[finite].max()) <= 100.0

    unambiguous = finite & (np.abs(fx - ix) < 0.49) & (np.abs(fy - iy) < 0.49)
    assert int(unambiguous.sum()) >= 150
    assert np.max(np.abs(values[unambiguous] - nearest[unambiguous])) == 0.0


# A Lambert azimuthal equal-area projection on GRS80, with the false easting
# and northing of the efas sample, described the way CF files describe a
# projection (same as the dataset of test_regrid_xarray_laea_dataset_to_ll).
LAEA_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


def _synthetic_laea_dataset():
    """A small synthetic Dataset on the plane of a LAEA projection.

    Returns
    -------
    xr.Dataset
        The dataset, with its data variable "a" over the 1D ``x``/``y`` axes of
        the projection plane (5 km, ``y`` descending), the 2D
        ``latitude``/``longitude`` coordinates pyproj computes for them, and
        the scalar CF grid mapping variable naming the projection.
    """
    from pyproj import CRS, Transformer

    x_values = 4321000.0 + 5000.0 * np.arange(6)
    y_values = 3210000.0 - 5000.0 * np.arange(4)

    transformer = Transformer.from_crs(CRS.from_proj4(LAEA_PROJ4), "EPSG:4326", always_xy=True)
    xv, yv = np.meshgrid(x_values, y_values)
    lon2d, lat2d = transformer.transform(xv, yv)

    ds = xr.Dataset(
        {
            "a": xr.Variable(
                ["y", "x"],
                np.arange(24, dtype=np.float64).reshape(4, 6),
                attrs={"grid_mapping": "mapping"},
            )
        },
        coords={
            "x": (["x"], x_values, {"standard_name": "projection_x_coordinate", "units": "m"}),
            "y": (["y"], y_values, {"standard_name": "projection_y_coordinate", "units": "m"}),
            "latitude": (["y", "x"], lat2d, {"standard_name": "latitude", "units": "degrees_north"}),
            "longitude": (["y", "x"], lon2d, {"standard_name": "longitude", "units": "degrees_east"}),
        },
    )
    ds["mapping"] = xr.DataArray(
        0.0,
        attrs={"grid_mapping_name": "lambert_azimuthal_equal_area", "proj4_params": LAEA_PROJ4},
    )
    return ds


@requires_mir
@skip_in_github
def test_regrid_xarray_projected_explicit_grid_transposed_dims():
    """A projected field ordered (x, y) regrids like its (y, x) original on an explicit grid.

    A structured explicit eckit grid (here the regular_xy one of the LAEA
    projection) iterates its points with its x axis fastest, so the loader
    reorders the (x, y) grid dimensions of the data to (y, x) before handing
    it to the backend; without that the transposed field is read point by
    point in the wrong order (see
    :func:`earthkit.geo.grids.xarray.loader.adjust_variable_dim_order`).
    """
    ds = _synthetic_laea_dataset()

    from earthkit.geo.grids.xarray.loader import variables

    # the detected projected grid, given explicitly as the input grid
    ek_grid = variables(ds)[0].ek_grid
    assert ek_grid.type == "regular_xy"

    # the data transposed onto the (x, y) dim order, the coordinates left as they are
    ds_t = ds.assign(a=ds["a"].transpose("x", "y"))

    assert list(variables(ds_t, user_ek_grid=ek_grid)[0].geo_dims) == ["y", "x"]

    out_grid = {"grid": [0.5, 0.5], "area": [60, -10, 45, 20]}

    ref = regrid(ds, in_grid=ek_grid, out_grid=out_grid, interpolation="nn")["a"].values
    r = regrid(ds_t, in_grid=ek_grid, out_grid=out_grid, interpolation="nn")["a"].values

    np.testing.assert_allclose(r, ref, equal_nan=True)


@requires_mir
@skip_in_github
def test_regrid_xarray_projected_unstructured_explicit_grid_transposed_dims():
    """A transposed field regrids identically on its detected and its explicit unstructured grid.

    An unstructured eckit grid holds its points in the order the data was
    laid out, not the x-fastest order of a structured one, so the loader must
    leave the (x, y) dimensions of the variable alone when that is the given
    input grid (see
    :func:`earthkit.geo.grids.xarray.loader.adjust_variable_dim_order`).
    """
    ds_t = _synthetic_laea_dataset().transpose("x", "y")

    from earthkit.geo.grids.xarray.loader import variables

    auto = variables(ds_t)[0]
    assert auto.ek_grid.type == "unstructured_ll"
    assert list(variables(ds_t, user_ek_grid=auto.ek_grid)[0].geo_dims) == ["x", "y"]

    out_grid = {"grid": [0.5, 0.5], "area": [60, -10, 45, 20]}

    ref = regrid(ds_t, out_grid=out_grid, interpolation="nn")["a"].values
    r = regrid(ds_t, in_grid=auto.ek_grid, out_grid=out_grid, interpolation="nn")["a"].values

    np.testing.assert_allclose(r, ref, equal_nan=True)


def _synthetic_rotated_dataset():
    """A small synthetic Dataset on a CF rotated lat/lon grid.

    The true latitudes and longitudes are those of the eckit-geo grid carrying
    the rotation, so its detected gridspec reproduces the points exactly (see
    the same construction in tests/xarray_grid/test_xarray_grid.py).

    Returns
    -------
    xr.Dataset
        The dataset, with its data variable "a" over the descending 1D
        ``rlat`` and ascending 1D ``rlon`` axes of the rotated frame, the 2D
        ``latitude``/``longitude`` coordinates of their rotated grid, and the
        scalar CF grid mapping variable naming the rotation.
    """
    from eckit.geo import Grid

    mapping = {
        "grid_mapping_name": "rotated_latitude_longitude",
        "grid_north_pole_latitude": 40.0,
        "grid_north_pole_longitude": -170.0,
    }
    rlat_values = np.array([30.0, 20.0, 10.0])  # descending, the default scan
    rlon_values = np.array([-20.0, -10.0, 0.0, 10.0])

    # eckit describes the rotation by its south pole, the antipode of the CF
    # north pole, as a [longitude, latitude] pair
    spec = {
        "grid": [10.0, 10.0],
        "area": [rlat_values[0], rlon_values[0], rlat_values[-1], rlon_values[-1]],
        "projection": {"type": "rotation", "south_pole": [10.0, -40.0]},
    }
    lat_flat, lon_flat = (np.asarray(a, dtype=np.float64) for a in Grid(spec).to_latlons())
    lat2d = lat_flat.reshape(3, 4)
    lon2d = lon_flat.reshape(3, 4)

    ds = xr.Dataset(
        {
            "a": xr.Variable(
                ["rlat", "rlon"],
                np.arange(12, dtype=np.float64).reshape(3, 4),
                attrs={"grid_mapping": "mapping"},
            )
        },
        coords={
            "rlat": (["rlat"], rlat_values, {"standard_name": "grid_latitude", "units": "degrees"}),
            "rlon": (["rlon"], rlon_values, {"standard_name": "grid_longitude", "units": "degrees"}),
            "latitude": (["rlat", "rlon"], lat2d, {"standard_name": "latitude", "units": "degrees_north"}),
            "longitude": (["rlat", "rlon"], lon2d, {"standard_name": "longitude", "units": "degrees_east"}),
        },
    )
    ds["mapping"] = xr.DataArray(0.0, attrs=mapping)
    return ds


@requires_mir
@skip_in_github
def test_regrid_xarray_rotated_explicit_outer_order_grid_transposed_dims():
    """A transposed rotated field regrids identically on its detected and its explicit grid.

    The auto path detects the transposed field's grid as a regular_ll with an
    order that scans the longitude (the x of the rotation) outermost: the data
    already lies in that point order, so the loader must leave its (rlon,
    rlat) dimensions alone for that explicit grid too (see
    :func:`earthkit.geo.grids.xarray.loader._scans_x_fastest`).
    """
    ds_t = _synthetic_rotated_dataset().transpose("rlon", "rlat")

    from earthkit.geo.grids.xarray.loader import variables

    auto = variables(ds_t)[0]
    assert auto.ek_grid.order == "j-i+"
    assert list(variables(ds_t, user_ek_grid=auto.ek_grid)[0].geo_dims) == ["rlon", "rlat"]

    out_grid = {"grid": [0.5, 0.5], "area": [65, -25, 20, 65]}

    ref = regrid(ds_t, out_grid=out_grid, interpolation="nn")["a"].values
    r = regrid(ds_t, in_grid=auto.ek_grid, out_grid=out_grid, interpolation="nn")["a"].values

    np.testing.assert_allclose(r, ref, equal_nan=True)


@skip_in_github
# @pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_regrid_xarray_laea_to_ll():
    """Regrid the LAEA projected EFAS field onto a global lat/lon grid.

    The input is a DataArray, which loses the scalar grid mapping variable on the way to
    the loader, so its grid is still detected from its 2D latitudes/longitudes as
    unstructured. The same file regridded as a whole Dataset is detected as a
    ``regular_xy`` grid in the plane of the file's PROJ projection, which no backend can
    read as an input grid: the xarray regrid handler converts it into its unstructured
    points before calling the backend (see
    :func:`test_regrid_xarray_laea_dataset_to_ll`).
    """
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
    """Regrid the CORDEX Dataset on a rotated lat/lon grid onto a lat/lon grid over its area.

    The output grid covers the true (unrotated) footprint of the sample instead of the
    whole globe: this MIR cannot regrid a regional gridded input onto a global grid (its
    matrix methods raise where the unstructured path would return all-the-NaN output).
    That is a MIR limitation of gridded inputs, unrelated to the rotation.
    """
    ds = _open_xr("cordex.nc")

    out_grid = {"grid": [0.25, 0.25], "area": [42.25, 10.0, 39.75, 13.0]}
    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    # MIR crops the requested area to the 0.25 grid points inside the input footprint
    out_dims = {"time": 2, "latitude": 9, "longitude": 9}
    compare_dims(r, out_dims, sizes=True)

    ref_lat = np.linspace(42.0, 40.0, 9)
    ref_lon = np.linspace(10.5, 12.5, 9)
    _assert_latlons(r, ref_lat, ref_lon)

    result = r["sfcWind"].to_numpy()

    # the middle (41 degree) row of the first time slice
    ref_data = np.array([
        np.nan,
        6.6690645,
        7.027647,
        6.617398,
        6.356368,
        5.8406057,
        5.609106,
        5.520202,
        5.442434,
    ])
    assert np.allclose(result[0][4], ref_data, equal_nan=True)

    # every finite output point carries the value of its nearest input point, computed
    # from the file's own 2D coordinates: this checks that the grid MIR was given from
    # the auto-detected rotated gridspec is the grid the file describes
    def _xyz(lats, lons):
        lats = np.radians(lats)
        lons = np.radians(lons)
        return np.stack([np.cos(lats) * np.cos(lons), np.cos(lats) * np.sin(lons), np.sin(lats)], axis=-1)

    in_lat = ds["lat"].values.reshape(-1)
    in_lon = ds["lon"].values.reshape(-1)
    grid_lat, grid_lon = np.meshgrid(r["latitude"].values, r["longitude"].values, indexing="ij")
    out_xyz = _xyz(grid_lat.reshape(-1), grid_lon.reshape(-1))[:, None]
    in_xyz = _xyz(in_lat, in_lon)[None, :]
    dist = np.linalg.norm(out_xyz - in_xyz, axis=2)
    nearest = dist.argmin(axis=1).reshape(grid_lat.shape)

    for time in range(2):
        values = result[time]
        finite = np.isfinite(values)
        assert int(finite.sum()) == 69
        expected = ds["sfcWind"].values[time].reshape(-1)[nearest[finite]]
        assert np.max(np.abs(values[finite] - expected)) == 0.0


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


# ---------------------------------------------------------------------------------------
# The lazy path: a chunked (dask-backed) input regrids into a result that is itself
# lazy, and whose dimensions, coordinates and gridspec already describe the grid the
# backend produces before any of it is computed. See issue #144.
# ---------------------------------------------------------------------------------------

# The output grids of test_regrid_xarray_dask_from_sample: the requested out_grid, the
# extra regrid() keyword arguments it needs, and the sizes of the result. The second
# case is the one the requested grid cannot answer for: MIR extends its area to the
# pole and regrids onto 11 latitudes where the eckit grid asked for has 10.
LAZY_REFS = [
    ({"grid": [5, 5]}, {}, {"step": 2, "latitude": 37, "longitude": 72}),
    ({"grid": [10, 10], "area": [80, -20, -10, 60]}, {}, {"step": 2, "latitude": 11, "longitude": 9}),
    (
        {
            "grid": [5, 5],
            "area": [65, -25, 20, 65],
            "projection": {"type": "rotation", "south_pole": [10.0, -40.0]},
        },
        {},
        {"step": 2, "y": 10, "x": 19},
    ),
]
LAZY_REFS_IDS = ["ll-5x5", "ll-10x10-area", "rotated-ll"]

# The grid the stub backend of the _recording_backend tests reports as the one it has
# regridded onto, whatever it was asked for, the grid those tests ask for, and the
# sizes of the result: the stub's 10 x 9 regional grid, twice as far from the 37 x 72
# global grid of the requested {"grid": [5, 5]} as it can be.
STUB_OUT_SPEC = {"grid": [10, 10], "area": [80, -20, -10, 60]}
STUB_OUT_GRID = {"grid": [5, 5]}
STUB_OUT_SIZES = {"level": 2, "latitude": 10, "longitude": 9}

# The geographical slice of _tiny_ll_dataset(), in the dimension order the handler
# hands it to a backend in, and the number of those slices the variable holds.
STUB_SLICE_SHAPE = (3, 4)
STUB_SLICES = 2


def _tiny_ll_dataset():
    """A small synthetic lat/lon Dataset the stub backend tests regrid.

    Returns
    -------
    xr.Dataset
        The dataset, with the single data variable ``a`` over the dimensions
        ``level`` (2), ``latitude`` (3) and ``longitude`` (4) -- the two geographical
        dimensions of different length, so the order they are passed to a backend in
        is visible in the shape of the array it gets -- and the values 1 .. 24, so
        that the sum of a slice identifies the slice.
    """
    n_level, n_lat, n_lon = STUB_SLICES, STUB_SLICE_SHAPE[0], STUB_SLICE_SHAPE[1]

    return _synthetic_dataset(
        {"level": n_level, "latitude": n_lat, "longitude": n_lon},
        {
            "level": np.arange(n_level),
            "latitude": np.linspace(50.0, 30.0, n_lat),
            "longitude": np.linspace(0.0, 20.0, n_lon),
        },
        np.arange(1.0, n_level * n_lat * n_lon + 1.0).reshape(n_level, n_lat, n_lon),
    )


def _recording_backend(out_spec):
    """A stub backend that records the slices it regrids and answers with a fixed grid.

    The lazy path of the xarray handler has to be tested against a backend a test can
    watch, and one that does not need MIR: this one reports ``out_spec`` as the grid it
    regridded onto, whichever grid it was asked for, and fills an array of that grid's
    shape with the sum of the slice it was given. That value is a pure function of the
    *whole* geographical slice, so a call that only saw part of the data (a core
    dimension split over blocks), or a slice of the wrong shape, shows up in it.

    Parameters
    ----------
    out_spec : dict
        The grid spec every call reports as its output grid.

    Returns
    -------
    Backend
        The backend, whose ``calls`` list holds the values array of every ``regrid()``
        call it has made, in the order they were made. A scheduler that forks the
        calls off into worker processes records their own copies, so the tests that
        read ``calls`` also ask for :func:`in_process_scheduler`.
    """
    import numpy as np

    from earthkit.geo.grids._regrid.backends import Backend
    from earthkit.geo.grids._regrid.data.utils import create_grid_object

    class _StubBackend(Backend):
        name = "recording-stub"

        def __init__(self, out_spec):
            self.out_spec = out_spec
            self.shape = create_grid_object(out_spec).shape
            self.calls = []

        def regrid(self, values, in_grid, out_grid, **kwargs):
            values = np.asarray(values)
            self.calls.append(values.copy())
            return np.full(self.shape, float(values.sum())), self.out_spec

    return _StubBackend(out_spec)


@pytest.fixture
def in_process_scheduler():
    """Run the dask tasks of a test in this process, whatever the configured scheduler is.

    The handler binds the input and output grids to the backend call, and an eckit
    ``Grid`` cannot be pickled, so a chunked object can only be regridded by a
    scheduler that runs those calls in this process. The tests of the call-recording
    stub need it for a second reason: a worker process would record its own calls.
    """
    with dask.config.set(scheduler="threads"):
        yield


def _assert_lazy_like_eager(lazy, eager, dims):
    """Assert that a lazy regridded object is the lazy twin of the eager result.

    Checks that ``lazy`` is dask-backed and that it describes ``dims`` before anything
    of it is computed -- the sizes of the graph, not of the values -- and that it
    agrees with ``eager`` on the dimensions, the coordinates and the attributes,
    including the grid spec earthkit-data results carry in their ``_earthkit``
    attribute. Computing it must then give the eager values.

    Parameters
    ----------
    lazy : xr.DataArray or xr.Dataset
        The regridded chunked object.
    eager : xr.DataArray or xr.Dataset
        The same object regridded from a non-chunked input.
    dims : dict
        The dimension sizes the result must have (see ``compare_dims``).
    """
    for name, var in lazy.data_vars.items() if isinstance(lazy, xr.Dataset) else [(lazy.name, lazy)]:
        assert isinstance(var.data, dask.array.Array), f"{name} is not dask-backed"

    if isinstance(lazy, xr.Dataset):
        assert list(lazy.data_vars) == list(eager.data_vars)
    else:
        assert lazy.dims == eager.dims

    compare_dims(lazy, dims, sizes=True)

    assert dict(lazy.sizes) == dict(eager.sizes)
    assert list(lazy.coords) == list(eager.coords)
    for name in lazy.coords:
        np.testing.assert_allclose(lazy[name].values, eager[name].values)

    # the attributes, with the _earthkit one of the earthkit-data inputs in them
    assert lazy.attrs == eager.attrs
    if isinstance(lazy, xr.Dataset):
        assert {k: v.attrs for k, v in lazy.data_vars.items()} == {k: v.attrs for k, v in eager.data_vars.items()}

    has_earthkit = "_earthkit" in lazy.attrs or (
        isinstance(lazy, xr.Dataset) and any("_earthkit" in v.attrs for v in lazy.data_vars.values())
    )
    if has_earthkit:
        assert lazy.earthkit.grid_spec.items() == eager.earthkit.grid_spec.items()

    xr.testing.assert_allclose(lazy.compute(), eager)


@requires_mir
@requires_ekd
@requires_dask
@pytest.mark.parametrize("out_grid,out_kwargs,dims", LAZY_REFS, ids=LAZY_REFS_IDS)
def test_regrid_xarray_dask_from_sample(in_process_scheduler, out_grid, out_kwargs, dims):
    """A chunked DataArray and a chunked Dataset regridded by MIR stay lazy and agree with the eager result."""
    ds = _sample_fields("O32_t2.grib2").to_xarray()
    # the geographical dimension of the sample is split over several chunks, which the
    # regridding has to put back together for every slice it hands to MIR
    chunks = {"values": 1000}

    eager = regrid(ds["2t"], out_grid=out_grid, interpolation="linear", **out_kwargs)
    lazy = regrid(ds["2t"].chunk(chunks), out_grid=out_grid, interpolation="linear", **out_kwargs)
    _assert_lazy_like_eager(lazy, eager, dims)

    # the same, for the whole Dataset
    eager = regrid(ds, out_grid=out_grid, interpolation="linear", **out_kwargs)
    lazy = regrid(ds.chunk(chunks), out_grid=out_grid, interpolation="linear", **out_kwargs)
    _assert_lazy_like_eager(lazy, eager, dims)


@requires_dask
def test_regrid_xarray_dask_backend_changes_grid(in_process_scheduler):
    """A lazy regridding takes its sizes, coordinates and grid from the backend's spec."""
    ds = _tiny_ll_dataset()
    eager_backend = _recording_backend(STUB_OUT_SPEC)
    lazy_backend = _recording_backend(STUB_OUT_SPEC)

    eager = regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=eager_backend)
    lazy = regrid(ds["a"].chunk({"level": 1}), out_grid=STUB_OUT_GRID, backend=lazy_backend)

    # the graph of the lazy result, before anything is computed, already has the size
    # of the grid the backend reported, not the 37 x 72 of the one that was requested
    compare_dims(lazy, STUB_OUT_SIZES, sizes=True)
    _assert_lazy_like_eager(lazy, eager, STUB_OUT_SIZES)

    # the geography of the result is the backend's grid: 10 latitudes and 9 longitudes
    # over its area, with the values the stub fills every point of
    np.testing.assert_allclose(lazy.latitude.values, np.linspace(80.0, -10.0, 10))
    np.testing.assert_allclose(lazy.longitude.values, np.linspace(-20.0, 60.0, 9))
    sums = np.array([float(s.sum()) for s in ds["a"]])
    np.testing.assert_allclose(lazy.to_numpy(), np.broadcast_to(sums[:, None, None], (STUB_SLICES, 10, 9)))


@requires_dask
def test_regrid_xarray_dask_transposed_dims(in_process_scheduler):
    """A chunked input whose dimensions are reordered is probed in the order the backend gets.

    The variable is held as (longitude, latitude, level) while the grid dimensions the
    loader detects for it are (latitude, longitude), so the array the backend sees, at
    the probe as well as in the graph, is the 3 x 4 geographical slice in that order.
    """
    from earthkit.geo.grids.xarray.loader import variables

    ds = _tiny_ll_dataset()
    ds = ds.assign(a=ds["a"].transpose("longitude", "latitude", "level"))
    assert list(variables(ds)[0].geo_dims) == ["latitude", "longitude"]

    eager_backend = _recording_backend(STUB_OUT_SPEC)
    lazy_backend = _recording_backend(STUB_OUT_SPEC)

    eager = regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=eager_backend)
    lazy = regrid(ds["a"].chunk({"level": 1}), out_grid=STUB_OUT_GRID, backend=lazy_backend)

    # the probe passes one geographical slice, in the grid dimension order and without
    # the level dimension, exactly like the eager calls and the graph do
    assert [c.shape for c in lazy_backend.calls] == [STUB_SLICE_SHAPE]
    assert [c.shape for c in eager_backend.calls] == [STUB_SLICE_SHAPE] * STUB_SLICES

    _assert_lazy_like_eager(lazy, eager, STUB_OUT_SIZES)
    assert [c.shape for c in lazy_backend.calls] == [STUB_SLICE_SHAPE] * (STUB_SLICES + 1)


@requires_dask
def test_regrid_xarray_dask_probe_only(in_process_scheduler):
    """Building the lazy graph costs one dummy slice; the data is regridded when it computes."""
    ds = _tiny_ll_dataset()
    backend = _recording_backend(STUB_OUT_SPEC)

    lazy = regrid(ds["a"].chunk({"level": 1}), out_grid=STUB_OUT_GRID, backend=backend)

    # while the graph is built the backend has only seen the probe: the values array is
    # still all zeros, because no data of the variable has reached it yet
    assert [c.shape for c in backend.calls] == [STUB_SLICE_SHAPE]
    assert float(backend.calls[0].sum()) == 0.0
    assert backend.calls[0].dtype == ds["a"].dtype

    lazy.compute()

    # one call per slice of the variable, on top of the probe, and every one of them
    # saw the whole slice it regrids. The scheduler may run the slices in any order,
    # so only the multiset of their sums is compared here; that they end up in the
    # right place of the result is what test_regrid_xarray_dask_backend_changes_grid
    # checks.
    assert [c.shape for c in backend.calls] == [STUB_SLICE_SHAPE] * (STUB_SLICES + 1)
    assert sorted(float(c.sum()) for c in backend.calls[1:]) == sorted(float(s.sum()) for s in ds["a"])


def test_regrid_xarray_eager_no_probe():
    """A non-chunked input regrids exactly as before: one backend call per slice, on the data.

    The probe the lazy path makes (see test_regrid_xarray_dask_probe_only) is not part
    of the eager path: it runs the backend once per slice of the variable, on the
    values themselves, and nothing else.
    """
    ds = _tiny_ll_dataset()
    backend = _recording_backend(STUB_OUT_SPEC)

    r = regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=backend)

    assert r.chunks is None
    # no probe: the eager path has never needed one, and the real calls tell it what
    # the lazy path can only learn in advance
    assert [c.shape for c in backend.calls] == [STUB_SLICE_SHAPE] * STUB_SLICES
    assert all(float(c.sum()) != 0.0 for c in backend.calls)

    compare_dims(r, STUB_OUT_SIZES, sizes=True)
    sums = np.array([float(s.sum()) for s in ds["a"]])
    np.testing.assert_allclose(r.to_numpy(), np.broadcast_to(sums[:, None, None], (STUB_SLICES, 10, 9)))


# ---------------------------------------------------------------------------------------
# Rotated (non-regular-mesh) output grids regridded without an out_dims argument: the
# output geography is made of 2D latitude/longitude coordinates on the (y, x) dimensions
# of the grid, so the dimensions the handler derives have to agree with it -- for the
# requested grid, by classifying its mesh strictly, and for the grid the backend really
# returns, by renaming the result's dimensions to its geography (see
# XarrayDataHandler.regrid).
# ---------------------------------------------------------------------------------------

# The rotated output grid of these tests: a 5 x 5 degree grid over the area
# [65, -25, 20, 65] of the frame whose rotation has its south pole at (10, -40),
# which regrids onto 10 rows and 19 columns.
ROT_OUT_GRID = {
    "grid": [5, 5],
    "area": [65, -25, 20, 65],
    "projection": {"type": "rotation", "south_pole": [10.0, -40.0]},
}
ROT_OUT_SIZES = {"y": 10, "x": 19}

# The request of the review reproduction: a rotated 10 x 10 degree grid over an area
# with two latitudes and two longitudes, whose 2-point rows and columns an equal-
# spacing test cannot tell apart from a regular mesh; MIR extends it to a 3 x 2 grid.
ROT_OUT_GRID_2X2 = {
    "grid": [10, 10],
    "area": [80, 0, 70, 10],
    "projection": {"type": "rotation", "south_pole": [10.0, -40.0]},
}


def _global_ll_dataarray():
    """The 19 x 36 global lat/lon DataArray the rotated-grid tests regrid."""
    return xr.DataArray(
        np.arange(19 * 36, dtype=np.float64).reshape(19, 36),
        dims=("latitude", "longitude"),
        coords={
            "latitude": np.linspace(90.0, -90.0, 19),
            "longitude": np.linspace(0.0, 350.0, 36),
        },
        name="a",
    )


def test_regrid_xarray_geo_dims_match_coords():
    """The derived output dimensions of the geography builder are those of its coordinates."""
    from earthkit.geo.grids.xarray.output import XarrayOutputGeographyBuilder

    # a regular lat/lon grid, a rotated one and a 1D unstructured grid
    specs = [
        ({"grid": [10.0, 10.0], "area": [50, -20, 20, 40]}, ["latitude", "longitude"]),
        (ROT_OUT_GRID, ["y", "x"]),
        ({"latitudes": [10.0, 20.0, 30.0], "longitudes": [0.0, 10.0, 20.0]}, ["values"]),
    ]
    for spec, expected in specs:
        geo = XarrayOutputGeographyBuilder(spec)

        assert geo.geo_dims() == expected

        # what coords() builds its coordinates on is exactly what geo_dims() names
        dims, coords, coords_dim = geo.coords()
        assert sorted(dims) == sorted(expected)
        for name in coords:
            assert set(coords_dim[name]) <= set(geo.geo_dims())


def test_regrid_xarray_geo_dims_computes_latlons_once():
    """geo_dims() and coords() share the one computation of the grid's latitudes/longitudes."""
    from earthkit.geo.grids.xarray.output import XarrayOutputGeographyBuilder

    class _CountingGrid:
        """Delegates to the wrapped grid, counting the to_latlons() calls it sees."""

        def __init__(self, grid):
            self._grid = grid
            self.calls = 0

        def to_latlons(self):
            self.calls += 1
            return self._grid.to_latlons()

        def __getattr__(self, name):
            return getattr(self._grid, name)

    # 2D grids of both kinds: the distinct lat/lon one and the 2D-coordinates one
    for spec in ({"grid": [10.0, 10.0], "area": [50, -20, 20, 40]}, ROT_OUT_GRID):
        geo = XarrayOutputGeographyBuilder(spec)
        counting = _CountingGrid(geo.grid._grid)
        geo.grid._grid = counting

        geo.geo_dims()
        geo.coords()

        assert counting.calls == 1


@requires_mir
def test_regrid_xarray_rotated_out_grid_default_out_dims():
    """A rotated output grid gets its (y, x) geography without an out_dims argument."""
    da = _global_ll_dataarray()

    r = regrid(da, out_grid=ROT_OUT_GRID, interpolation="linear")
    ref = regrid(da, out_grid=ROT_OUT_GRID, interpolation="linear", out_dims=["y", "x"])

    # the result is on the dimensions of its 2D geography, as it was with the out_dims
    # that had to be passed before geo_dims() knew about them; passing them explicitly
    # still wins over the derived ones
    assert list(r.dims) == ["y", "x"]
    compare_dims(r, ROT_OUT_SIZES, sizes=True)

    for name in ("latitude", "longitude"):
        assert r[name].dims == ("y", "x")
        assert r[name].shape == (ROT_OUT_SIZES["y"], ROT_OUT_SIZES["x"])
        np.testing.assert_allclose(r[name].values, ref[name].values)

    xr.testing.assert_allclose(r, ref)


@requires_mir
@requires_dask
def test_regrid_xarray_rotated_out_grid_default_out_dims_lazy(in_process_scheduler):
    """A chunked input regridded onto a rotated grid stays lazy without an out_dims argument."""
    da = _global_ll_dataarray()

    eager = regrid(da, out_grid=ROT_OUT_GRID, interpolation="linear")
    lazy = regrid(da.chunk({"latitude": 5}), out_grid=ROT_OUT_GRID, interpolation="linear")

    assert list(lazy.dims) == ["y", "x"]
    _assert_lazy_like_eager(lazy, eager, ROT_OUT_SIZES)


def test_regrid_xarray_geo_dims_small_rotated_grid_not_regular_mesh():
    """A rotated grid is no regular mesh, not even when its rows and columns hold 2 points."""
    from earthkit.geo.grids.xarray.output import XarrayOutputGeographyBuilder

    # its true latitudes differ along every row, but rows and columns of two points
    # trivially have equal spacing: only a constancy test classifies it correctly
    geo = XarrayOutputGeographyBuilder(ROT_OUT_GRID_2X2)

    assert geo.geo_dims() == ["y", "x"]

    dims, coords, coords_dim = geo.coords()
    assert sorted(dims) == ["x", "y"]
    for name in coords:
        assert coords_dim[name] == ("y", "x")


def test_regrid_xarray_out_dims_follow_backend_grid():
    """Without out_dims the result's dimensions follow the grid the backend returned."""
    ds = _tiny_ll_dataset()

    # the requested grid is a regular mesh, so the regridding runs onto the
    # latitude/longitude dimensions; the stub backend answers with the rotated grid,
    # whose 2D coordinates make the handler rename those dimensions to (y, x)
    r = regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=_recording_backend(ROT_OUT_GRID))
    ref = regrid(
        ds["a"],
        out_grid=STUB_OUT_GRID,
        backend=_recording_backend(ROT_OUT_GRID),
        out_dims=["y", "x"],
    )

    assert list(r.dims) == ["level", "y", "x"]
    compare_dims(r, {"level": 2, "y": 10, "x": 19}, sizes=True)

    for name in ("latitude", "longitude"):
        assert r[name].dims == ("y", "x")
        assert r[name].shape == (ROT_OUT_SIZES["y"], ROT_OUT_SIZES["x"])

    xr.testing.assert_identical(r, ref)


@requires_dask
def test_regrid_xarray_dask_out_dims_follow_backend_grid(in_process_scheduler):
    """The rename onto the backend grid's geography dims applies to a lazy result too."""
    ds = _tiny_ll_dataset()

    eager = regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=_recording_backend(ROT_OUT_GRID))
    lazy = regrid(ds["a"].chunk({"level": 1}), out_grid=STUB_OUT_GRID, backend=_recording_backend(ROT_OUT_GRID))

    assert list(lazy.dims) == ["level", "y", "x"]
    _assert_lazy_like_eager(lazy, eager, {"level": 2, "y": 10, "x": 19})


def test_regrid_xarray_out_dims_follow_backend_grid_clash():
    """A rename onto dimensions the result already has is refused with a clear error."""
    # a batch dimension named "y" clashes with the (y, x) geography of the stub's grid
    ds = _tiny_ll_dataset().rename({"level": "y"})

    with pytest.raises(ValueError, match="out_dims"):
        regrid(ds["a"], out_grid=STUB_OUT_GRID, backend=_recording_backend(ROT_OUT_GRID))


@requires_mir
def test_regrid_xarray_rotated_small_grid_default_out_dims():
    """The review reproduction: a rotated 2 x 2 request regrids fine without out_dims."""
    da = _global_ll_dataarray()

    r = regrid(da, out_grid=ROT_OUT_GRID_2X2, interpolation="linear")
    ref = regrid(da, out_grid=ROT_OUT_GRID_2X2, interpolation="linear", out_dims=["y", "x"])

    # neither the requested 2 x 2 mesh nor the 3 x 2 grid MIR returns from it is a
    # regular lat/lon mesh: the geography is 2D on (y, x) whatever the classification
    # of the request got wrong before
    assert list(r.dims) == ["y", "x"]
    compare_dims(r, {"y": 3, "x": 2}, sizes=True)

    for name in ("latitude", "longitude"):
        assert r[name].dims == ("y", "x")
        assert r[name].shape == (3, 2)
        np.testing.assert_allclose(r[name].values, ref[name].values)

    xr.testing.assert_allclose(r, ref)


@requires_mir
@requires_dask
def test_regrid_xarray_rotated_small_grid_default_out_dims_lazy(in_process_scheduler):
    """The chunked twin of the review reproduction stays lazy on the (y, x) grid."""
    da = _global_ll_dataarray()

    eager = regrid(da, out_grid=ROT_OUT_GRID_2X2, interpolation="linear")
    lazy = regrid(da.chunk({"latitude": 5}), out_grid=ROT_OUT_GRID_2X2, interpolation="linear")

    assert list(lazy.dims) == ["y", "x"]
    _assert_lazy_like_eager(lazy, eager, {"y": 3, "x": 2})
