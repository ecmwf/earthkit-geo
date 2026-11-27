# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest

from earthkit.geo.grids.xarray.loader import variables as get_variables
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
    # the points form a complete regular lat/lon grid, so the grid is now
    # auto-detected as "regular_ll" rather than treated as "unstructured_ll"
    assert ek_grid.type == "regular_ll"
    assert ek_grid.shape == (7, 12)

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
@pytest.mark.parametrize(
    "filename,order,spec_order",
    [
        # 0.2 degree ocean grid excluding the poles, with the latitudes
        # ascending (-89.8 to 89.8) and descending (89.8 to -89.8)
        ("latlon_02x02_ocean_nopole_jplus.nc", "i+j+", "i+j+"),
        ("latlon_02x02_ocean_nopole_jminus.nc", "i+j-", None),
    ],
)
def test_xarray_grid_netcdf_ll_nopole(filename, order, spec_order):
    path = get_test_data(filename, subfolder="xr")
    ds = xr.open_dataset(path)

    lat_num, lon_num = 899, 1800
    variable_grid_dims = ("latitude", "longitude")
    points_num = lat_num * lon_num

    # auto grid-spec detection
    variables = get_variables(ds)

    assert len(variables) == 1
    v = variables[0]

    # the points form a complete regular lat/lon grid, so the grid is
    # detected as "regular_ll" rather than treated as "unstructured_ll"
    ek_grid = v.ek_grid
    assert ek_grid.type == "regular_ll"
    assert ek_grid.order == order
    assert ek_grid.shape == (lat_num, lon_num)

    # the detected gridspec: the increment is the average of the float32 steps
    # and the area is the corner values of the coordinate arrays
    spec = v.xr_grid.grid_spec
    assert spec["grid"] == pytest.approx([0.2, 0.2])
    assert spec["area"] == pytest.approx([89.8, -180.0, -89.8, 179.8])

    if spec_order is None:
        # the default, north to south latitude order is not repeated in the spec
        assert "order" not in spec
    else:
        assert spec["order"] == spec_order

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_xarray_grid_netcdf_ll_twelfths():
    # a 1/12 degree float32 ocean grid: its increment has more than 6 significant
    # digits, so it is detected by snapping the axis onto the 1/12 lattice rather
    # than by rounding it (see _snapped_axis)
    path = get_test_data("latlon_1_per_12_ocean_nopole_jplus.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    lat_num, lon_num = 2041, 4320
    variable_grid_dims = ("latitude", "longitude")
    points_num = lat_num * lon_num

    # auto grid-spec detection
    variables = get_variables(ds)

    assert len(variables) == 1
    v = variables[0]

    ek_grid = v.ek_grid
    assert ek_grid.type == "regular_ll"
    assert ek_grid.order == "i+j+"
    assert ek_grid.shape == (lat_num, lon_num)

    # the increment is carried at the 12 decimals of the gridspec, not the
    # 6 significant digits of the raw float32 cleaning
    spec = v.xr_grid.grid_spec
    assert spec["grid"] == pytest.approx([0.083333333333, 0.083333333333], abs=1e-12)
    assert spec["area"] == pytest.approx([90.0, -180.0, -80.0, 179.916666666667], abs=1e-9)
    assert spec["order"] == "i+j+"
    assert "reference" not in spec

    # the grid built from the spec reproduces the file's coordinates: the file
    # stepped them in float32, so they sit up to a couple of ulps off the
    # 1/12 lattice (an ulp at 180 degrees is about 1.5e-5)
    assert np.allclose(ek_grid.distinct_latitudes(), ds["latitude"].values, rtol=0, atol=1e-4)
    assert np.allclose(ek_grid.distinct_longitudes(), ds["longitude"].values, rtol=0, atol=1e-4)

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

    # the CF lambert_azimuthal_equal_area grid mapping and the regular 1D x/y
    # axes, in metres, are recognised, so the grid is detected as a regular grid
    # in the plane of that projection rather than treated as "unstructured_ll"
    projection_source = ds["lambert_azimuthal_equal_area"].attrs["proj4_params"]

    for variable in variables:
        assert variable.ek_grid.type == "regular_xy"
        assert variable.ek_grid.shape == (950, 1000)

        spec = variable.xr_grid.grid_spec
        assert spec["type"] == "regular_xy"
        # the ranges are the file's own axes, [start, stop, step], y descending
        assert spec["x"] == [2502500.0, 7497500.0, 5000.0]
        assert spec["y"] == [5497500.0, 752500.0, -5000.0]
        assert spec["projection"] == {"type": "proj", "source": projection_source}

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)
    assert np.allclose(lat, ds["latitude"].values.flatten())
    assert np.allclose(lon, ds["longitude"].values.flatten())

    # the eckit grid built from the spec reproduces the file's own coordinates,
    # in the same (row major, x fastest) point order as the latlons above
    ek_lat, ek_lon = v.ek_grid.to_latlons()
    dlat = np.max(np.abs(np.asarray(ek_lat) - ds["latitude"].values.flatten()))
    dlon = np.max(np.abs((np.asarray(ek_lon) - ds["longitude"].values.flatten() + 180.0) % 360.0 - 180.0))
    assert dlat <= 1e-3
    assert dlon <= 1e-3


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_xarray_grid_polar_stereo_north():
    # the AMSR2 polar stereographic sample: a projected grid held as its two 1D
    # axes only, with no 2D latitude/longitude field to check a gridspec against,
    # and a grid mapping describing the projection by its crs_wkt alone
    path = get_test_data("polar_stereo_north_amsr2.nc", subfolder="xr")
    ds = xr.open_dataset(path)

    variable_grid_dims = ("y", "x")
    points_num = 1792 * 1216

    variables = get_variables(ds)

    assert len(variables) == 1
    v = variables[0]
    assert v.name == "z"

    # the axes and the grid mapping are enough: the grid is detected as a regular
    # grid in the plane of that projection, not treated as "unstructured_ll"
    assert v.ek_grid.type == "regular_xy"
    assert v.ek_grid.shape == (1792, 1216)

    spec = v.xr_grid.grid_spec
    assert spec["type"] == "regular_xy"
    # the ranges are the file's own axes, [start, stop, step]. Its y axis runs
    # north, so the spec's does too and the rows of the data keep their order
    assert spec["x"] == [-3846875.0, 3746875.0, 6250.0]
    assert spec["y"] == [-5346875.0, 5846875.0, 6250.0]

    from earthkit.geo.grids.xarray.grid_spec import proj_projection

    projection = proj_projection(dict(ds["polar_stereographic"].attrs))
    assert projection["type"] == "proj"
    assert "+proj=stere" in projection["source"]
    assert spec["projection"] == projection

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims
    assert xr_grid.grid_type == "regular_xy"
    assert list(v.geo_dims) == list(variable_grid_dims)

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)

    # the eckit grid built from the spec reproduces the points the projection
    # computes for the file's axes, in the same (row major, x fastest) order
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in v.ek_grid.to_latlons())
    assert np.max(np.abs(ek_lat - lat)) <= 1e-3
    dlon = np.max(np.abs((ek_lon - lon + 180.0) % 360.0 - 180.0))
    assert dlon <= 1e-3


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_xarray_grid_reduced_gaussian_O1280():
    # a real IFS sample on an O1280 octahedral reduced Gaussian grid: the data
    # variable carries no latitudes or longitudes of its own, only the
    # "reduced_gaussian" grid mapping, so the grid has to be detected from that
    path = get_test_data("t2m_O1280.nc", subfolder="xr")

    with xr.open_dataset(path) as ds:
        points_num = 6599680

        variables = get_variables(ds)

        # only the data variable is returned: the "pl" Gaussian latitude coordinate
        # and the scalar "reduced_gaussian" grid mapping are not
        assert len(variables) == 1
        v = variables[0]
        assert v.name == "air_temperature"

        assert v.xr_grid.grid_spec == {"grid": "O1280"}

        assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
        assert v.geo_dims == ("reduced_gaussian_index",)

        assert v.ek_grid.type == "reduced_gg"
        assert v.ek_grid.shape == (points_num,)

        lat, lon = v.xr_grid.latlons
        assert lat.shape == (points_num,)
        assert lon.shape == (points_num,)

        # the northernmost Gaussian parallel, held as its full row of 20 points
        assert lat[:20] == pytest.approx(89.946, abs=1e-3)
        assert lat[20] < lat[0]
        assert lat[-1] == pytest.approx(-89.946, abs=1e-3)
        assert lon[:20] == pytest.approx(18.0 * np.arange(20), abs=1e-3)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_xarray_grid_reduced_gaussian_O32():
    # a small synthetic file on an octahedral (grid_subtype "octahedral") reduced
    # Gaussian grid, laid out like the real O1280 sample: the grid has to be
    # detected from the "reduced_gaussian" grid mapping alone
    path = get_test_data("t2m_O32.nc", subfolder="xr")

    points_num = 5248

    with xr.open_dataset(path) as ds:
        variables = get_variables(ds)

        # only the data variable is returned: the "pl" table and the scalar
        # "reduced_gaussian" grid mapping are not
        assert [v.name for v in variables] == ["air_temperature"]

        v = variables[0]
        assert v.xr_grid.grid_spec == {"grid": "O32"}
        assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
        assert v.geo_dims == ("reduced_gaussian_index",)

        assert v.ek_grid.type == "reduced_gg"
        assert v.ek_grid.shape == (points_num,)

        lat, lon = v.xr_grid.latlons
        assert lat.shape == (points_num,)
        assert lon.shape == (points_num,)

        # the northernmost Gaussian parallel, held as its full row of 20 points
        assert np.allclose(lat[:20], lat[0])
        assert lat[20] < lat[0]
        assert np.allclose(lon[:20], 18.0 * np.arange(20))

        # the points are the ones of the eckit grid built from the same spec
        ek_lat, ek_lon = (np.asarray(a, dtype=np.float64) for a in v.ek_grid.to_latlons())
        assert np.array_equal(lat, ek_lat)
        assert np.array_equal(lon, ek_lon)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
@pytest.mark.parametrize("n, points_num", [(16, 1688), (80, 35718)])
def test_xarray_grid_reduced_gaussian_normal_netcdf(n, points_num):
    # a synthetic file on a classical (grid_subtype "normal") reduced Gaussian
    # grid, written the way the real O1280 sample is: its attributes, its integer
    # table of points per latitude and its scalar grid mapping come through the
    # netCDF encoding, and the grid has to be detected from them
    path = get_test_data(f"t2m_N{n}.nc", subfolder="xr")

    with xr.open_dataset(path) as ds:
        variables = get_variables(ds)

        # the grid mapping and the point table are skipped, as before: only the
        # data variable is left, on the grid the mapping describes
        assert [v.name for v in variables] == ["air_temperature"]

        v = variables[0]
        assert v.xr_grid.grid_spec == {"grid": f"N{n}"}
        assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
        assert v.geo_dims == ("reduced_gaussian_index",)

        assert v.ek_grid.type == "reduced_gg"
        assert v.ek_grid.shape == (points_num,)

        lat, lon = v.xr_grid.latlons
        assert lat.shape == (points_num,)
        assert lon.shape == (points_num,)

        # the northernmost Gaussian parallel, held as its full row of points: one
        # latitude, longitudes spaced all the way round the globe from Greenwich
        first_row = int(ds["pl"].values[0])
        assert np.allclose(lat[:first_row], lat[0])
        assert np.allclose(lon[:first_row], 360.0 / first_row * np.arange(first_row))

        # the points are the ones of the eckit grid built from the same spec
        ek_lat, ek_lon = (np.asarray(a, dtype=np.float64) for a in v.ek_grid.to_latlons())
        assert np.array_equal(lat, ek_lat)
        assert np.array_equal(lon, ek_lon)


# the twelve CF reduced Gaussian sub-area files of the download server, cut by
# tools/test-data/make_cf_reduced_gaussian_files.py from the real t2m_O1280 sample
# and from the fixed global t2m_O32 and t2m_N80 files: the file, the global grid it
# cuts, the case it holds, the number of points it keeps, and the [north, west,
# south, east] area those points sit in -- None for the selective subset, which is
# a scatter of points no area of the grid describes
T2M_SUBAREA_FILES = [
    pytest.param("t2m_O1280_area_east.nc", "O1280", "east", 10184, [50.0, 10.0, 40.0, 20.0], id="O1280-area-east"),
    pytest.param("t2m_O1280_area_cross0.nc", "O1280", "cross0", 9058, [55.0, -5.0, 45.0, 5.0], id="O1280-area-cross0"),
    pytest.param("t2m_O1280_band.nc", "O1280", "band", 10344, [45.0, 0.0, 44.7, 360.0], id="O1280-band"),
    pytest.param("t2m_O1280_selective.nc", "O1280", "selective", 3101, None, id="O1280-selective"),
    pytest.param("t2m_O32_area_east.nc", "O32", "east", 75, [60.0, 10.0, 30.0, 40.0], id="O32-area-east"),
    pytest.param("t2m_O32_area_cross0.nc", "O32", "cross0", 123, [60.0, -10.0, 30.0, 40.0], id="O32-area-cross0"),
    pytest.param("t2m_O32_band.nc", "O32", "band", 728, [40.0, 0.0, 20.0, 360.0], id="O32-band"),
    pytest.param("t2m_O32_selective.nc", "O32", "selective", 852, None, id="O32-selective"),
    pytest.param("t2m_N80_area_east.nc", "N80", "east", 571, [60.0, 10.0, 30.0, 40.0], id="N80-area-east"),
    pytest.param("t2m_N80_area_cross0.nc", "N80", "cross0", 947, [60.0, -10.0, 30.0, 40.0], id="N80-area-cross0"),
    pytest.param("t2m_N80_band.nc", "N80", "band", 5424, [40.0, 0.0, 20.0, 360.0], id="N80-band"),
    pytest.param("t2m_N80_selective.nc", "N80", "selective", 5155, None, id="N80-selective"),
]

# the tolerance of the checks against the area a cut was made with: its points are
# the exact ones of the global grid, so only their float representation can differ
SUBAREA_TOL = 1e-6


def _assert_in_area(lat, lon, area, tol=SUBAREA_TOL):
    """Check flat points against the [north, west, south, east] area they were cut with.

    Their latitudes are bounded by those of the area. Their longitudes are counted
    from its western edge, modulo 360, so an area which crosses the Greenwich
    meridian -- and whose points the file numbers from that meridian round the
    globe -- holds the ones east of it and the ones west of it alike.
    """
    north, west, south, east = area
    lat, lon = np.asarray(lat, dtype=np.float64), np.asarray(lon, dtype=np.float64)

    assert lat.min() >= south - tol and lat.max() <= north + tol, (float(lat.min()), float(lat.max()))
    offset = (lon - west) % 360.0
    assert float(offset.max()) <= east - west + tol, (float(offset.max()), east - west)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
@pytest.mark.parametrize("filename,grid,kind,points,area", T2M_SUBAREA_FILES)
def test_xarray_grid_reduced_gaussian_subarea_files(filename, grid, kind, points, area):
    # a cut of a CF reduced Gaussian file keeps the global mapping and table of its
    # source and holds a strictly increasing subset of the global index: the reader
    # finds the area of the grid those points are, or, for a subset no area describes,
    # the unstructured scatter of points the file holds
    path = get_test_data(filename, subfolder="xr")

    with xr.open_dataset(path) as ds:
        variables = get_variables(ds)
        values = ds["air_temperature"].values

        assert [v.name for v in variables] == ["air_temperature"]
        (v,) = variables

        assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
        assert v.geo_dims == ("reduced_gaussian_index",)
        assert v.ek_grid.shape == (points,)
        assert values.shape == (1, 1, points)

        # the file holds the temperatures of its source, within the range the
        # generator checked them against, and not the zeros some test files used to
        assert np.all(np.isfinite(values))
        assert float(values.min()) >= 200.0 and float(values.max()) <= 330.0

        lat, lon = v.xr_grid.latlons
        spec = v.xr_grid.grid_spec
        order = v.xr_grid.point_order

        if kind == "selective":
            # the points are a scatter, and the grid is built from them
            assert spec is None
            assert v.xr_grid.grid_type == "unstructured"
            assert v.ek_grid.type == "unstructured_ll"
            assert order is None

            # the file's points are the global grid's points its CF index names:
            # the same latitudes, and the same longitudes modulo 360
            index = np.asarray(ds["reduced_gaussian_index"].values, dtype=np.int64)
            assert index.size == points
            assert np.array_equal(index, np.unique(index))
            _assert_points_at_index(grid, index, lat, lon)
        else:
            assert v.xr_grid.grid_type == "reduced_gg"
            assert v.ek_grid.type == "reduced_gg"
            assert set(spec) == {"grid", "area"}
            assert spec["grid"] == grid

            north, west, south, east = spec["area"]
            assert south <= north
            # the area the detection names is the one the points were cut with, its
            # latitudes the outermost lines the file holds, and a band of full lines
            # the whole longitude circle from 0 to 360
            if kind == "band":
                assert west == pytest.approx(0.0, abs=SUBAREA_TOL)
                assert east == pytest.approx(360.0, abs=SUBAREA_TOL)
            else:
                assert west % 360.0 == pytest.approx(area[1] % 360.0, abs=SUBAREA_TOL)
                assert east % 360.0 == pytest.approx(area[3] % 360.0, abs=SUBAREA_TOL)
            assert south >= area[2] - SUBAREA_TOL and north <= area[0] + SUBAREA_TOL

            # an area crossing the Greenwich meridian is not the order eckit iterates
            # its points in, which the grid carries as a permutation; the others are
            if kind == "cross0":
                assert order is not None
                assert not np.array_equal(order, np.arange(order.size))
            else:
                assert order is None

            _assert_in_area(lat, lon, area)

            # the eckit grid of the spec holds the file's points, and its own point p
            # sits at position order[p] of the file
            ek_lat, ek_lon = _eckit_points(v.ek_grid)
            index = np.arange(points) if order is None else np.asarray(order)
            assert np.array_equal(np.asarray(lat)[index], ek_lat)
            assert np.array_equal(np.asarray(lon)[index], ek_lon)


@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
@pytest.mark.parametrize("grid,points", [("O32", 5248), ("N16", 1688), ("N80", 35718)])
def test_xarray_grid_reduced_gaussian_t2m_global_values(grid, points):
    # the three small global files of the server used to hold all-zero values: the
    # fixed ones hold the temperature of the t2m_O1280 sample regridded onto their
    # grid, which they are still detected as, mapping and table unchanged
    path = get_test_data(f"t2m_{grid}.nc", subfolder="xr")

    with xr.open_dataset(path) as ds:
        variables = get_variables(ds)
        values = ds["air_temperature"].values

        assert [v.name for v in variables] == ["air_temperature"]
        (v,) = variables

        assert v.xr_grid.grid_spec == {"grid": grid}
        assert v.ek_grid.type == "reduced_gg"
        assert v.ek_grid.shape == (points,)
        assert values.shape == (1, 1, points)

        assert np.all(np.isfinite(values))
        assert int(np.count_nonzero(values)) == values.size  # no point left at zero
        assert float(values.min()) >= 200.0 and float(values.max()) <= 330.0


# a Lambert azimuthal equal-area projection on GRS80, with the false easting and
# northing of the efas sample, described the way CF files describe a projection
LAEA_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


def _projected_dataset(grid_mapping=True):
    """A small synthetic dataset on the plane of a LAEA projection.

    Parameters
    ----------
    grid_mapping : bool
        Whether the variable refers to its projection by a CF ``grid_mapping``
        attribute naming the scalar mapping variable held in the dataset.

    Returns
    -------
    xr.Dataset
        The dataset, with its data variable "a" over the 1D ``x``/``y`` axes of
        the projection plane (5 km, ``y`` descending) and the 2D
        ``latitude``/``longitude`` coordinates pyproj computes for them.
    """
    from pyproj import CRS, Transformer

    x_values = 4321000.0 + 5000.0 * np.arange(6)
    y_values = 3210000.0 - 5000.0 * np.arange(4)

    transformer = Transformer.from_crs(CRS.from_proj4(LAEA_PROJ4), "EPSG:4326", always_xy=True)
    xv, yv = np.meshgrid(x_values, y_values)
    lon2d, lat2d = transformer.transform(xv, yv)

    attrs = {"grid_mapping": "mapping"} if grid_mapping else {}
    ds = xr.Dataset(
        {"a": xr.Variable(["y", "x"], np.arange(24, dtype=np.float64).reshape(4, 6), attrs=attrs)},
        coords={
            "x": (["x"], x_values, {"standard_name": "projection_x_coordinate", "units": "m"}),
            "y": (["y"], y_values, {"standard_name": "projection_y_coordinate", "units": "m"}),
            "latitude": (["y", "x"], lat2d, {"standard_name": "latitude", "units": "degrees_north"}),
            "longitude": (["y", "x"], lon2d, {"standard_name": "longitude", "units": "degrees_east"}),
        },
    )
    if grid_mapping:
        ds["mapping"] = xr.DataArray(
            0.0,
            attrs={"grid_mapping_name": "lambert_azimuthal_equal_area", "proj4_params": LAEA_PROJ4},
        )
    return ds


def test_xarray_grid_projected_xy():
    # the projected grid of the efas sample, at a size a test can build on the fly
    ds = _projected_dataset()

    variables = get_variables(ds)
    assert len(variables) == 1

    v = variables[0]
    assert v.name == "a"

    # the CF grid mapping names a projection PROJ knows, and the 1D x/y axes are
    # regular: the grid is a regular_xy one in the plane of that projection
    ek_grid = v.ek_grid
    assert ek_grid.type == "regular_xy"
    assert ek_grid.shape == (4, 6)

    spec = v.xr_grid.grid_spec
    assert spec["x"] == pytest.approx([4321000.0, 4346000.0, 5000.0])
    assert spec["y"] == pytest.approx([3210000.0, 3195000.0, -5000.0])
    assert spec["projection"] == {"type": "proj", "source": LAEA_PROJ4}

    assert v.xr_grid.variable_dims == ("y", "x")

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (24,)
    assert lon.shape == (24,)
    assert np.allclose(lat, ds["latitude"].values.flatten())

    # the eckit grid built from the spec reproduces the dataset's own coordinates
    ek_lat, ek_lon = ek_grid.to_latlons()
    dlat = np.max(np.abs(np.asarray(ek_lat) - ds["latitude"].values.flatten()))
    dlon = np.max(np.abs((np.asarray(ek_lon) - ds["longitude"].values.flatten() + 180.0) % 360.0 - 180.0))
    assert dlat <= 1e-3
    assert dlon <= 1e-3


def test_xarray_grid_projected_without_grid_mapping():
    # the same coordinates, with no CF grid mapping to say what projection they
    # are in: nothing is guessed, and the grid is still the unstructured one
    ds = _projected_dataset(grid_mapping=False)

    variables = get_variables(ds)
    assert len(variables) == 1

    v = variables[0]
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (24,)

    assert v.xr_grid.grid_spec is None
    assert v.xr_grid.variable_dims == ("y", "x")

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (24,)
    assert lon.shape == (24,)
    assert np.allclose(lat, ds["latitude"].values.flatten())


def _projected_pair_dataset(mapping_first: bool = True) -> xr.Dataset:
    """A synthetic projected dataset with two variables sharing one lat/lon pair.

    Parameters
    ----------
    mapping_first : bool
        Whether the variable referring to the grid mapping comes first among the
        dataset's data variables. The two are otherwise the same.

    Returns
    -------
    xr.Dataset
        The dataset, with the data variables "mapped", carrying a CF
        ``grid_mapping`` attribute naming the scalar mapping variable, and
        "unmapped", without one, both defined on the very same 2D
        ``latitude``/``longitude`` coordinates of the projection plane.
    """
    ds = _projected_dataset().rename({"a": "mapped"})
    ds["unmapped"] = xr.Variable(["y", "x"], np.zeros((4, 6), dtype=np.float64))

    order = ["mapped", "unmapped"] if mapping_first else ["unmapped", "mapped"]
    return ds[order + ["mapping"]]


@pytest.mark.parametrize("mapping_first", [False, True])
def test_xarray_grid_projected_grid_mapping_per_variable(mapping_first):
    # which of two variables sharing the same coordinates carries the grid mapping
    # must not depend on the order they are declared in: the guesser caches the
    # grid it derives for a coordinate pair, and that grid now depends on the
    # variable's mapping, so the cache keeps the mappings apart
    ds = _projected_pair_dataset(mapping_first)

    variables = {v.name: v for v in get_variables(ds)}
    assert sorted(variables) == ["mapped", "unmapped"]

    assert variables["mapped"].ek_grid.type == "regular_xy"
    assert variables["mapped"].xr_grid.grid_spec["projection"] == {"type": "proj", "source": LAEA_PROJ4}

    assert variables["unmapped"].ek_grid.type == "unstructured_ll"
    assert variables["unmapped"].ek_grid.shape == (24,)
    assert variables["unmapped"].xr_grid.grid_spec is None


def _plane_latlons(x_values, y_values):
    """The ``(lat, lon)`` fields the projection computes for the mesh of two plane axes.

    Parameters
    ----------
    x_values : Any
        The x axis of the plane, in metres.
    y_values : Any
        The y axis of the plane, in metres.

    Returns
    -------
    Tuple[Any, Any]
        The latitudes and longitudes, shaped ``(y, x)``: the x axis fastest, the
        order a variable indexed ``(y, x)`` is flattened into.
    """
    from pyproj import CRS, Transformer

    transformer = Transformer.from_crs(CRS.from_proj4(LAEA_PROJ4), "EPSG:4326", always_xy=True)
    xv, yv = np.meshgrid(np.asarray(x_values, dtype=np.float64), np.asarray(y_values, dtype=np.float64))
    lon, lat = transformer.transform(xv, yv)

    return lat, lon


def _projected_axes_only_dataset(
    ascending_y: bool = False, dims=("y", "x"), mapping: str = "proj", irregular_x: bool = False
):
    """A small synthetic dataset held on the two 1D axes of a projection plane.

    Unlike :func:`_projected_dataset` it carries no 2D ``latitude``/``longitude``
    coordinates: the grid can only come from the axes and the CF grid mapping, as
    in the polar stereographic AMSR2 sample.

    Parameters
    ----------
    ascending_y : bool
        Whether the y axis runs north, from its lowest value up.
    dims : tuple
        The dimension order of the data variable.
    mapping : str
        ``"proj"`` for a grid mapping carrying the projection as a CF ``crs_wkt``
        and ``proj4_params``, or ``"cf"`` for one carrying the individual
        projection parameters only, which is no PROJ description at all.
    irregular_x : bool
        Whether to move one of the x values off the step of the axis.

    Returns
    -------
    xr.Dataset
        The dataset: the regular 5 km ``x``/``y`` axes, the scalar grid mapping
        variable and the data variable "a" indexed by them.
    """
    from pyproj import CRS

    x_values = 4321000.0 + 5000.0 * np.arange(6)
    if irregular_x:
        x_values = x_values.copy()
        x_values[3] += 1000.0

    y_values = 3210000.0 - 5000.0 * np.arange(4)
    if ascending_y:
        y_values = y_values[::-1]

    if mapping == "proj":
        attrs = {
            "grid_mapping_name": "lambert_azimuthal_equal_area",
            "proj4_params": LAEA_PROJ4,
            "crs_wkt": CRS.from_proj4(LAEA_PROJ4).to_wkt(),
        }
    elif mapping == "cf":
        attrs = {
            "grid_mapping_name": "lambert_azimuthal_equal_area",
            "latitude_of_projection_origin": 52.0,
            "longitude_of_projection_origin": 10.0,
            "false_easting": 4321000.0,
            "false_northing": 3210000.0,
        }
    else:
        raise ValueError(f"unknown mapping {mapping}")

    values = np.arange(24, dtype=np.float64).reshape(4, 6)
    if dims == ("x", "y"):
        values = values.T

    ds = xr.Dataset(
        {"a": xr.Variable(list(dims), values, attrs={"grid_mapping": "mapping"})},
        coords={
            "x": (["x"], x_values, {"standard_name": "projection_x_coordinate", "units": "m"}),
            "y": (["y"], y_values, {"standard_name": "projection_y_coordinate", "units": "m"}),
        },
    )
    ds["mapping"] = xr.DataArray(0.0, attrs=attrs)
    return ds


@pytest.mark.parametrize("ascending_y", [False, True])
def test_xarray_grid_projected_axes_only(ascending_y):
    # a projected file holding nothing but its two 1D axes: the axes and the grid
    # mapping are enough to detect the regular_xy grid, whichever way y runs
    ds = _projected_axes_only_dataset(ascending_y=ascending_y)

    v = get_variables(ds)[0]

    ek_grid = v.ek_grid
    assert ek_grid.type == "regular_xy"
    assert ek_grid.shape == (4, 6)

    spec = v.xr_grid.grid_spec
    assert spec["type"] == "regular_xy"
    assert spec["x"] == pytest.approx([4321000.0, 4346000.0, 5000.0])
    # the y range runs in the direction of the axis, so the data keeps its rows
    y_low, y_high = 3195000.0, 3210000.0
    assert spec["y"] == pytest.approx([y_low, y_high, 5000.0] if ascending_y else [y_high, y_low, -5000.0])
    assert spec["projection"] == {"type": "proj", "source": LAEA_PROJ4}

    assert v.xr_grid.variable_dims == ("y", "x")
    assert v.xr_grid.grid_type == "regular_xy"
    assert list(v.geo_dims) == ["y", "x"]

    lat, lon = v.xr_grid.latlons
    ref_lat, ref_lon = (a.reshape(-1) for a in _plane_latlons(ds["x"].values, ds["y"].values))
    assert lat.shape == (24,)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-3
    assert np.max(np.abs(lon - ref_lon)) <= 1e-3

    # the eckit grid built from the spec holds the same points, in the same order
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in ek_grid.to_latlons())
    assert np.max(np.abs(ek_lat - lat)) <= 1e-3
    assert np.max(np.abs((ek_lon - lon + 180.0) % 360.0 - 180.0)) <= 1e-3


def test_xarray_grid_projected_axes_only_x_major():
    # the data varies fastest along y: eckit scans a regular_xy grid with its x
    # axis fastest and cannot be told otherwise, so the spec is rejected by the
    # point check and the grid falls back to unstructured points, laid out in the
    # order the data holds them
    ds = _projected_axes_only_dataset(dims=("x", "y"))

    v = get_variables(ds)[0]

    assert v.xr_grid.variable_dims == ("x", "y")
    assert v.xr_grid.grid_spec is None
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (24,)
    assert list(v.geo_dims) == ["x", "y"]

    lat, lon = v.xr_grid.latlons
    ref_lat, ref_lon = (a.T.reshape(-1) for a in _plane_latlons(ds["x"].values, ds["y"].values))
    assert np.max(np.abs(lat - ref_lat)) <= 1e-3
    assert np.max(np.abs(lon - ref_lon)) <= 1e-3

    # the points of the unstructured grid are those of the mesh, y fastest, so the
    # values of the (x, y) data flattened that way match them
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in v.ek_grid.to_latlons())
    assert np.max(np.abs(ek_lat - lat)) <= 1e-9
    assert np.max(np.abs(ek_lon - lon)) <= 1e-9


@pytest.mark.parametrize("dims", [("y", "x"), ("x", "y")])
def test_xarray_grid_projected_axes_only_renamed_dims(dims):
    # the variable is indexed by dimensions renamed away from the names of the
    # coordinates carrying them (rename_dims leaves those under their own
    # names): the point order must be read off the axes' dimensions, never
    # their names, so the points keep the order the data is flattened in
    ds = _projected_axes_only_dataset(dims=dims).rename_dims({"x": "i", "y": "j"})

    v = get_variables(ds)[0]

    # the x axis sits on the dimension i, the y axis on j
    data_dims = ("j", "i") if dims == ("y", "x") else ("i", "j")
    assert v.xr_grid.variable_dims == data_dims
    assert list(v.geo_dims) == list(data_dims)

    ref = _plane_latlons(ds["x"].values, ds["y"].values)  # the plane mesh, shaped (y, x)
    if dims == ("x", "y"):
        # the data varies fastest along y: no regular_xy spec can hold it, and
        # the fallback points follow the data's own (i, j) flattening
        assert v.xr_grid.grid_spec is None
        assert v.ek_grid.type == "unstructured_ll"
        assert v.ek_grid.shape == (24,)
        ref_lat, ref_lon = (a.T.reshape(-1) for a in ref)
    else:
        assert v.xr_grid.grid_spec["type"] == "regular_xy"
        assert v.ek_grid.type == "regular_xy"
        assert v.ek_grid.shape == (4, 6)
        ref_lat, ref_lon = (a.reshape(-1) for a in ref)

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (24,)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-3
    assert np.max(np.abs((lon - ref_lon + 180.0) % 360.0 - 180.0)) <= 1e-3

    # the grid which takes the data's values holds the same points, in the same order
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in v.ek_grid.to_latlons())
    assert np.max(np.abs(ek_lat - lat)) <= 1e-3
    assert np.max(np.abs((ek_lon - lon + 180.0) % 360.0 - 180.0)) <= 1e-3


def test_xarray_grid_projected_axes_only_irregular_x():
    # an axis which does not step regularly: no gridspec, and the grid is the
    # unstructured points of the mesh, as before
    ds = _projected_axes_only_dataset(irregular_x=True)

    v = get_variables(ds)[0]

    assert v.xr_grid.variable_dims == ("y", "x")
    assert v.xr_grid.grid_spec is None
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (24,)

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (24,)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))


def test_xarray_grid_projected_axes_only_without_proj_string():
    # a grid mapping carrying its projection in the CF parameters alone: pyproj
    # computes the points from them, but there is no PROJ description to build a
    # gridspec from, so the grid keeps its unstructured treatment
    ds = _projected_axes_only_dataset(mapping="cf")

    v = get_variables(ds)[0]

    assert v.xr_grid.grid_spec is None
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (24,)

    lat, lon = v.xr_grid.latlons
    ref_lat, ref_lon = (a.reshape(-1) for a in _plane_latlons(ds["x"].values, ds["y"].values))
    assert np.max(np.abs(lat - ref_lat)) <= 1e-3
    assert np.max(np.abs((lon - ref_lon + 180.0) % 360.0 - 180.0)) <= 1e-3


# the plane coordinates of the scattered points of _projected_scatter_dataset:
# pairs of a LAEA plane, one point per index, as a file holding its irregular
# points or stations in a single dimension spells them
SCATTER_X = [4321000.0, 4326000.0, 4331000.0, 4321000.0, 4341000.0, 4351000.0]
SCATTER_Y = [3210000.0, 3205000.0, 3200000.0, 3190000.0, 3215000.0, 3195000.0]


def _projected_scatter_dataset() -> xr.Dataset:
    """A small synthetic dataset of projected points paired along one dimension.

    Unlike :func:`_projected_axes_only_dataset` its ``x``/``y`` are not the two
    axes of a mesh: they are two 1D coordinates over the *same* dimension, so
    index ``i`` holds one point of the plane, and the grid is a scatter of those
    points. The mapping carries its projection as a CF ``crs_wkt``, which is
    what pyproj reads, and there are no 2D ``latitude``/``longitude`` to fall
    back on.

    Returns
    -------
    xr.Dataset
        The dataset: the paired 1D ``x``/``y`` coordinates over ``values``, the
        data variable "a" indexed by it and the scalar grid mapping variable
        "mapping".
    """
    from pyproj import CRS

    ds = xr.Dataset(
        {"a": xr.Variable(["values"], np.arange(6, dtype=np.float64), attrs={"grid_mapping": "mapping"})},
        coords={
            "x": (["values"], np.asarray(SCATTER_X), {"standard_name": "projection_x_coordinate", "units": "m"}),
            "y": (["values"], np.asarray(SCATTER_Y), {"standard_name": "projection_y_coordinate", "units": "m"}),
        },
    )
    ds["mapping"] = xr.DataArray(
        0.0,
        attrs={"grid_mapping_name": "lambert_azimuthal_equal_area", "crs_wkt": CRS.from_proj4(LAEA_PROJ4).to_wkt()},
    )
    return ds


def test_xarray_grid_projected_scatter_points():
    # the projected points of a variable which is not a mesh: the pairs of its
    # x/y coordinates, transformed with the projection, and no gridspec to claim
    ds = _projected_scatter_dataset()

    variables = get_variables(ds)
    assert len(variables) == 1

    v = variables[0]
    assert v.name == "a"

    assert v.xr_grid.variable_dims == ("values",)
    assert list(v.geo_dims) == ["values"]

    # a scatter of points has no regular spec: the eckit grid is built from the
    # latitudes and longitudes themselves, one per point
    assert v.xr_grid.grid_spec is None
    assert v.xr_grid.grid_type == "unstructured"
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (6,)

    # the points are the pairs of the coordinates, not the mesh of them
    x, y = v.xr_grid.xys
    assert x.shape == (6,)
    assert np.max(np.abs(x - np.asarray(SCATTER_X))) == 0.0
    assert np.max(np.abs(y - np.asarray(SCATTER_Y))) == 0.0

    from pyproj import CRS, Transformer

    transformer = Transformer.from_crs(CRS.from_proj4(LAEA_PROJ4), "EPSG:4326", always_xy=True)
    ref_lon, ref_lat = transformer.transform(np.asarray(SCATTER_X), np.asarray(SCATTER_Y))

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (6,)
    assert np.max(np.abs(lat - ref_lat)) <= 1e-6
    assert np.max(np.abs((lon - ref_lon + 180.0) % 360.0 - 180.0)) <= 1e-6

    # the eckit grid which takes the data's values holds the same points
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64).reshape(-1) for a in v.ek_grid.to_latlons())
    assert np.max(np.abs(ek_lat - lat)) <= 1e-6
    assert np.max(np.abs((ek_lon - lon + 180.0) % 360.0 - 180.0)) <= 1e-6


def test_xarray_grid_projected_scatter_with_level():
    # the same scatter, with the data variable carrying a leading level: the
    # grid dimension is the one its coordinates are defined on, and the points
    # keep the order the values are held in
    ds = _projected_scatter_dataset()
    ds["a"] = xr.Variable(
        ["level", "values"],
        np.stack([ds["a"].values, -ds["a"].values]),
        attrs={"grid_mapping": "mapping"},
    )
    ds = ds.assign_coords(level=("level", np.array([700.0, 500.0])))

    v = get_variables(ds)[0]

    assert v.xr_grid.variable_dims == ("values",)
    assert list(v.geo_dims) == ["values"]
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (6,)

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (6,)
    assert np.all(np.isfinite(lat)) and np.all(np.isfinite(lon))

    # the bounding box of the points, which the grid is only usable with once
    # its latlons are
    north, south, east, west = v.xr_grid.bbox
    assert north > south
    assert east > west


def _rotated_dataset() -> xr.Dataset:
    """A small synthetic dataset on a CF rotated lat/lon grid.

    Returns
    -------
    xr.Dataset
        The dataset, with its data variable "a" over the 1D ``rlat``/``rlon``
        axes of the rotated frame, the 2D ``latitude``/``longitude``
        coordinates of the eckit-geo grid carrying that rotation (the
        convention a detected rotated gridspec is checked against), and the
        scalar CF grid mapping variable naming the rotation. Its dimension
        axes play the roles of the projection plane axes of
        :func:`_projected_dataset`: the loader classifies them as an x/y pair.
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


@pytest.mark.parametrize(
    "build,yx",
    [
        (_projected_dataset, ("y", "x")),
        (_rotated_dataset, ("rlat", "rlon")),
    ],
)
def test_xarray_grid_explicit_structured_grid_dim_order(build, yx):
    # an explicit eckit grid that iterates its points with its x axis (the
    # rotated longitude for a rotation) fastest -- a regular_xy or regular_ll
    # one -- takes the data with its (x, y) dimensions swapped to (y, x); a
    # grid detected from the dataset is instead built to match the variable's
    # own dimension order, which must therefore be kept
    x, y = yx[1], yx[0]

    ds = build()
    ds_t = ds.transpose(x, y)  # the data and the 2D coordinates in (x, y) order

    ordered = get_variables(ds)[0]
    assert ordered.ek_grid.type in ("regular_xy", "regular_ll")

    # the auto-detected grids keep the variable's own dimension order
    assert list(ordered.geo_dims) == [y, x]
    assert list(get_variables(ds_t)[0].geo_dims) == [x, y]

    # the same grid, given explicitly: both orders come out as (y, x)
    ek_grid = ordered.ek_grid
    assert list(get_variables(ds, user_ek_grid=ek_grid)[0].geo_dims) == [y, x]
    variable = get_variables(ds_t, user_ek_grid=ek_grid)[0]
    assert list(variable.geo_dims) == [y, x]
    assert variable.ek_grid is ek_grid


def test_xarray_grid_rotated_grid_detection():
    # the grid the rotated dataset's coordinates are detected as: a regular_ll
    # spec in the rotated frame, which the auto path does not reorder data
    # onto -- it builds the spec to match the variable's own dimension order,
    # with the order that scans it with the latitude fastest when the variable
    # transposed onto (rlon, rlat) does (see _ROTATED_ORDER_SWAP)
    ds = _rotated_dataset()

    ordered = get_variables(ds)[0]
    assert ordered.ek_grid.type == "regular_ll"
    assert "order" not in ordered.ek_grid.spec

    transposed = get_variables(ds.transpose("rlon", "rlat"))[0]
    assert transposed.ek_grid.type == "regular_ll"
    assert transposed.ek_grid.spec["order"] == "j-i+"
    assert list(transposed.geo_dims) == ["rlon", "rlat"]


def test_xarray_grid_explicit_outer_order_grid_keeps_dim_order():
    # the transposed rotated dataset's auto-detected grid scans its longitude
    # (the x of the rotation) axis outermost, so its points already follow the
    # data's (rlon, rlat) order: giving that grid explicitly must not swap the
    # dimensions, or the values are matched to the wrong points
    from eckit.geo import Grid

    ds_t = _rotated_dataset().transpose("rlon", "rlat")

    auto = get_variables(ds_t)[0]
    assert auto.ek_grid.type == "regular_ll"
    assert auto.ek_grid.order == "j-i+"

    variable = get_variables(ds_t, user_ek_grid=auto.ek_grid)[0]
    assert list(variable.geo_dims) == ["rlon", "rlat"]

    # the grid's points are the coordinates flattened in that same order, the
    # rotated latitude varying fastest
    lat, lon = auto.xr_grid.latlons
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64) for a in auto.ek_grid.to_latlons())
    assert np.allclose(ek_lat, lat)
    assert np.allclose(ek_lon, lon)

    # the same grid scanned with its longitude fastest instead (an i-order,
    # ascending rotated latitudes) does swap the dimensions to (rlat, rlon)
    lon_fastest = Grid(dict(auto.ek_grid.spec, order="i+j+"))
    assert list(get_variables(ds_t, user_ek_grid=lon_fastest)[0].geo_dims) == ["rlat", "rlon"]


def test_xarray_grid_explicit_unstructured_keeps_dim_order():
    # an explicit UNSTRUCTURED grid holds its points in the order the data was
    # laid out onto it: taking the (x, y)-ordered projected dataset's own
    # detected unstructured grid and giving it explicitly must not swap the
    # dimensions, or the values are matched to the wrong points
    ds_t = _projected_dataset().transpose("x", "y")

    auto = get_variables(ds_t)[0]
    assert auto.ek_grid.type == "unstructured_ll"

    variable = get_variables(ds_t, user_ek_grid=auto.ek_grid)[0]
    assert list(variable.geo_dims) == ["x", "y"]

    # the grid's points are the coordinates flattened in that same (x, y) order
    lat, lon = auto.xr_grid.latlons
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64) for a in auto.ek_grid.to_latlons())
    assert np.allclose(ek_lat, lat)
    assert np.allclose(ek_lon, lon)


def test_xarray_grid_explicit_grid_from_accessor(monkeypatch):
    # the earthkit-data accessor's grid_spec is an explicit grid too: the
    # loader reorders (x, y) data onto it when it scans x fastest, and keeps
    # the order when it does not (the accessor is the loader's
    # grid_from_earthkit, replaced here so no earthkit-data is needed)
    from earthkit.geo.grids.xarray import loader

    ds = _projected_dataset()
    ds_t = ds.transpose("x", "y")
    regular = get_variables(ds)[0].ek_grid  # the detected regular_xy grid
    unstructured = get_variables(ds_t)[0].ek_grid  # the detected unstructured one
    assert (regular.type, unstructured.type) == ("regular_xy", "unstructured_ll")

    monkeypatch.setattr(loader, "grid_from_earthkit", lambda ds: regular)
    assert list(get_variables(ds_t)[0].geo_dims) == ["y", "x"]

    monkeypatch.setattr(loader, "grid_from_earthkit", lambda ds: unstructured)
    variable = get_variables(ds_t)[0]
    assert list(variable.geo_dims) == ["x", "y"]
    assert variable.ek_grid is unstructured


def _cf_area_index(name: str, area) -> np.ndarray:
    """The global CF indices, increasing, of the points of an eckit area grid.

    Parameters
    ----------
    name : str
        The name of the global grid, e.g. ``"O32"``.
    area : Any
        The ``[north, west, south, east]`` area of it.

    Returns
    -------
    np.ndarray
        The index of every point of the grid ``{"grid": name, "area": area}``,
        sorted into the increasing order a CF file must hold them in.
    """
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    encoding = cf_reduced_gaussian_encoding(Grid({"grid": name, "area": area}))
    assert encoding is not None, (name, area)

    indices = np.asarray(encoding.indices, dtype=np.int64)
    if encoding.order is None:
        return indices
    return indices[np.asarray(encoding.order)]


def _global_grid_points(name: str) -> tuple:
    """The flat ``(lat, lon)`` points of the global grid ``name``, in its index order.

    A reduced Gaussian grid built without an area iterates exactly the points
    the CF global index numbers, west to east and north to south, so its point
    ``i`` is the one the index ``i`` names.
    """
    from eckit.geo import Grid

    lat, lon = Grid({"grid": name}).to_latlons()
    return np.asarray(lat, dtype=np.float64).reshape(-1), np.asarray(lon, dtype=np.float64).reshape(-1)


def _selective_index(name: str = "O32") -> np.ndarray:
    """The increasing global indices of a subset of ``name`` no area describes.

    Two patches far apart: the latitude lines between them hold no point of the
    file, and no area which keeps the lines of the one keeps the lines of the
    other without the rows in between, so the subset is a selective one.
    """
    return np.union1d(_cf_area_index(name, [80, -170, 75, -130]), _cf_area_index(name, [-75, 10, -80, 50]))


def _reduced_gaussian_dataset(n: int, **changes) -> xr.Dataset:
    """A small synthetic dataset on a CF reduced Gaussian grid, ``O<n>`` or ``N<n>``.

    Mirrors the structure of a real ``t2m_O1280.nc`` file: the data variable has
    no per-point latitudes and longitudes, only a ``reduced_gaussian_index``
    dimension numbering its points, and the grid those points sit on is
    described by the scalar ``reduced_gaussian`` grid mapping variable and the
    table of points per latitude it names.

    Parameters
    ----------
    n : int
        The Gaussian resolution of the grid: the dataset describes ``O<n>``, or
        ``N<n>`` when ``classical``.
    **changes
        How to depart from the grid that is detected:

        - ``classical``: describe a classical reduced Gaussian grid, with the
          subtype "normal" and the table of points per latitude eckit builds the
          ``N<n>`` grid from, instead of an octahedral one.
        - ``ascending_lat``: run the table of latitudes south to north instead
          of the usual north to south.
        - ``accumulated``: hold the running total of the points per latitude,
          named by ``accumulated_points_per_latitude``, instead of the counts.
        - ``index_offset``: number the points from that value instead of 0, as a
          subset of the global grid does.
        - ``index_values``: number the points with these global indices instead
          of the canonical numbering; the data variable then holds one value per
          index, whichever part of the grid the file covers.
        - ``index_reversed``: number them from the last point to the first.
        - ``with_points``: hold the latitude and longitude of every point as 1D
          auxiliary coordinates of the data variable along the point dimension,
          as CF allows next to the mapping and the index.
        - ``points_offset``: add that value to every auxiliary longitude (only
          meaningful with ``with_points``).
        - ``extra_index``: the name of an unrelated data variable which holds the
          canonical numbering this variable's points do not, named as a reduced
          Gaussian index; it comes first among the dataset's variables.
        - ``second_index``: the name of a second coordinate of the data variable on
          its point dimension, named as a reduced Gaussian index.
        - ``second_variable``: the name of a copy of the data variable, adding a
          second variable on the same grid to the dataset.
        - ``pl_offset``: add that value to every line's point count.
        - ``pl_swap``: swap the first two lines of the table, and with them the
          last two: a table which still mirrors itself and adds up to the points
          of the grid, but is the table of no grid eckit can name.
        - ``grid_subtype``: the subtype the mapping declares, e.g. "normal".
        - ``longitude_of_first_point``: the meridian the first point of every
          line sits on.
        - ``no_mapping``: drop the grid mapping variable.
        - ``no_attr``: drop the data variable's ``grid_mapping`` attribute.

    Returns
    -------
    xr.Dataset
        The dataset.
    """
    from eckit.geo import Grid

    if changes.get("classical"):
        # a classical grid follows no rule a table could be computed from: its
        # points per latitude, and its latitudes, are read off the grid eckit
        # builds under the N<n> name the CF subtype "normal" stands for
        lat_flat = np.asarray(Grid({"grid": f"N{n}"}).to_latlons()[0], dtype=np.float64)
        starts = np.concatenate(([0], np.flatnonzero(np.diff(lat_flat) != 0) + 1))
        counts = np.diff(np.concatenate((starts, [lat_flat.size]))).astype(np.int32)
        latitudes = lat_flat[starts]
    else:
        counts = np.concatenate((20 + 4 * np.arange(n), (20 + 4 * np.arange(n))[::-1])).astype(np.int32)

        # the Gaussian latitudes, one per line, north to south, as eckit's own grid gives them
        lat_flat = np.asarray(Grid({"grid": f"O{n}"}).to_latlons()[0], dtype=np.float64)
        lines = [lat_flat[0]]
        for value in lat_flat[1:]:
            if value != lines[-1]:
                lines.append(value)
        latitudes = np.asarray(lines)

    points = int(counts.sum())

    if changes.get("pl_swap"):
        counts = counts[[1, 0, *range(2, counts.size - 2), counts.size - 1, counts.size - 2]]

    if changes.get("pl_offset"):
        counts = counts + changes["pl_offset"]

    if changes.get("accumulated"):
        pl_name, pl_values = "apl", np.cumsum(counts)
    else:
        pl_name, pl_values = "pl", counts

    mapping = {
        "grid_mapping_name": "reduced_gaussian",
        "grid_subtype": changes.get("grid_subtype", "normal" if changes.get("classical") else "octahedral"),
        "latitude_dimension": "lat",
        "points_per_latitude" if pl_name == "pl" else "accumulated_points_per_latitude": pl_name,
    }
    if "longitude_of_first_point" in changes:
        mapping["longitude_of_first_point"] = changes["longitude_of_first_point"]

    index_values = changes.get("index_values")
    if index_values is None:
        index = np.arange(points, dtype=np.int32) + changes.get("index_offset", 0)
    else:
        index = np.asarray(index_values, dtype=np.int64)
    if changes.get("index_reversed"):
        index = index[::-1]
    points = int(index.size)

    index_attrs = {"standard_name": "reduced_gaussian_index", "units": "1"}

    data_vars = {
        "air_temperature": xr.Variable(
            ("time", "height", "reduced_gaussian_index"),
            np.zeros((1, 1, points), dtype=np.float32),
            attrs={"standard_name": "air_temperature", "units": "K", "grid_mapping": "reduced_gaussian"},
        ),
        pl_name: xr.Variable(("lat",), pl_values, attrs={"long_name": "number of points per latitude", "units": "1"}),
    }
    if changes.get("extra_index"):
        # a variable of the dataset which is not one of the data variable's
        # coordinates, holding the numbering its points do not have
        data_vars[changes["extra_index"]] = xr.Variable(
            "reduced_gaussian_index", np.arange(points, dtype=np.int32), attrs=index_attrs
        )

    coords = {
        "time": ("time", np.array(["2024-12-16"], dtype="datetime64[ns]")),
        "height": ("height", np.array([2.0]), {"units": "m"}),
        "reduced_gaussian_index": ("reduced_gaussian_index", index, index_attrs),
        "lat": ("lat", latitudes, {"standard_name": "latitude", "units": "degrees_north"}),
    }
    if changes.get("second_index"):
        # a second coordinate numbering the same points, so that the one the data
        # follows cannot be told
        coords[changes["second_index"]] = ("reduced_gaussian_index", np.arange(points, dtype=np.int32), index_attrs)

    if changes.get("with_points"):
        # CF lets a reduced Gaussian file hold the latitude and longitude of
        # every point as 1D auxiliary coordinates, alongside the index and the
        # mapping which place them just as well
        grid_name = f"N{n}" if changes.get("classical") else f"O{n}"
        point_lat, point_lon = _global_grid_points(grid_name)
        coords["latitude"] = (
            "reduced_gaussian_index",
            point_lat[index],
            {"standard_name": "latitude", "units": "degrees_north"},
        )
        coords["longitude"] = (
            "reduced_gaussian_index",
            point_lon[index] + changes.get("points_offset", 0.0),
            {"standard_name": "longitude", "units": "degrees_east"},
        )

    ds = xr.Dataset(data_vars, coords=coords)
    ds["reduced_gaussian"] = xr.DataArray(0, attrs=mapping)

    if changes.get("ascending_lat"):
        ds["lat"] = ("lat", latitudes[::-1])

    if changes.get("no_mapping"):
        del ds["reduced_gaussian"]
    if changes.get("no_attr"):
        del ds["air_temperature"].attrs["grid_mapping"]

    if changes.get("second_variable"):
        # another data variable of the dataset on the same grid, indexed by the
        # same coordinate: a copy with the final attributes of the first
        ds[changes["second_variable"]] = ds["air_temperature"]

    return ds


@pytest.mark.parametrize("n", [16, 80])
def test_xarray_grid_reduced_gaussian(n):
    from earthkit.geo.grids.xarray import ReducedGaussianGrid

    ds = _reduced_gaussian_dataset(n)
    points = ds.sizes["reduced_gaussian_index"]

    variables = get_variables(ds)

    # the grid mapping variable is skipped, and the point table has no grid
    # coordinate of its own: only the data variable is left
    assert [v.name for v in variables] == ["air_temperature"]

    v = variables[0]
    assert isinstance(v.xr_grid, ReducedGaussianGrid)
    assert v.xr_grid.grid_spec == {"grid": f"O{n}"}
    assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
    assert v.geo_dims == ("reduced_gaussian_index",)

    # the eckit grid is the octahedral one the spec names, built by eckit rather
    # than gathered from coordinates the dataset does not hold
    assert v.ek_grid.type == "reduced_gg"
    assert v.ek_grid.shape == (points,)

    lat, lon = v.xr_grid.latlons
    assert lat.shape == (points,)
    assert lon.shape == (points,)

    # the first latitude line, the one closest to the north pole: points at the
    # same latitude, spaced all the way round the globe from Greenwich
    first_row = int(ds["pl"].values[0])
    assert np.allclose(lat[:first_row], lat[0])
    assert np.allclose(lon[:first_row], 360.0 / first_row * np.arange(first_row))


def test_xarray_grid_reduced_gaussian_ascending_latitudes():
    # the direction the table of latitudes runs in says nothing about the points,
    # which the index numbers: the same grid, with the same points in the same order
    descending = _reduced_gaussian_dataset(16)
    ascending = _reduced_gaussian_dataset(16, ascending_lat=True)

    v = get_variables(ascending)[0]
    assert v.xr_grid.grid_spec == get_variables(descending)[0].xr_grid.grid_spec
    assert "order" not in v.ek_grid.spec

    lat, lon = v.xr_grid.latlons
    ek_lat, ek_lon = (np.asarray(a, dtype=np.float64) for a in v.ek_grid.to_latlons())
    assert np.array_equal(lat, ek_lat)
    assert np.array_equal(lon, ek_lon)


def test_xarray_grid_reduced_gaussian_accumulated():
    # CF spells the table of points per latitude either as the counts themselves
    # or as their running total along the latitude dimension
    v = get_variables(_reduced_gaussian_dataset(16, accumulated=True))[0]
    assert v.xr_grid.grid_spec == {"grid": "O16"}


@pytest.mark.parametrize("changes", [{}, {"accumulated": True}, {"ascending_lat": True}])
def test_xarray_grid_reduced_gaussian_normal(changes):
    # a classical grid, subtype "normal", is detected from the table of points per
    # latitude eckit builds the N<n> grid from: the same grid, the same points in
    # the same order, whichever way its table is spelled
    n = 16
    reference = get_variables(_reduced_gaussian_dataset(n, classical=True))[0]
    v = get_variables(_reduced_gaussian_dataset(n, classical=True, **changes))[0]

    assert v.xr_grid.grid_spec == {"grid": f"N{n}"}
    assert "order" not in v.ek_grid.spec

    lat, lon = v.xr_grid.latlons
    ref_lat, ref_lon = reference.xr_grid.latlons
    assert np.array_equal(lat, ref_lat)
    assert np.array_equal(lon, ref_lon)


def test_xarray_grid_reduced_gaussian_normal_netcdf_ascending(tmp_path):
    # the direction the table of latitudes runs in says nothing about the points,
    # even once the file has been written and read back: both spellings describe
    # the same grid, with the same points in the same order
    pytest.importorskip("netCDF4")

    paths = []
    for name, changes in (("t2m_N16.nc", {}), ("t2m_N16_ascending.nc", {"ascending_lat": True})):
        ds = _reduced_gaussian_dataset(16, classical=True, **changes)
        path = tmp_path / name
        ds.to_netcdf(path)
        paths.append(path)

    specs, points = [], []
    for path in paths:
        with xr.open_dataset(path) as reopened:
            v = get_variables(reopened)[0]
            specs.append(v.xr_grid.grid_spec)
            points.append(v.xr_grid.latlons)

    assert specs[0] == specs[1] == {"grid": "N16"}
    assert np.array_equal(points[0][0], points[1][0])
    assert np.array_equal(points[0][1], points[1][1])


def test_xarray_grid_reduced_gaussian_normal_netcdf_not_detected(tmp_path):
    # a normal mapping whose table belongs to no classical grid is left alone once
    # on disk too: reading the file back finds no variable and raises nothing
    pytest.importorskip("netCDF4")

    ds = _reduced_gaussian_dataset(16, classical=True, pl_swap=True)
    path = tmp_path / "t2m_not_a_classical_grid.nc"
    ds.to_netcdf(path)

    with xr.open_dataset(path) as reopened:
        assert get_variables(reopened) == []


@pytest.mark.parametrize(
    "changes",
    [
        {"grid_subtype": "normal"},  # the table of an O<n> grid is no classical table
        {"classical": True, "pl_swap": True},  # a mirrored table of the right count which is no grid's
        {"classical": True, "grid_subtype": "octahedral"},  # a classical table is no octahedral one
        {"index_offset": 1},  # the points keep the numbers they have on the global grid
        {"index_reversed": True},  # the points are numbered from the last one to the first
        {"pl_offset": 4},  # a table which is not an octahedral one
        {"longitude_of_first_point": 3.0},  # a grid which does not start at Greenwich
        {"no_mapping": True},  # the mapping variable is gone
        {"no_attr": True},  # the variable does not name it
    ],
)
def test_xarray_grid_reduced_gaussian_not_detected(changes):
    # a grid the mapping does not describe as a global reduced Gaussian one,
    # octahedral or classical, is skipped: a variable with fewer than two grid
    # coordinates is no more a variable than it was before, and nothing is raised
    # for the metadata which cannot be used
    assert get_variables(_reduced_gaussian_dataset(16, **changes)) == []


def test_xarray_grid_reduced_gaussian_index_from_the_variable():
    # the numbering a grid is detected from is the coordinate of the variable
    # itself: another variable of the dataset which carries the same standard name,
    # and the canonical index the variable's own coordinate does not, cannot stand
    # in for it, whichever of the two the dataset lists first
    ds = _reduced_gaussian_dataset(16, index_reversed=True, extra_index="global_index")

    assert list(ds.variables).index("global_index") < list(ds.variables).index("reduced_gaussian_index")
    assert list(ds["air_temperature"].coords) == ["time", "height", "reduced_gaussian_index"]

    assert get_variables(ds) == []


def test_xarray_grid_reduced_gaussian_index_ambiguous():
    # two coordinates of the variable number the same points, each one the index of
    # a global O<n> grid: which of them the data follows cannot be told, so the grid
    # is not detected rather than guessed at
    ds = _reduced_gaussian_dataset(16, second_index="other_index")

    named = [
        c for c, v in ds["air_temperature"].coords.items() if v.attrs.get("standard_name") == "reduced_gaussian_index"
    ]
    assert sorted(named) == ["other_index", "reduced_gaussian_index"]

    assert get_variables(ds) == []


# the areas of a global grid a CF index can hold: ones eckit iterates in the
# increasing index order the file holds its points in, and ones crossing the
# Greenwich meridian or wrapping round the globe, which it does not
SUBAREA_IDENTITY = [("O32", [60, 10, 30, 40]), ("N32", [60, 10, 30, 40])]
SUBAREA_ORDERED = [("O32", [60, -10, 30, 40]), ("N32", [60, -10, 30, 40]), ("O32", [60, 80, 30, 438])]


def _subarea_dataset(name: str, area, **changes) -> xr.Dataset:
    """A CF dataset on ``name`` holding exactly the points of one of its areas."""
    index = _cf_area_index(name, area)
    return _reduced_gaussian_dataset(int(name[1:]), classical=name[0] == "N", index_values=index, **changes)


def _eckit_points(grid) -> tuple:
    """The flat ``(lat, lon)`` points of an eckit grid, in its own order."""
    lat, lon = grid.to_latlons()
    return np.asarray(lat, dtype=np.float64).reshape(-1), np.asarray(lon, dtype=np.float64).reshape(-1)


def _assert_points_at_index(name: str, index: np.ndarray, lat, lon):
    """Check flat point arrays against the global grid's points an increasing index names."""
    global_lat, global_lon = _global_grid_points(name)
    assert np.allclose(lat, global_lat[index], atol=1e-9)
    dlon = np.abs((np.asarray(lon) - global_lon[index] + 180.0) % 360.0 - 180.0)
    assert np.allclose(dlon, 0.0, atol=1e-6)


@pytest.mark.parametrize("name,area", SUBAREA_IDENTITY)
def test_xarray_grid_reduced_gaussian_subarea_rectangular(name, area):
    # a CF index holding the points of an area of the global grid is detected as
    # that grid with that area, and eckit builds the file's points from the spec
    from earthkit.geo.grids.xarray import ReducedGaussianGrid

    index = _cf_area_index(name, area)
    ds = _reduced_gaussian_dataset(int(name[1:]), classical=name[0] == "N", index_values=index)

    # a proper subset of the global grid
    assert index.size < int(ds["pl"].values.sum(dtype=np.int64))

    (v,) = get_variables(ds)
    assert v.name == "air_temperature"

    assert isinstance(v.xr_grid, ReducedGaussianGrid)
    assert v.xr_grid.grid_type == "reduced_gg"
    assert set(v.xr_grid.grid_spec) == {"grid", "area"}
    assert v.xr_grid.grid_spec["grid"] == name
    assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
    assert v.geo_dims == ("reduced_gaussian_index",)

    # this area eckit iterates the way the file holds it: no permutation
    assert v.xr_grid.point_order is None

    assert v.ek_grid.type == "reduced_gg"
    assert v.ek_grid.shape == (index.size,)

    # the grid's latlons are the eckit grid's own, the file's points being its
    # first ones, and the global grid's points the index names
    lat, lon = v.xr_grid.latlons
    ek_lat, ek_lon = _eckit_points(v.ek_grid)
    assert np.array_equal(lat, ek_lat)
    assert np.array_equal(lon, ek_lon)

    _assert_points_at_index(name, index, lat, lon)


@pytest.mark.parametrize("name,area", SUBAREA_ORDERED)
def test_xarray_grid_reduced_gaussian_subarea_ordered(name, area):
    # an area crossing the Greenwich meridian, or wrapping round the globe: the
    # eckit grid of the spec runs each of its lines from the west edge of the
    # area, another order than the file's, and the grid carries the permutation
    # between the two
    from earthkit.geo.grids.xarray import ReducedGaussianGrid

    index = _cf_area_index(name, area)
    ds = _reduced_gaussian_dataset(int(name[1:]), classical=name[0] == "N", index_values=index)

    (v,) = get_variables(ds)

    assert isinstance(v.xr_grid, ReducedGaussianGrid)
    assert v.xr_grid.grid_spec["grid"] == name
    assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
    assert v.geo_dims == ("reduced_gaussian_index",)

    assert v.ek_grid.type == "reduced_gg"
    assert v.ek_grid.shape == (index.size,)

    order = v.xr_grid.point_order
    assert order is not None
    assert not np.array_equal(order, np.arange(order.size))
    assert np.array_equal(np.sort(order), np.arange(order.size))

    # the eckit grid's point p sits at position order[p] of the file
    lat, lon = v.xr_grid.latlons
    ek_lat, ek_lon = _eckit_points(v.ek_grid)
    assert np.array_equal(ek_lat, lat[order])
    assert np.array_equal(ek_lon, lon[order])

    # and the file holds the global grid's points the index names, in its order
    _assert_points_at_index(name, index, lat, lon)


def test_xarray_grid_reduced_gaussian_subarea_two_variables():
    # two variables indexed by one coordinate describe one grid: each is
    # returned, both sit on the very same grid object, and a variable which
    # looks its points up through another coordinate of another name, i.e. of
    # other contents, is never served this one from the cache
    name, area = "O32", [60, -10, 30, 40]
    ds = _subarea_dataset(name, area, second_variable="geopotential")

    variables = get_variables(ds)
    assert [v.name for v in variables] == ["air_temperature", "geopotential"]

    air_temperature, geopotential = variables
    assert air_temperature.xr_grid is geopotential.xr_grid
    assert air_temperature.geo_dims == geopotential.geo_dims == ("reduced_gaussian_index",)
    assert air_temperature.xr_grid.grid_spec == geopotential.xr_grid.grid_spec
    assert geopotential.ek_grid.type == "reduced_gg"
    assert geopotential.ek_grid.shape == air_temperature.ek_grid.shape


def test_xarray_grid_reduced_gaussian_subarea_selective():
    # a subset no area of the grid describes: the points keep the global grid's
    # latitudes and longitudes their indices name them with, and the loader
    # builds an unstructured eckit grid out of them
    from earthkit.geo.grids.xarray.grid import SelectiveReducedGaussianGrid

    index = _selective_index()
    ds = _reduced_gaussian_dataset(32, index_values=index)

    assert index.size < int(ds["pl"].values.sum(dtype=np.int64))

    (v,) = get_variables(ds)

    assert isinstance(v.xr_grid, SelectiveReducedGaussianGrid)
    assert v.xr_grid.grid_spec is None
    assert v.xr_grid.grid_type == "unstructured"
    assert v.xr_grid.point_order is None
    assert v.xr_grid.variable_dims == ("reduced_gaussian_index",)
    assert v.geo_dims == ("reduced_gaussian_index",)

    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (index.size,)

    lat, lon = v.xr_grid.latlons
    _assert_points_at_index("O32", index, lat, lon)

    ek_lat, ek_lon = _eckit_points(v.ek_grid)
    assert np.allclose(lat, ek_lat)
    assert np.allclose(lon, ek_lon)


def _assert_same_grid(ds_without_points, ds_with_points):
    """Check that per-point coordinates added to a dataset change nothing detected."""
    # the coordinates make the variable one with two grid coordinates: read
    # without the mapping it would go through the lat/lon path now
    assert "latitude" in ds_with_points["air_temperature"].coords
    assert "longitude" in ds_with_points["air_temperature"].coords

    (plain,) = get_variables(ds_without_points)
    (points,) = get_variables(ds_with_points)

    assert type(points.xr_grid) is type(plain.xr_grid)
    assert points.xr_grid.grid_spec == plain.xr_grid.grid_spec
    assert points.geo_dims == plain.geo_dims
    assert points.ek_grid.type == plain.ek_grid.type
    assert points.ek_grid.shape == plain.ek_grid.shape

    order, plain_order = points.xr_grid.point_order, plain.xr_grid.point_order
    assert (order is None) == (plain_order is None)
    if order is not None:
        assert np.array_equal(order, plain_order)

    lat, lon = points.xr_grid.latlons
    plain_lat, plain_lon = plain.xr_grid.latlons
    assert np.array_equal(lat, plain_lat)
    assert np.array_equal(lon, plain_lon)


def test_xarray_grid_reduced_gaussian_with_points():
    # a global grid which also holds its points as auxiliary coordinates is
    # still detected from the mapping: the grid, the points and the eckit grid
    # are the ones the mapping and the index describe (decision A5)
    _assert_same_grid(_reduced_gaussian_dataset(16), _reduced_gaussian_dataset(16, with_points=True))


@pytest.mark.parametrize("name,area", SUBAREA_IDENTITY + SUBAREA_ORDERED)
def test_xarray_grid_reduced_gaussian_subarea_with_points(name, area):
    # the same for a sub-area: the auxiliary coordinates, or the absence of
    # them, make no difference to the area spec and the point order
    _assert_same_grid(
        _subarea_dataset(name, area),
        _subarea_dataset(name, area, with_points=True),
    )


def test_xarray_grid_reduced_gaussian_subarea_selective_with_points():
    _assert_same_grid(
        _reduced_gaussian_dataset(32, index_values=_selective_index()),
        _reduced_gaussian_dataset(32, index_values=_selective_index(), with_points=True),
    )


@pytest.mark.parametrize("name,area", SUBAREA_IDENTITY + SUBAREA_ORDERED)
def test_xarray_grid_reduced_gaussian_subarea_points_do_not_count(name, area):
    # the auxiliary coordinates are not what the grid is built from: shifted a
    # tenth of a degree away from the points the mapping names, the detection
    # does not notice — the mapping wins (decision A5)
    _assert_same_grid(
        _subarea_dataset(name, area),
        _subarea_dataset(name, area, with_points=True, points_offset=0.1),
    )


@pytest.mark.parametrize("change", ["bad_index", "bad_table", "bad_meridian"])
def test_xarray_grid_reduced_gaussian_subarea_invalid_with_points_falls_back(change):
    # a subset the metadata cannot resolve to no reduced Gaussian grid is not
    # detected through it: with per-point coordinates the variable goes on
    # through the lat/lon path, as the unstructured scatter of points it holds,
    # and without them it is skipped, as it was before any of this
    from earthkit.geo.grids.xarray import ReducedGaussianGrid
    from earthkit.geo.grids.xarray.grid import SelectiveReducedGaussianGrid

    name, area = "O32", [60, 10, 30, 40]
    index = _cf_area_index(name, area)

    changes = {}
    if change == "bad_index":
        index = index.copy()
        index[40], index[41] = index[41], index[40]  # not strictly increasing
    elif change == "bad_table":
        changes = {"pl_offset": 4}  # the table of no grid
    else:
        changes = {"longitude_of_first_point": 3.0}  # not the Greenwich meridian

    ds = _reduced_gaussian_dataset(32, index_values=index, with_points=True, **changes)
    (v,) = get_variables(ds)

    assert not isinstance(v.xr_grid, (ReducedGaussianGrid, SelectiveReducedGaussianGrid))
    assert v.geo_dims == ("reduced_gaussian_index",)
    assert v.ek_grid.type == "unstructured_ll"
    assert v.ek_grid.shape == (index.size,)

    assert get_variables(_reduced_gaussian_dataset(32, index_values=index, **changes)) == []


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

    # the CF rotated_latitude_longitude grid mapping and the regular rlat/rlon
    # axes are recognised, so the grid is detected as a regular lat/lon grid
    # in a rotated frame rather than treated as "unstructured_ll"
    ek_grid = v.ek_grid
    assert ek_grid.type == "regular_ll"
    assert ek_grid.shape == (19, 15)

    projection = ek_grid.spec["projection"]
    assert projection["type"] == "rotation"
    # the antipode of the grid north pole (39.25, -162.0)
    assert projection["south_pole"] == pytest.approx([18.0, -39.25])

    xr_grid = v.xr_grid
    assert xr_grid.variable_dims == variable_grid_dims

    lat, lon = xr_grid.latlons
    assert lat.shape == (points_num,)
    assert lon.shape == (points_num,)
    assert np.allclose(lat, ds["lat"].values.flatten())
    assert np.allclose(lon, ds["lon"].values.flatten())

    # the eckit grid built from the spec reproduces the file's own coordinates
    ek_lat, ek_lon = ek_grid.to_latlons()
    dlat = np.max(np.abs(np.asarray(ek_lat) - ds["lat"].values.flatten()))
    dlon = np.max(np.abs((np.asarray(ek_lon) - ds["lon"].values.flatten() + 180.0) % 360.0 - 180.0))
    assert dlat <= 1e-3
    assert dlon <= 1e-3


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
