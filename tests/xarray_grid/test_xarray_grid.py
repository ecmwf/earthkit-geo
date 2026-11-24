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
