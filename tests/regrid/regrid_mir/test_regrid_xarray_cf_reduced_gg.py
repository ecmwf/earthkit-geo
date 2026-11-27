# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""MIR regridding of xarray data on CF-encoded reduced Gaussian grids.

The INPUT half of the plan (task T2b): data on a global CF reduced Gaussian grid,
on a rectangular sub-area -- whose points eckit may iterate in another order than
the file holds them, which the handler corrects before the backend runs
(decision A2) -- and on a selective subset, handled as unstructured points, all
regridded through ``earthkit.geo.regrid`` and compared against the array API
(``earthkit.geo.grids.array.regrid``) on the very same grids. A variable taken
out of a CF dataset without its mapping and points-per-latitude tables must say
so rather than fail cryptically (decision A7), unless the caller gives
``in_grid`` and ``in_dims``, which regrids index-only data as it stands
(decision A8).

The builders near the top make a synthetic CF dataset from an eckit grid name
with an optional area or index subset, holding a smooth analytic field and a
``time`` dimension of size 2. The OUTPUT section at the bottom (task T3) checks
that regridding onto a plain reduced Gaussian target writes that same encoding
back: the CF ``reduced_gaussian`` mapping, its ``pl`` table and the increasing
``reduced_gaussian_index`` of the points, with the values gathered into that
order (decisions D2, D3, D4 and A3). The last section (task T4) runs the whole
end-to-end matrix -- every input form against every target, with the CF
encoding of the reduced Gaussian outputs checked alongside -- round-trips the
result of each sub-area input through NetCDF and regrids the reference file
itself. The section after it takes the twelve CF reduced Gaussian sub-area files
of the download server, cuts of the real ``t2m_O1280`` sample and of the fixed
global ``t2m_O32`` and ``t2m_N80`` files, and regrids each of them onto a regular
lat/lon grid and onto a reduced Gaussian area. The datasets mirror the reference
file ``t2m_O32.nc``; the three index
helpers are copied from ``tests/xarray_grid/test_xarray_grid.py``, where the
same structure is used to test the detection (folders do not import each
other's test modules, so they are duplicated here on purpose).
"""

import numpy as np
import pytest
from regrid_common import data_path, requires_mir

from earthkit.geo import regrid
from earthkit.geo.grids.array import regrid as aregrid
from earthkit.geo.grids.xarray.loader import variables as get_variables
from earthkit.geo.utils.testing import get_test_data

xr = pytest.importorskip("xarray")

# Dask, like xarray, is an optional dependency; test_regrid_xarray.py looks for
# it the same way.
try:
    import dask.array

    NO_DASK = False
except ImportError:
    NO_DASK = True

requires_dask = pytest.mark.skipif(NO_DASK, reason="No dask available")


@pytest.fixture
def in_process_scheduler():
    """Run the dask tasks of a test in this process, whatever the configured scheduler is.

    The handler binds the grids to the backend call and an eckit ``Grid`` cannot
    be pickled, so a chunked object can only be regridded by a scheduler that
    runs the calls in this process (as in ``test_regrid_xarray.py``).
    """
    with dask.config.set(scheduler="threads"):
        yield


# ---------------------------------------------------------------------------------------
# Synthetic CF reduced Gaussian datasets
# ---------------------------------------------------------------------------------------


def _field(lat, lon):
    """The smooth analytic field every dataset built here holds, on degrees."""
    return np.sin(np.deg2rad(2.0 * lat)) + np.cos(np.deg2rad(3.0 * lon))


def _global_grid_points(name: str) -> tuple:
    """The flat ``(lat, lon)`` points of the global grid ``name``, in its index order.

    A reduced Gaussian grid built without an area iterates exactly the points
    the CF global index numbers, west to east and north to south, so its point
    ``i`` is the one the index ``i`` names. (Copied from
    ``tests/xarray_grid/test_xarray_grid.py``.)
    """
    from eckit.geo import Grid

    lat, lon = Grid({"grid": name}).to_latlons()
    return np.asarray(lat, dtype=np.float64).reshape(-1), np.asarray(lon, dtype=np.float64).reshape(-1)


def _cf_area_index(name: str, area) -> np.ndarray:
    """The global CF indices, increasing, of the points of an eckit area grid.

    (Copied from ``tests/xarray_grid/test_xarray_grid.py``.)
    """
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    encoding = cf_reduced_gaussian_encoding(Grid({"grid": name, "area": area}))
    assert encoding is not None, (name, area)

    indices = np.asarray(encoding.indices, dtype=np.int64)
    if encoding.order is None:
        return indices
    return indices[np.asarray(encoding.order)]


def _selective_index(name: str = "O32") -> np.ndarray:
    """The increasing global indices of a subset of ``name`` no area describes.

    Two wide patches far apart: the latitude lines between them hold no point of
    the file, so no area which keeps the lines of the one keeps the lines of the
    other without the rows in between. (Adapted from
    ``tests/xarray_grid/test_xarray_grid.py``: the patches are wide enough for
    the points of a 5 by 5 degree output grid to fall well inside them.)
    """
    return np.union1d(
        _cf_area_index(name, [80, -170, 60, -130]),
        _cf_area_index(name, [-60, 10, -80, 50]),
    )


def _band_index(name: str, north: float = 50.0, south: float = 20.0) -> np.ndarray:
    """The increasing global indices of every point of ``name`` between two latitudes.

    A latitude line is either whole or absent, so this is a band of full rows:
    a rectangular subset which the detection must read as the area
    ``[north, 0, south, 360]`` (decision A1, which lets a full-longitude band
    keep W=0 and E=360).
    """
    lat, _ = _global_grid_points(name)
    return np.where((lat >= south) & (lat <= north))[0]


def _cf_reduced_gaussian_dataset(
    name: str,
    area=None,
    index=None,
    *,
    mapping_name: str = "reduced_gaussian",
    pl_name: str = "pl",
    with_points: bool = False,
) -> xr.Dataset:
    """A synthetic CF dataset on a reduced Gaussian grid, with an analytic field.

    The dataset mirrors the structure of a real ``t2m_O32.nc`` file: the data
    variable ``air_temperature`` is indexed by a ``reduced_gaussian_index``
    dimension and a ``time`` dimension of size 2, and the grid is described by
    the scalar ``reduced_gaussian`` grid mapping variable and the table of
    global points per latitude it names. The two time slices hold ``f`` and
    ``1.5 f + 1``, ``f`` being :func:`_field` on the points the index names.

    Parameters
    ----------
    name : str
        The eckit name of the global grid, e.g. ``"O32"`` or ``"N32"``.
    area : Optional[List]
        The ``[north, west, south, east]`` area whose eckit points the dataset
        holds. None for a global grid.
    index : Optional[np.ndarray]
        The increasing global indices the dataset holds; a subset no area
        describes yields a selective dataset. Mutually exclusive with ``area``.
    mapping_name : str
        The name of the grid mapping variable; the data variable refers to it
        by its ``grid_mapping`` attribute.
    pl_name : str
        The name of the points-per-latitude table the mapping refers to.
    with_points : bool
        Whether to hold the latitude and longitude of every point as 1D
        auxiliary coordinates of the data variable, as CF allows next to the
        index and the mapping.

    Returns
    -------
    xr.Dataset
        The dataset, its points in the increasing index order a CF file must
        hold them in.
    """
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    assert not (area is not None and index is not None), "area and index are mutually exclusive"

    # the mapping always describes the GLOBAL grid: its name, its subtype, its
    # table of points per latitude and its latitude lines come from the global grid
    global_encoding = cf_reduced_gaussian_encoding(Grid({"grid": name}))
    assert global_encoding is not None, name

    if index is None:
        grid = Grid({"grid": name, "area": area} if area is not None else {"grid": name})
        encoding = cf_reduced_gaussian_encoding(grid)
        assert encoding is not None, (name, area)
        indices = np.asarray(encoding.indices, dtype=np.int64)
        if encoding.order is not None:
            # the area grid holds its points in eckit order; a CF file must not
            indices = indices[np.asarray(encoding.order)]
    else:
        indices = np.sort(np.asarray(index, dtype=np.int64))

    assert np.array_equal(indices, np.unique(indices)), "the index must number distinct points"

    point_lat, point_lon = _global_grid_points(name)
    lat = point_lat[indices]
    lon = point_lon[indices]
    f = _field(lat, lon)
    values = np.stack([f, 1.5 * f + 1.0])

    index_attrs = {"standard_name": "reduced_gaussian_index", "units": "1"}

    coords = {
        "time": ("time", np.array(["2026-01-01", "2026-01-02"], dtype="datetime64[ns]")),
        "reduced_gaussian_index": ("reduced_gaussian_index", indices, index_attrs),
        "lat": ("lat", np.asarray(global_encoding.latitudes), {"standard_name": "latitude", "units": "degrees_north"}),
    }
    if with_points:
        coords["latitude"] = ("reduced_gaussian_index", lat, {"standard_name": "latitude", "units": "degrees_north"})
        coords["longitude"] = ("reduced_gaussian_index", lon, {"standard_name": "longitude", "units": "degrees_east"})

    ds = xr.Dataset(
        {
            "air_temperature": xr.Variable(
                ("time", "reduced_gaussian_index"),
                values,
                attrs={"standard_name": "air_temperature", "units": "K", "grid_mapping": mapping_name},
            ),
            pl_name: xr.Variable(
                ("lat",),
                np.asarray(global_encoding.pl),
                attrs={"long_name": "number of points per latitude", "units": "1"},
            ),
        },
        coords=coords,
    )
    ds[mapping_name] = xr.DataArray(
        0,
        attrs={
            "grid_mapping_name": "reduced_gaussian",
            "grid_subtype": global_encoding.subtype,
            "latitude_dimension": "lat",
            "points_per_latitude": pl_name,
        },
    )
    return ds


def _eckit_points(grid) -> tuple:
    """The flat ``(lat, lon)`` points of an eckit grid, in its own order."""
    lat, lon = grid.to_latlons()
    return np.asarray(lat, dtype=np.float64).reshape(-1), np.asarray(lon, dtype=np.float64).reshape(-1)


def _assert_differ(left, right, at_least: float = 0.5) -> None:
    """Assert that two regridded arrays really disagree, beyond their NaN positions.

    ``np.allclose`` is False as soon as either array holds a NaN, which would
    make a "these two results differ" assertion pass on identical arrays; the
    comparison therefore runs over the positions where both are finite, after
    checking the two share one NaN mask.
    """
    nan_left, nan_right = np.isnan(left), np.isnan(right)
    np.testing.assert_array_equal(nan_left, nan_right)
    both = ~nan_left
    assert int(both.sum()) > 0, "the arrays have no finite values to compare"
    assert float(np.abs(left[both] - right[both]).max()) > at_least


def _cf_gather(values, out_spec) -> np.ndarray:
    """The array-API output of ``out_spec``, in the order a CF reduced Gaussian holds it.

    A regridded plain reduced Gaussian target is written with the CF encoding
    (task T3): its values sit on the increasing global index, which for a grid
    the backend returned with an area is another order than the grid's own
    (decision A3). It is the grid the backend really produced which orders the
    result: MIR clips a global target to the area of a sub-area input and
    reports the clipped spec, so the caller hands this helper the spec the
    array API returned, not the one that was asked for.
    """
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    encoding = cf_reduced_gaussian_encoding(Grid(out_spec))
    if encoding is None or encoding.order is None:
        return values

    return values[..., np.asarray(encoding.order)]


def _cf_points(spec) -> tuple:
    """The flat ``(lat, lon)`` points of a plain reduced Gaussian ``spec``, in CF index order."""
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    grid = Grid(spec)
    lat, lon = _eckit_points(grid)
    encoding = cf_reduced_gaussian_encoding(grid)
    if encoding is not None and encoding.order is not None:
        order = np.asarray(encoding.order)
        lat, lon = lat[order], lon[order]
    return lat, lon


def _point_set(grid) -> set:
    """The ``(lat, lon)`` points of an eckit grid as a set, order and -0/0 modulo 360 aside."""
    lat, lon = _eckit_points(grid)
    return {(round(float(la), 6), round(float(lo) % 360.0, 6)) for la, lo in zip(lat, lon)}


def _ll_shape(spec) -> tuple:
    """The ``(rows, columns)`` of a regular lat/lon grid ``spec``: the grid step over its area.

    A full longitude circle holds ``360 / step`` columns, its eastern edge
    standing for the point west of the northern one; any other span closes
    both ends.
    """
    north, west, south, east = spec.get("area", [90.0, 0.0, -90.0, 360.0])
    dy, dx = spec["grid"]
    rows = int(round((north - south) / dy)) + 1
    columns = int(round((east - west) / dx)) if abs(east - west - 360.0) < 1e-9 else int(round((east - west) / dx)) + 1
    return rows, columns


# the interior of the two patches of :func:`_selective_index`, comfortably away
# from their edges so that linear interpolation between the O32 points sitting
# on them comes out near the analytic field (latitudes, longitudes in [0, 360))
SELECTIVE_INSIDE = ((65.0, 195.0, 75.0, 225.0), (-75.0, 15.0, -65.0, 45.0))


def _selective_inside(lat, lon) -> np.ndarray:
    """Whether each point lies inside one of the patches of :func:`_selective_index`."""
    lat, lon = np.asarray(lat), np.asarray(lon) % 360.0
    mask = np.zeros(np.shape(lat), dtype=bool)
    for south, west, north, east in SELECTIVE_INSIDE:
        mask |= (lat >= south) & (lat <= north) & (lon >= west) & (lon <= east)
    return mask


# ---------------------------------------------------------------------------------------
# Regridding a detected CF reduced Gaussian input (decisions A2, D1)
# ---------------------------------------------------------------------------------------

# The full input matrix of task T4: global grids, the sub-areas of the plan
# (a plain one, one crossing the Greenwich meridian, one wide one whose points
# nearly fill a latitude line), a band of full rows on a classical grid given
# by its index, and a selective subset which no area describes. The kind says
# how the detection must read the file: as the global grid, as a rectangular
# sub-area (decision A1) or as unstructured points (decision D1).
CF_MATRIX_INPUTS = [
    pytest.param("O32", {}, "global", id="global-O32"),
    pytest.param("N32", {}, "global", id="global-N32"),
    pytest.param("O32", {"area": [60, 10, 30, 40]}, "rect", id="area-east"),
    pytest.param("O32", {"area": [60, -10, 30, 40]}, "rect", id="area-crossing"),
    pytest.param("N32", {"area": [60, -10, 30, 40]}, "rect", id="area-crossing-normal"),
    pytest.param("O32", {"area": [60, 80, 30, 438]}, "rect", id="area-wide"),
    pytest.param("N32", {"index": _band_index("N32")}, "rect", id="band-N32"),
    pytest.param("O32", {"index": _selective_index()}, "selective", id="selective-O32"),
]

# The full target matrix: a regular lat/lon grid, the global octahedral grid
# and two sub-areas which MIR returns with a snapped, clipped spec. Every
# combination has a non-empty output: even the area targets, which reach
# neither patch of the selective input, come back holding their (extrapolated)
# points -- the equality with the array API below pins those values too.
CF_MATRIX_TARGETS = [
    pytest.param({"grid": [5, 5]}, id="ll-5x5"),
    pytest.param({"grid": "O16"}, id="o16"),
    pytest.param({"grid": "O16", "area": [60, -10, 30, 40]}, id="o16-area"),
    pytest.param({"grid": "N16", "area": [50, 0, 20, 30]}, id="n16-area"),
]


@requires_mir
@pytest.mark.parametrize("out_grid", CF_MATRIX_TARGETS)
@pytest.mark.parametrize("name,kwargs,kind", CF_MATRIX_INPUTS)
def test_regrid_xarray_cf_reduced_gaussian_matrix(name, kwargs, kind, out_grid):
    """The end-to-end matrix: every CF input on every target equals the array API it is detected to.

    A structured input (global grid or rectangular sub-area) is regridded on
    exactly the grid and in exactly the order the array API sees it, the
    handler's reorder (decision A2) included; a selective one is regridded as
    the unstructured points the detection settles on (decision D1), with the
    analytic field recovered wherever the target reaches into its patches. A
    reduced Gaussian target is written back with the CF encoding of task T3,
    and the result covers the grid the backend really returned -- MIR clips a
    global target to the area of a sub-area input and reports the clipped
    spec, so the comparison follows that spec, not the one asked for.
    """
    ds = _cf_reduced_gaussian_dataset(name, **kwargs)

    (variable,) = get_variables(ds)
    spec = variable.xr_grid.grid_spec
    order = variable.xr_grid.point_order

    if kind == "global":
        assert spec == {"grid": name}
        assert order is None
    elif kind == "rect":
        # a rectangular subset is detected as an area of the named grid; the
        # area numbers may differ from the ones the file was made with (A1)
        assert spec is not None and spec["grid"] == name and "area" in spec
        # eckit starts every row of an area at its western edge, the file's
        # increasing global index at longitude 0: the orders differ exactly
        # when the area wraps past that meridian, its east beyond 360 degrees
        assert (order is None) == (spec["area"][3] <= 360.0 + 1e-9)
    else:
        assert spec is None  # the unstructured treatment (decision D1)
        assert order is None

    # the file holds the values in the order of its increasing index
    lat, lon = variable.xr_grid.latlons
    np.testing.assert_allclose(ds["air_temperature"].values[0], _field(lat, lon))

    values = ds["air_temperature"].values

    result = regrid(ds, out_grid=out_grid, interpolation="linear")
    regridded = result["air_temperature"]
    assert regridded.dims[0] == "time"

    reduced_target = isinstance(out_grid["grid"], str)
    if reduced_target:
        # the CF encoding of the plain target (decisions D2, D3): the increasing
        # global index, the tables and the mapping the reference file holds
        index = result["reduced_gaussian_index"]
        assert index.dims == ("reduced_gaussian_index",)
        assert np.all(np.diff(index.values) > 0)
        assert index.attrs == {"standard_name": "reduced_gaussian_index", "units": "1"}
        assert {"pl", "reduced_gaussian"} <= set(result.data_vars)
        assert "lat" in result.variables
        assert regridded.attrs["grid_mapping"] == "reduced_gaussian"
    else:
        # a regular target holds the field on its own coordinates; the input's
        # mapping described the input grid and is gone (decision A6)
        assert "pl" not in result.data_vars
        assert "grid_mapping" not in regridded.attrs

    # the reference: the array API on the very grid the detection settled on,
    # fed the values as the handler hands them to it -- reordered for a sub-area
    # whose eckit order differs from the file's, unchanged for the selective
    # subset which the coordinates themselves locate
    in_grid = variable.ek_grid if kind == "selective" else spec
    reorder = None if kind == "selective" else order
    reference = []
    returned = None
    for field in values:
        eckit_order = field if reorder is None else field[np.asarray(reorder)]
        regridded_field, returned = aregrid(eckit_order, in_grid=in_grid, out_grid=out_grid, interpolation="linear")
        reference.append(regridded_field)
    # the result of a reduced Gaussian target holds its values in the CF order
    reference = _cf_gather(np.stack(reference).reshape(regridded.shape), returned)

    np.testing.assert_allclose(regridded.values, reference)

    # the shape of the result is the shape of the grid the backend returned
    if reduced_target:
        assert regridded.shape == (values.shape[0], len(_cf_indices_of(returned)))
    else:
        rows, columns = _ll_shape(returned)
        assert (regridded.sizes["latitude"], regridded.sizes["longitude"]) == (rows, columns)

    if kind == "selective":
        # where the target reaches into one of the two patches its scattered
        # points do define the field there; the area targets of this matrix
        # reach neither, and hold only extrapolated (or NaN) values, which the
        # equality with the array API above pins all the same
        if reduced_target:
            out_lat, out_lon = _cf_points(returned)
        else:
            out_lat, out_lon = np.meshgrid(result["latitude"].values, result["longitude"].values, indexing="ij")
        inside = _selective_inside(out_lat, out_lon)
        if "area" not in out_grid:
            # a global target certainly reaches into the patches; an area one
            # of this matrix lies away from both and defines no field value
            assert int(inside.sum()) > 0
        if int(inside.sum()) > 0:
            first = regridded.values[0]  # the time slice which holds f itself
            assert np.isfinite(first[inside]).all()
            np.testing.assert_allclose(first[inside], _field(out_lat[inside], out_lon[inside]), atol=0.5)


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_subarea_order_matters():
    """Without the reorder a crossing sub-area would be regridded onto the wrong points."""
    ds = _cf_reduced_gaussian_dataset("O32", area=[60, -10, 30, 40])

    (variable,) = get_variables(ds)
    order = variable.xr_grid.point_order
    assert order is not None

    out_grid = {"grid": "O16"}
    result = regrid(ds, out_grid=out_grid, interpolation="linear")["air_temperature"]

    # the values as the file holds them, on the grid of the detected spec: the
    # result must NOT be that, and the reorder must move it
    wrong = []
    returned = None
    for field in ds["air_temperature"].values:
        regridded, returned = aregrid(
            field, in_grid=variable.xr_grid.grid_spec, out_grid=out_grid, interpolation="linear"
        )
        wrong.append(regridded)

    _assert_differ(result.values, _cf_gather(np.stack(wrong).reshape(result.shape), returned))

    # the points the values sit on after the reorder are the eckit grid's:
    # gather the file's coordinates the same way and they match the backend's
    file_lat, file_lon = variable.xr_grid.latlons
    ek_lat, ek_lon = _eckit_points(variable.ek_grid)
    np.testing.assert_allclose(ek_lat, file_lat[np.asarray(order)])
    np.testing.assert_allclose(ek_lon, file_lon[np.asarray(order)])


@requires_mir
@requires_dask
@pytest.mark.parametrize("chunks", [{"time": 1}, {"reduced_gaussian_index": 50}], ids=["chunk-time", "chunk-index"])
def test_regrid_xarray_cf_reduced_gaussian_dask(chunks, in_process_scheduler):
    """A chunked input on a crossing sub-area regrids lazily, to the eager values."""
    ds = _cf_reduced_gaussian_dataset("O32", area=[60, -10, 30, 40])

    out_grid = {"grid": [5, 5]}
    eager = regrid(ds, out_grid=out_grid, interpolation="linear")["air_temperature"]
    lazy = regrid(ds.chunk(chunks), out_grid=out_grid, interpolation="linear")["air_temperature"]

    assert isinstance(lazy.data, dask.array.Array)
    assert dict(lazy.sizes) == dict(eager.sizes)
    np.testing.assert_allclose(lazy.compute().values, eager.values)


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_dask_matches_array_api(in_process_scheduler):
    """The lazy result of a crossing sub-area also equals the array API on the reordered values."""
    ds = _cf_reduced_gaussian_dataset("O32", area=[60, -10, 30, 40])

    (variable,) = get_variables(ds)
    order = np.asarray(variable.xr_grid.point_order)
    out_grid = {"grid": [5, 5]}

    lazy = regrid(ds.chunk({"time": 1}), out_grid=out_grid, interpolation="linear")["air_temperature"]

    reference = []
    for field in ds["air_temperature"].values:
        regridded, _ = aregrid(
            field[order], in_grid=variable.xr_grid.grid_spec, out_grid=out_grid, interpolation="linear"
        )
        reference.append(regridded)

    np.testing.assert_allclose(lazy.compute().values, np.stack(reference).reshape(lazy.shape))


# ---------------------------------------------------------------------------------------
# The in_grid + in_dims escape hatch (decision A8)
# ---------------------------------------------------------------------------------------


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_in_grid_in_dims():
    """A 1D index-only DataArray regrids on an explicit in_grid + in_dims, as the dataset does."""
    ds = _cf_reduced_gaussian_dataset("O32")
    da = ds["air_temperature"].isel(time=0)

    out_grid = {"grid": [5, 5]}
    result = regrid(
        da, in_grid={"grid": "O32"}, in_dims=["reduced_gaussian_index"], out_grid=out_grid, interpolation="linear"
    )

    reference = regrid(ds, out_grid=out_grid, interpolation="linear")["air_temperature"].isel(time=0)

    assert result.name == "air_temperature"
    assert dict(result.sizes) == dict(reference.sizes)
    np.testing.assert_allclose(result.values, reference.values)


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_in_grid_no_reordering():
    """An explicit in_grid on a crossing sub-area takes the values as they stand (decision A8)."""
    ds = _cf_reduced_gaussian_dataset("O32", area=[60, -10, 30, 40])

    (variable,) = get_variables(ds)
    spec = variable.xr_grid.grid_spec
    order = variable.xr_grid.point_order
    assert order is not None

    out_grid = {"grid": "O16"}
    result = regrid(ds, in_grid=spec, in_dims=["reduced_gaussian_index"], out_grid=out_grid, interpolation="linear")[
        "air_temperature"
    ]

    # the user stated the grid: the values are regridded in the order the file
    # holds them, with no permutation, which the array API on the same data shows
    reference = []
    returned = None
    for field in ds["air_temperature"].values:
        regridded, returned = aregrid(field, in_grid=spec, out_grid=out_grid, interpolation="linear")
        reference.append(regridded)

    assert result.dims == ("time", "reduced_gaussian_index")
    np.testing.assert_allclose(result.values, _cf_gather(np.stack(reference), returned))

    # and that is not what the detected grid gives: the detected path reorders
    detected = regrid(ds, out_grid=out_grid, interpolation="linear")["air_temperature"]
    _assert_differ(result.values, detected.values)


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_in_grid_errors():
    """An in_grid + in_dims which cannot describe the data is refused, clearly."""
    ds = _cf_reduced_gaussian_dataset("O32")
    da = ds["air_temperature"].isel(time=0)

    # a grid with another number of points than the named dimensions hold
    with pytest.raises(ValueError, match="do not describe the same grid"):
        regrid(da, in_grid={"grid": "O16"}, in_dims=["reduced_gaussian_index"], out_grid={"grid": [5, 5]})

    # a grid without the dimensions the values are laid out along
    with pytest.raises(ValueError, match="'in_dims'"):
        regrid(da, in_grid={"grid": "O32"}, out_grid={"grid": [5, 5]})

    # dimensions the variable does not have
    with pytest.raises(ValueError, match="'in_dims'"):
        regrid(da, in_grid={"grid": "O32"}, in_dims=["values"], out_grid={"grid": [5, 5]})


# ---------------------------------------------------------------------------------------
# The clear error for a variable out of reach of its CF tables (decision A7)
# ---------------------------------------------------------------------------------------

MAPPING_NAMES = [
    pytest.param("reduced_gaussian", "pl", id="default-names"),
    pytest.param("my_grid_mapping", "my_points_per_latitude", id="other-names"),
]


@requires_mir
@pytest.mark.parametrize("mapping_name,pl_name", MAPPING_NAMES)
def test_regrid_xarray_cf_reduced_gaussian_index_without_tables(mapping_name, pl_name):
    """A DataArray out of a CF dataset names both ways out when its tables are gone."""
    ds = _cf_reduced_gaussian_dataset("O32", mapping_name=mapping_name, pl_name=pl_name)
    da = ds["air_temperature"].isel(time=0)
    assert da.attrs["grid_mapping"] == mapping_name

    with pytest.raises(ValueError, match="reduced Gaussian") as error:
        regrid(da, out_grid={"grid": [5, 5]}, interpolation="linear")

    message = str(error.value)
    assert "dataset" in message  # the way out that says: pass the dataset
    assert "in_grid" in message and "in_dims" in message  # ... or give the grid yourself

    # the same for a Dataset slice which dropped the mapping and the table
    without_tables = ds[[da.name]]
    assert mapping_name not in without_tables
    with pytest.raises(ValueError, match="reduced Gaussian"):
        regrid(without_tables, out_grid={"grid": [5, 5]}, interpolation="linear")


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_one_unresolved_variable_in_a_dataset():
    """A Dataset which regrids one field but cannot place another indexed one says so (A7)."""
    ds = _cf_reduced_gaussian_dataset("O32")

    # a second field on the same point index, naming a mapping the dataset does not hold
    ds["q"] = xr.Variable(
        ("time", "reduced_gaussian_index"),
        ds["air_temperature"].values + 1.0,
        attrs={"grid_mapping": "another_mapping"},
    )

    # the resolvable air_temperature does not excuse the unresolvable q: it is named
    with pytest.raises(ValueError, match="'q'") as error:
        regrid(ds, out_grid={"grid": [5, 5]}, interpolation="linear")
    message = str(error.value)
    assert "air_temperature" not in message  # only the variable in the way is named
    assert "reduced Gaussian" in message
    assert "dataset" in message and "in_grid" in message and "in_dims" in message

    # with the mapping of both fields in the dataset, the same Dataset regrids;
    # the mapping and the points-per-latitude variables stay skipped, as always
    ds["q"].attrs["grid_mapping"] = "reduced_gaussian"
    r = regrid(ds, out_grid={"grid": [5, 5]}, interpolation="linear")
    assert sorted(r.data_vars) == ["air_temperature", "q"]

    # given in_grid + in_dims, an indexed field the mapping cannot place is the
    # user's, taken as it stands (decision A8): no error names it
    ds["q"].attrs["grid_mapping"] = "another_mapping"
    with_grid = regrid(
        ds,
        in_grid={"grid": "O32"},
        in_dims=["reduced_gaussian_index"],
        out_grid={"grid": [5, 5]},
        interpolation="linear",
    )
    assert sorted(with_grid.data_vars) == ["air_temperature", "q"]

    # q numbers the very points air_temperature does, shifted by one, so the
    # linear interpolation of both is the same shift on the output
    np.testing.assert_allclose(
        with_grid["q"].values, with_grid["air_temperature"].values + 1.0, equal_nan=True, atol=1e-9
    )


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_auxiliary_on_the_index_is_not_unresolved():
    """A variable the dataset declares auxiliary is skipped as data, tables or none (A7).

    Cell bounds and auxiliary coordinates sit on the point index like the field
    beside them, but another variable names them: the loader does not regrid
    them ON PURPOSE, and the error for a field whose tables cannot place it
    must not report them. A valid CF file whose field points at its
    auxiliaries must regrid like the file without them. (A data variable named
    like its own single dimension becomes a coordinate in xarray, so the
    auxiliary coordinate here carries another name.)
    """
    ds = _cf_reduced_gaussian_dataset("O16")
    index = ds["reduced_gaussian_index"].values
    lon = _global_grid_points("O16")[1][index]

    # a per-point variable with cell bounds, and an auxiliary coordinate, all
    # indexed by the point dimension and all named by an attribute of the field
    # or of one of those (as CF names the bounds of its longitude variable)
    ds["air_pressure"] = xr.Variable(
        ("reduced_gaussian_index",),
        lon,
        attrs={"standard_name": "air_pressure", "units": "Pa", "bounds": "air_pressure_bounds"},
    )
    ds["air_pressure_bounds"] = xr.Variable(
        ("reduced_gaussian_index", "bounds"), np.stack([lon - 5.0, lon + 5.0], axis=1)
    )
    ds["qc"] = xr.Variable(("reduced_gaussian_index",), np.arange(lon.size, dtype=np.float64))
    ds["air_temperature"].attrs["coordinates"] = "air_pressure qc"

    # the loader skips exactly the auxiliaries it was told to, and nothing else
    assert sorted(v.name for v in get_variables(ds)) == ["air_temperature"]

    # regrid does not mistake the auxiliaries for unresolved fields
    r = regrid(ds, out_grid={"grid": [5, 5]}, interpolation="linear")
    assert sorted(r.data_vars) == ["air_temperature"]


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_unresolved_named_beside_auxiliaries():
    """Beside its auxiliaries the clear error still names the real field it cannot place (A7)."""
    ds = _cf_reduced_gaussian_dataset("O16")
    index = ds["reduced_gaussian_index"].values
    lon = _global_grid_points("O16")[1][index]

    # the same auxiliaries as the test above
    ds["air_pressure"] = xr.Variable(
        ("reduced_gaussian_index",),
        lon,
        attrs={"standard_name": "air_pressure", "units": "Pa", "bounds": "air_pressure_bounds"},
    )
    ds["air_pressure_bounds"] = xr.Variable(
        ("reduced_gaussian_index", "bounds"), np.stack([lon - 5.0, lon + 5.0], axis=1)
    )
    ds["qc"] = xr.Variable(("reduced_gaussian_index",), np.arange(lon.size, dtype=np.float64))
    ds["air_temperature"].attrs["coordinates"] = "air_pressure qc"

    # and a genuine field on the point index whose mapping the dataset does not hold
    ds["other"] = xr.Variable(
        ("time", "reduced_gaussian_index"),
        ds["air_temperature"].values + 1.0,
        attrs={"grid_mapping": "another_mapping"},
    )

    with pytest.raises(ValueError, match="'other'") as error:
        regrid(ds, out_grid={"grid": [5, 5]}, interpolation="linear")
    message = str(error.value)
    assert "air_temperature" not in message  # the field whose mapping is there is fine
    # none of the auxiliaries is in the way: they are named by an attribute, and skipped
    assert "air_pressure_bounds" not in message and "qc" not in message


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_index_without_in_grid():
    """A DataArray with no grid information at all gives a clear error, not a KeyError."""
    da = xr.DataArray(np.arange(12.0), dims=("values",), name="a")

    with pytest.raises(ValueError, match="Could not determine the grid"):
        regrid(da, out_grid={"grid": [5, 5]}, interpolation="linear")


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_with_points_without_in_grid():
    """A DataArray with an index coordinate and per-point latitudes/longitudes regrids anyway."""
    ds = _cf_reduced_gaussian_dataset("O32")
    da = ds["air_temperature"].isel(time=0)
    index = da["reduced_gaussian_index"].values.astype(int)
    point_lat, point_lon = _global_grid_points("O32")

    da = da.assign_coords(
        latitude=("reduced_gaussian_index", point_lat[index], {"standard_name": "latitude", "units": "degrees_north"}),
        longitude=("reduced_gaussian_index", point_lon[index], {"standard_name": "longitude", "units": "degrees_east"}),
    )
    del da.attrs["grid_mapping"]  # the partial CF form of decision D4 names no mapping

    out_grid = {"grid": [5, 5]}
    result = regrid(da, out_grid=out_grid, interpolation="linear")

    # the coordinates locate the points wherever the mapping tables are: the
    # grid the detection settles on must sit on the coordinates the variable
    # holds, in their order, with no reordering
    (variable,) = get_variables(da.to_dataset())
    assert variable.point_order is None
    ek_lat, ek_lon = _eckit_points(variable.ek_grid)
    np.testing.assert_allclose(ek_lat, da["latitude"].values, atol=1e-6)
    dlon = np.abs((ek_lon - da["longitude"].values + 180.0) % 360.0 - 180.0)
    np.testing.assert_allclose(dlon, 0.0, atol=1e-6)

    reference, _ = aregrid(da.values, in_grid=variable.ek_grid, out_grid=out_grid, interpolation="linear")

    assert result.dims == ("latitude", "longitude")
    np.testing.assert_allclose(result.values, reference)


@requires_mir
def test_regrid_xarray_cf_reduced_gaussian_subarea_with_points_without_in_grid():
    """The same on a crossing sub-area: the per-point coordinates keep their own order."""
    ds = _cf_reduced_gaussian_dataset("O32", area=[60, -10, 30, 40], with_points=True)
    da = ds["air_temperature"].isel(time=0)
    del da.attrs["grid_mapping"]

    (variable,) = get_variables(da.to_dataset())
    assert variable.ek_grid.type == "unstructured_ll"  # no mapping: the scatter locates it
    assert variable.point_order is None

    out_grid = {"grid": "O16"}
    result = regrid(da, out_grid=out_grid, interpolation="linear")

    from eckit.geo import Grid

    reference, _ = aregrid(
        da.values,
        in_grid=Grid({"latitudes": da["latitude"].values, "longitudes": da["longitude"].values}),
        out_grid=out_grid,
        interpolation="linear",
    )

    assert result.dims == ("reduced_gaussian_index",)
    np.testing.assert_allclose(result.values, reference)


# ---------------------------------------------------------------------------------------
# Regridding ONTO a plain reduced Gaussian grid writes a CF encoding
# (decisions D2, D3, D4, A3, A6, A9, A10)
# ---------------------------------------------------------------------------------------


# The plain reduced Gaussian targets: global octahedral and classical grids, and a
# crossing sub-area of one. The tuple is (grid name, area or None, CF subtype).
CF_OUT_CASES = [
    pytest.param("O32", None, "octahedral", id="o32"),
    pytest.param("N32", None, "normal", id="n32"),
    pytest.param("O32", [60, -10, 30, 40], "octahedral", id="o32-crossing"),
]


def _out_grid(name: str, area) -> dict:
    """The target spec of one ``CF_OUT_CASES`` triple."""
    return {"grid": name} if area is None else {"grid": name, "area": area}


def _ll_input() -> xr.Dataset:
    """The two-variable global 10 by 10 degree dataset the output tests regrid.

    The data variable "f" carries a ``grid_mapping`` attribute of a mapping the
    input holds, which no regridded result may keep (decision A6).
    """
    lats = np.linspace(90.0, -90.0, 19)
    lons = np.linspace(0.0, 350.0, 36)
    da = xr.DataArray(
        _field(lats[:, None], lons[None, :]),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
        name="f",
        attrs={"grid_mapping": "regular_ll"},
    )
    ds = da.to_dataset()
    ds["g"] = ds["f"] * 2
    return ds


def _cf_indices_of(grid_spec) -> np.ndarray:
    """The increasing global CF index of every point of ``grid_spec``, from the T1 encoding."""
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    encoding = cf_reduced_gaussian_encoding(Grid(grid_spec))
    assert encoding is not None, grid_spec
    indices = np.asarray(encoding.indices, dtype=np.int64)

    return indices if encoding.order is None else indices[np.asarray(encoding.order)]


def _global_encoding(name: str):
    """The CF encoding of the global grid ``name``: its latitude lines and ``pl`` table."""
    from eckit.geo import Grid

    from earthkit.geo.grids.xarray.reduced_gaussian import cf_reduced_gaussian_encoding

    encoding = cf_reduced_gaussian_encoding(Grid({"grid": name}))
    assert encoding is not None, name

    return encoding


@requires_mir
@pytest.mark.parametrize("name,area,subtype", CF_OUT_CASES)
def test_regrid_xarray_cf_output_dataset(name, area, subtype):
    """A Dataset regridded onto a plain reduced Gaussian target is written as a CF one (D2, D3, A3)."""
    ds = _ll_input()
    out_grid = _out_grid(name, area)

    r = regrid(ds, out_grid=out_grid, interpolation="linear")

    # the tables of the encoding sit beside the regridded variables, and the
    # per-point coordinates the old output wrote are gone (decision D3)
    assert list(r.data_vars) == ["f", "g", "pl", "reduced_gaussian"]
    assert "latitude" not in r.variables and "longitude" not in r.variables

    index = r["reduced_gaussian_index"]
    assert index.dims == ("reduced_gaussian_index",)
    assert index.dtype == np.int32  # as in the reference file
    assert int(index[0]) >= 0
    assert np.all(np.diff(index.values) > 0)
    assert index.attrs == {"standard_name": "reduced_gaussian_index", "units": "1"}
    # the points are the ones of the grid which was asked for, numbered by their
    # global index in that increasing order, crossing 0 degrees west or not
    np.testing.assert_array_equal(index.values, _cf_indices_of(out_grid))

    assert r["lat"].dims == ("lat",)
    np.testing.assert_allclose(r["lat"].values, _global_encoding(name).latitudes)
    assert r["lat"].attrs == {"standard_name": "latitude", "units": "degrees_north"}

    assert r["pl"].dtype == np.int32
    np.testing.assert_array_equal(r["pl"].values, _global_encoding(name).pl)
    assert r["pl"].attrs == {"long_name": "number of points per latitude", "units": "1"}

    assert r["reduced_gaussian"].ndim == 0
    assert r["reduced_gaussian"].attrs == {
        "grid_mapping_name": "reduced_gaussian",
        "grid_subtype": subtype,
        "latitude_dimension": "lat",
        "points_per_latitude": "pl",
    }

    # the stale mapping of the input is gone, the new one names the target (A6)
    assert r["f"].attrs == {"grid_mapping": "reduced_gaussian"}
    assert r["g"].attrs == {"grid_mapping": "reduced_gaussian"}

    # the values number the index: the array API on the same grids, gathered into
    # the increasing index order of the grid the backend reported (decision A3)
    reference, returned = aregrid(ds["f"].values, in_grid={"grid": [10, 10]}, out_grid=out_grid, interpolation="linear")
    reference = _cf_gather(reference, returned)
    np.testing.assert_allclose(r["f"].values, reference)
    np.testing.assert_allclose(r["g"].values, 2.0 * reference)


@requires_mir
@pytest.mark.parametrize("name,area,subtype", CF_OUT_CASES)
def test_regrid_xarray_cf_output_dataarray(name, area, subtype):
    """A DataArray onto the same target holds the index and per-point coordinates, no tables (D4, A4)."""
    ds = _ll_input()
    r = regrid(ds["f"], out_grid=_out_grid(name, area), interpolation="linear")

    assert r.dims == ("reduced_gaussian_index",)
    assert list(r.coords) == ["reduced_gaussian_index", "latitude", "longitude"]
    assert "pl" not in r.coords and "pl" not in r.attrs

    # a DataArray can name no mapping variable to point at, and carries none (D4)
    assert "grid_mapping" not in r.attrs

    np.testing.assert_array_equal(r["reduced_gaussian_index"].values, _cf_indices_of(_out_grid(name, area)))

    # the per-point coordinates are the points of the global grid which the index
    # numbers, in the order the file holds them
    glats, glons = _global_grid_points(name)
    idx = r["reduced_gaussian_index"].values
    np.testing.assert_allclose(r["latitude"].values, glats[idx])
    np.testing.assert_allclose(r["longitude"].values % 360.0, glons[idx] % 360.0, atol=1e-9)
    assert r["latitude"].attrs["standard_name"] == "latitude"
    assert r["longitude"].attrs["standard_name"] == "longitude"


@requires_mir
def test_regrid_xarray_cf_output_user_out_dims():
    """A caller who names the point dimension keeps its name; the tables hold theirs."""
    ds = _ll_input()
    r = regrid(ds, out_grid={"grid": "O32"}, out_dims=["points"], interpolation="linear")

    assert r["f"].dims == ("points",)
    assert r["points"].dims == ("points",)
    assert r["points"].attrs["standard_name"] == "reduced_gaussian_index"
    np.testing.assert_array_equal(r["points"].values, _cf_indices_of({"grid": "O32"}))
    assert list(r.data_vars) == ["f", "g", "pl", "reduced_gaussian"]


@requires_mir
@requires_dask
def test_regrid_xarray_cf_output_dask(in_process_scheduler):
    """A chunked input regridded onto a CF target stays lazy and computes to the eager result."""
    ds = _ll_input()
    out_grid = {"grid": "O32", "area": [60, -10, 30, 40]}

    eager = regrid(ds, out_grid=out_grid, interpolation="linear")
    lazy = regrid(ds.chunk({"latitude": 5}), out_grid=out_grid, interpolation="linear")

    for name in ("f", "g"):
        assert isinstance(lazy[name].data, dask.array.Array), f"{name} was computed"
    assert dict(lazy.sizes) == dict(eager.sizes)

    xr.testing.assert_allclose(lazy.compute(), eager)


@requires_mir
@pytest.mark.parametrize("name,area,subtype", CF_OUT_CASES)
def test_regrid_xarray_cf_output_netcdf_round_trip(name, area, subtype, tmp_path):
    """The CF output written to NetCDF and read back is the same grid, detected and regridable.

    The tables and the index come back unchanged, the loader reads the file as
    the target it was written for, its points are the points of that target,
    and regridding the reopened data gives what regridding the in-memory result
    gives.
    """
    ds = _ll_input()
    out_grid = _out_grid(name, area)
    r = regrid(ds, out_grid=out_grid, interpolation="linear")

    path = tmp_path / "cf_output.nc"
    r.to_netcdf(path)
    back = xr.open_dataset(path)
    try:
        np.testing.assert_array_equal(back["reduced_gaussian_index"].values, r["reduced_gaussian_index"].values)
        np.testing.assert_allclose(back["lat"].values, r["lat"].values)
        np.testing.assert_array_equal(back["pl"].values, r["pl"].values)
        assert back["reduced_gaussian"].attrs == r["reduced_gaussian"].attrs
        assert back["f"].attrs == r["f"].attrs

        # the loader reads the written file as the grid it was written for
        variables = get_variables(back)
        assert sorted(v.name for v in variables) == ["f", "g"]
        for v in variables:
            assert v.ek_grid is not None and v.ek_grid.type == "reduced_gg"
            assert v.ek_grid.spec["grid"] == name
        # the file holds its points in the increasing index order, whatever
        # order the eckit grid of its area iterates them in (the loader's own
        # permutation, decisions A1/A2, puts the values back from there)
        assert np.all(np.diff(back["reduced_gaussian_index"].values) > 0)

        # the points the file's tables describe are the points of the target
        from eckit.geo import Grid

        def _point_set(grid):
            lat, lon = _eckit_points(grid)
            return {(round(float(la), 6), round(float(lo) % 360.0, 6)) for la, lo in zip(lat, lon)}

        assert _point_set(Grid(out_grid)) == _point_set(variables[0].ek_grid)

        # and the field it carries regrids on like the one in memory
        again_back = regrid(back, out_grid={"grid": [5, 5]}, interpolation="linear")
        again_direct = regrid(r, out_grid={"grid": [5, 5]}, interpolation="linear")
        xr.testing.assert_allclose(again_back, again_direct)
    finally:
        back.close()


@requires_mir
def test_regrid_xarray_cf_output_matches_reference_file():
    """The tables the output writes are the tables of the reference NetCDF file."""
    r = regrid(_ll_input(), out_grid={"grid": "O32"}, interpolation="linear")

    ref = xr.open_dataset(data_path("legacy", "t2m_O32.nc", subfolder="xr"))
    try:
        np.testing.assert_array_equal(r["reduced_gaussian_index"].values, ref["reduced_gaussian_index"].values)
        np.testing.assert_allclose(r["lat"].values, ref["lat"].values)
        np.testing.assert_array_equal(r["pl"].values, ref["pl"].values)
        assert r["pl"].dtype == ref["pl"].dtype
        assert r["reduced_gaussian"].attrs == ref["reduced_gaussian"].attrs
    finally:
        ref.close()


@requires_mir
def test_regrid_xarray_cf_output_rotated_reduced_gg_keeps_points():
    """A reduced Gaussian target the encoding does not cover keeps the old output form (A10)."""
    ds = _ll_input()
    out_grid = {"grid": "O32", "projection": {"type": "rotation", "south_pole": [18.0, -39.25]}}

    r = regrid(ds, out_grid=out_grid, interpolation="linear")

    assert r["f"].dims == ("values",)
    assert "latitude" in r.coords and "longitude" in r.coords
    assert "pl" not in r.variables and "reduced_gaussian" not in r.variables
    # and the stale mapping of the input is gone here too (decision A6)
    assert "grid_mapping" not in r["f"].attrs
    assert "grid_mapping" not in r["g"].attrs


@requires_mir
def test_regrid_xarray_cf_output_from_cf_input():
    """A CF dataset regridded away loses its mapping, and onto another CF grid gets that grid's (A6)."""
    ds = _cf_reduced_gaussian_dataset("O32")
    assert ds["air_temperature"].attrs["grid_mapping"] == "reduced_gaussian"

    plain = regrid(ds, out_grid={"grid": [5, 5]}, interpolation="linear")
    assert "grid_mapping" not in plain["air_temperature"].attrs
    assert "pl" not in plain.variables  # the input's tables described the input grid

    n16 = regrid(ds, out_grid={"grid": "N16"}, interpolation="linear")
    assert n16["air_temperature"].attrs["grid_mapping"] == "reduced_gaussian"
    np.testing.assert_array_equal(n16["pl"].values, _global_encoding("N16").pl)
    assert n16["reduced_gaussian"].attrs["grid_subtype"] == "normal"


@requires_mir
def test_regrid_xarray_cf_output_name_clash():
    """A target whose encoding finds its names taken by the data is refused, clearly (A9)."""
    ds = _ll_input()

    # the point dimension would sit on the name of the latitude table
    with pytest.raises(ValueError, match="point dimension"):
        regrid(ds, out_grid={"grid": "O32"}, out_dims=["lat"], interpolation="linear")

    # a surviving non-geographical coordinate would be overwritten by the table
    with_pl = ds.copy()
    with_pl["f"] = with_pl["f"].assign_coords(pl=0.5)
    with pytest.raises(ValueError, match="'pl'"):
        regrid(with_pl, out_grid={"grid": "O32"}, interpolation="linear")


@requires_mir
def test_regrid_xarray_cf_output_variable_named_like_the_point_dimension():
    """A regridded variable which takes the name of the CF point dimension is refused (A9).

    Assigning a 1D variable named like the point dimension into the result makes
    it that dimension's coordinate, and the index coordinate the encoding writes
    would overwrite the regridded values with it, in silence: the clash with the
    dimension's own name is refused like a clash with a table name.
    """
    # the default dimension name of the encoding
    ds = _ll_input()
    ds["reduced_gaussian_index"] = ds["f"] * 2.0
    with pytest.raises(ValueError, match="point dimension") as error:
        regrid(ds, out_grid={"grid": "O32"}, interpolation="linear")
    assert "'reduced_gaussian_index'" in str(error.value)

    # the same under a name the caller chose for the points
    ds = _ll_input()
    ds["points"] = ds["f"] * 2.0
    with pytest.raises(ValueError, match="'points'"):
        regrid(ds, out_grid={"grid": "O32"}, out_dims=["points"], interpolation="linear")

    # a DataArray which cannot hold the index coordinate beside its own name is
    # refused too
    with pytest.raises(ValueError, match="point dimension"):
        regrid(ds["points"], out_grid={"grid": "O32"}, out_dims=["points"], interpolation="linear")

    # the clash is refused, the use of the name is not: a variable named like
    # nothing of the encoding goes through, on the dimension which does
    ds = _ll_input()
    ds["index_like"] = ds["f"] * 2.0
    r = regrid(ds, out_grid={"grid": "O32"}, interpolation="linear")
    assert sorted(r.data_vars) == ["f", "g", "index_like", "pl", "reduced_gaussian"]
    assert r["index_like"].dims == ("reduced_gaussian_index",)


# ---------------------------------------------------------------------------------------
# End-to-end (task T4): the sub-area outputs survive NetCDF, the partial-CF
# DataArray is an input again, and the reference file regrids onto a sub-area
# ---------------------------------------------------------------------------------------

# The sub-area inputs of the round trip: every kind the detection makes of an
# index subset -- a plain area, one crossing the Greenwich meridian which
# eckit reorders, one wide area which leaves only a narrow gap on each line, a
# band of full rows on a classical grid given by its index, and the selective
# subset of unstructured points.
CF_ROUND_TRIP_INPUTS = [
    pytest.param("O32", {"area": [60, 10, 30, 40]}, id="area-east"),
    pytest.param("O32", {"area": [60, -10, 30, 40]}, id="area-crossing"),
    pytest.param("O32", {"area": [60, 80, 30, 438]}, id="area-wide"),
    pytest.param("N32", {"index": _band_index("N32")}, id="band-N32"),
    pytest.param("O32", {"index": _selective_index()}, id="selective-O32"),
]


@requires_mir
@pytest.mark.parametrize("name,kwargs", CF_ROUND_TRIP_INPUTS)
def test_regrid_xarray_cf_reduced_gaussian_subarea_round_trip(name, kwargs, tmp_path):
    """A regrid of every sub-area input survives NetCDF as a CF file and regrids on.

    MIR clips the global O16 target to the area of the input and reports the
    clipped spec; written out as a CF dataset and read back, the result must
    still describe exactly those points: the loader detects a reduced Gaussian
    grid whose eckit grid holds the SAME point set as the grid the backend
    produced -- with an inferred area which may number the points differently
    than the backend's own (decision A1) -- and regridding the reopened data
    gives what regridding the in-memory result gives.
    """
    from eckit.geo import Grid

    ds = _cf_reduced_gaussian_dataset(name, **kwargs)
    (variable,) = get_variables(ds)

    out_grid = {"grid": "O16"}
    r = regrid(ds, out_grid=out_grid, interpolation="linear")
    assert "reduced_gaussian" in r.data_vars  # the CF form of task T3
    assert np.all(np.diff(r["reduced_gaussian_index"].values) > 0)

    # the grid the backend really produced: the same array-API call on the
    # values and grid the handler feeds to it
    spec, order = variable.xr_grid.grid_spec, variable.xr_grid.point_order
    field = ds["air_temperature"].values[0]
    if spec is None:
        _, returned = aregrid(field, in_grid=variable.ek_grid, out_grid=out_grid, interpolation="linear")
    else:
        eckit_order = field if order is None else field[np.asarray(order)]
        _, returned = aregrid(eckit_order, in_grid=spec, out_grid=out_grid, interpolation="linear")

    path = tmp_path / "cf_subarea.nc"
    r.to_netcdf(path)
    back = xr.open_dataset(path)
    try:
        variables = get_variables(back)
        assert sorted(v.name for v in variables) == ["air_temperature"]
        (read,) = variables

        assert read.ek_grid.type == "reduced_gg"
        assert read.xr_grid.grid_spec["grid"] == "O16"
        np.testing.assert_array_equal(back["reduced_gaussian_index"].values, r["reduced_gaussian_index"].values)

        # the area numbers of the detection may differ from the ones the
        # backend reported; the points they describe are one set
        assert _point_set(read.ek_grid) == _point_set(Grid(returned))

        # and it regrids on: the reopened file gives what the in-memory result gives
        again_back = regrid(back, out_grid={"grid": [5, 5]}, interpolation="linear")
        again_direct = regrid(r, out_grid={"grid": [5, 5]}, interpolation="linear")
        xr.testing.assert_allclose(again_back, again_direct)
    finally:
        back.close()


@requires_mir
def test_regrid_xarray_cf_partial_dataarray_regrids_again():
    """The partial-CF DataArray output of a regrid is an input again, with no in_grid (D4, A7)."""
    da = regrid(_ll_input()["f"], out_grid={"grid": "O16", "area": [60, -10, 30, 40]}, interpolation="linear")

    # the partial form: the increasing index and per-point coordinates, and no
    # mapping it could name (a DataArray holds no table for one)
    assert da.dims == ("reduced_gaussian_index",)
    assert list(da.coords) == ["reduced_gaussian_index", "latitude", "longitude"]
    assert "grid_mapping" not in da.attrs
    assert np.all(np.diff(da["reduced_gaussian_index"].values) > 0)

    # its own coordinates locate it: the detection goes through them rather
    # than raising the error an index-only variable meets (decision A7)
    (variable,) = get_variables(da.to_dataset())
    assert variable.ek_grid.type == "unstructured_ll"
    assert variable.point_order is None

    # onto a regular target, equal to the array API on the very grid
    again = regrid(da, out_grid={"grid": [5, 5]}, interpolation="linear")
    reference, _ = aregrid(da.values, in_grid=variable.ek_grid, out_grid={"grid": [5, 5]}, interpolation="linear")
    assert again.dims == ("latitude", "longitude")
    np.testing.assert_allclose(again.values, reference)

    # and onto another reduced Gaussian target, where the index comes back CF-encoded
    again_rg = regrid(da, out_grid={"grid": "N16"}, interpolation="linear")
    assert again_rg.dims == ("reduced_gaussian_index",)
    index = again_rg["reduced_gaussian_index"].values
    assert np.all(np.diff(index) > 0)
    np.testing.assert_array_equal(index, _cf_indices_of({"grid": "N16"}))
    reference, returned = aregrid(da.values, in_grid=variable.ek_grid, out_grid={"grid": "N16"}, interpolation="linear")
    np.testing.assert_allclose(again_rg.values, _cf_gather(reference, returned))


@requires_mir
@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(90)
def test_regrid_xarray_cf_reduced_gaussian_t2m_O32_to_area():
    """The reference file ``t2m_O32.nc`` regrids onto an O16 sub-area, written back as a CF file."""
    path = get_test_data("t2m_O32.nc", subfolder="xr")

    out_grid = {"grid": "O16", "area": [60, -10, 30, 40]}
    with xr.open_dataset(path) as ds:
        (variable,) = get_variables(ds)
        assert variable.xr_grid.grid_spec == {"grid": "O32"}  # the file is global

        r = regrid(ds, out_grid=out_grid, interpolation="linear")

        # the file holds its field with leading singletons: (time, height, index)
        leading = ds["air_temperature"].shape[:-1]

        # the array API on the same grid: the values the CF output must hold
        reference, returned = aregrid(
            ds["air_temperature"].values, in_grid={"grid": "O32"}, out_grid=out_grid, interpolation="linear"
        )

    values = r["air_temperature"].values
    index = r["reduced_gaussian_index"].values

    # the encoding of the target: an increasing index over the grid the backend
    # returned, the global tables and the octahedral mapping the data names
    assert np.all(np.diff(index) > 0)
    np.testing.assert_array_equal(index, _cf_indices_of(returned))
    assert list(r.data_vars) == ["air_temperature", "pl", "reduced_gaussian"]
    assert r["air_temperature"].attrs["grid_mapping"] == "reduced_gaussian"
    assert r["reduced_gaussian"].attrs["grid_subtype"] == "octahedral"
    assert "latitude" not in r.variables and "longitude" not in r.variables
    np.testing.assert_array_equal(r["pl"].values, _global_encoding("O16").pl)
    np.testing.assert_allclose(r["lat"].values, _global_encoding("O16").latitudes)

    # the input is global, so nothing is clipped away from the target area:
    # every point holds the field, equal to the array API in the CF order
    assert values.shape == leading + (index.size,)
    assert np.isfinite(values).all()
    reference = _cf_gather(np.broadcast_to(reference, values.shape).copy(), returned)
    np.testing.assert_allclose(values, reference)


# ---------------------------------------------------------------------------------------
# The sub-area files on the download server: the twelve cuts of the real t2m_O1280
# sample and of the fixed global t2m_O32 and t2m_N80 files, made by
# tools/test-data/make_cf_reduced_gaussian_files.py, each regridded onto a regular
# lat/lon grid and onto a reduced Gaussian area
# ---------------------------------------------------------------------------------------

# The file, its regular target and its reduced Gaussian target. Every target here is
# coarser than the input it takes and overlaps the area the cut holds, so that MIR
# clips it to that area and returns the clipped spec. The band of the O1280 file is
# the exception: it is a pair of latitude lines 0.21 degrees tall, and MIR crops a
# target coarser than that until no line of it is left in the band, whereupon it
# asserts on the empty interpolation ("FiniteElement: non-empty interpolation"). Its
# targets are therefore the 0.1 degree grid, which puts two lines inside the band,
# and the O320 lines of it.
CF_T2M_REGRID_INPUTS = [
    pytest.param(
        "t2m_O1280_area_east.nc", {"grid": [1, 1]}, {"grid": "O32", "area": [60, -10, 30, 40]}, id="O1280-area-east"
    ),
    pytest.param(
        "t2m_O1280_area_cross0.nc",
        {"grid": [1, 1]},
        {"grid": "O32", "area": [60, -10, 30, 40]},
        id="O1280-area-cross0",
    ),
    pytest.param(
        "t2m_O1280_band.nc", {"grid": [0.1, 0.1]}, {"grid": "O320", "area": [46, 0, 44, 360]}, id="O1280-band"
    ),
    pytest.param(
        "t2m_O1280_selective.nc",
        {"grid": [1, 1]},
        {"grid": "O32", "area": [60, -10, 30, 40]},
        id="O1280-selective",
    ),
    pytest.param(
        "t2m_O32_area_east.nc", {"grid": [1, 1]}, {"grid": "O16", "area": [65, -20, 25, 50]}, id="O32-area-east"
    ),
    pytest.param(
        "t2m_O32_area_cross0.nc",
        {"grid": [1, 1]},
        {"grid": "O16", "area": [65, -20, 25, 50]},
        id="O32-area-cross0",
    ),
    pytest.param("t2m_O32_band.nc", {"grid": [1, 1]}, {"grid": "O16", "area": [50, 0, 10, 360]}, id="O32-band"),
    pytest.param(
        "t2m_O32_selective.nc", {"grid": [1, 1]}, {"grid": "O16", "area": [65, -20, 25, 50]}, id="O32-selective"
    ),
    pytest.param(
        "t2m_N80_area_east.nc", {"grid": [1, 1]}, {"grid": "N32", "area": [65, -20, 25, 50]}, id="N80-area-east"
    ),
    pytest.param(
        "t2m_N80_area_cross0.nc",
        {"grid": [1, 1]},
        {"grid": "N32", "area": [65, -20, 25, 50]},
        id="N80-area-cross0",
    ),
    pytest.param("t2m_N80_band.nc", {"grid": [1, 1]}, {"grid": "N32", "area": [50, 0, 10, 360]}, id="N80-band"),
    pytest.param(
        "t2m_N80_selective.nc", {"grid": [1, 1]}, {"grid": "N32", "area": [65, -20, 25, 50]}, id="N80-selective"
    ),
]


@requires_mir
@pytest.mark.long_test
@pytest.mark.download
@pytest.mark.timeout(300)
@pytest.mark.parametrize("filename,ll_target,rg_target", CF_T2M_REGRID_INPUTS)
def test_regrid_xarray_cf_reduced_gaussian_t2m_subarea(filename, ll_target, rg_target):
    """Every sub-area file of the server regrids like the array API on the grid it is detected as.

    The files hold their temperatures as float32, under a ``time`` and a
    ``height`` singleton: the xarray result keeps both of them, the array API
    regrids one field and answers without them, and the two agree to within an
    ulp of float32. Whatever value either of them places is a temperature of the
    area it came out of.
    """
    path = get_test_data(filename, subfolder="xr")

    with xr.open_dataset(path) as ds:
        (variable,) = get_variables(ds)

        # what the detection settles on: an area of the grid the mapping names, with
        # the permutation of a crossing one, or unstructured points to give to eckit
        spec = variable.xr_grid.grid_spec
        order = variable.xr_grid.point_order
        in_grid = spec if spec is not None else variable.ek_grid

        values = ds["air_temperature"].values
        assert values.dtype == np.float32 and values.ndim == 3
        assert values.min() >= 200.0 and values.max() <= 330.0

        # the values in the order the handler passes them to the backend
        eckit_values = values if order is None else values[..., np.asarray(order)]

        for out_grid in (ll_target, rg_target):
            result = regrid(ds, out_grid=out_grid, interpolation="linear")
            regridded = result["air_temperature"]
            got = regridded.values

            # the reference: the array API on the very same grid and values. It
            # regrids one field and answers without the two singleton dimensions the
            # file holds before its points, which the xarray result keeps; a
            # reference for a reduced Gaussian target is then gathered into the CF
            # order the encoding writes it in
            reference, returned = aregrid(eckit_values, in_grid=in_grid, out_grid=out_grid, interpolation="linear")
            reference = _cf_gather(np.broadcast_to(reference, got.shape), returned)

            # an ulp of float32 at these temperatures is 3e-5: the two paths differ
            # by no more than the rounding of the values they write out
            np.testing.assert_allclose(got, reference, rtol=0, atol=1e-4)

            finite = np.isfinite(got)
            assert int(finite.sum()) > 0, (filename, out_grid)
            assert float(got[finite].min()) >= 200.0 and float(got[finite].max()) <= 330.0, (filename, out_grid)

            if spec is None:
                # a selective subset of scattered points covers part of the target
                # only, which MIR fills with extrapolated values; the comparison
                # above holds for the whole of it
                continue

            # a rectangular input covers the target MIR clips it to: only the outer
            # edge of the clip, which no input point surrounds, is left without a value
            assert float(finite.mean()) > 0.5, (filename, out_grid)

            if isinstance(out_grid["grid"], str):
                # the CF encoding of the plain target: the increasing global index
                # the comparison follows, the tables and the mapping the data names
                index = result["reduced_gaussian_index"].values
                assert np.all(np.diff(index) > 0)
                assert regridded.dims == ("time", "height", "reduced_gaussian_index")
                assert result["air_temperature"].attrs["grid_mapping"] == "reduced_gaussian"
                assert {"pl", "reduced_gaussian"} <= set(result.data_vars)
                np.testing.assert_array_equal(index, _cf_indices_of(returned))
