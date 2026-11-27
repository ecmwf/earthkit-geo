# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


"""Unit tests for the CF reduced Gaussian helpers.

No downloads: every grid is built from eckit, and the CF metadata a file
would hold (the grid mapping, the global pl table, the point index) is
computed from those grids.
"""

import numpy as np
import pytest

from earthkit.geo.grids.xarray.reduced_gaussian import (
    GLOBAL,
    RECTANGULAR,
    SELECTIVE,
    cf_reduced_gaussian_encoding,
    cf_reduced_gaussian_grid_name,
    cf_reduced_gaussian_subset,
)

eckit_geo = pytest.importorskip("eckit.geo")

NAMES = ["O16", "O32", "N16", "N32"]
SUBTYPE = {"O": "octahedral", "N": "normal"}


def _mapping(subtype, **kwargs):
    return {"grid_mapping_name": "reduced_gaussian", "grid_subtype": subtype, **kwargs}


def _octahedral_pl(n):
    north = 20 + 4 * np.arange(n, dtype=np.int64)
    return np.concatenate((north, north[::-1]))


def _row_sizes(lat):
    lat = np.asarray(lat, dtype=np.float64)
    bounds = np.flatnonzero(np.diff(lat) != 0) + 1
    return np.diff(np.concatenate(([0], bounds, [lat.size])))


def _normal_pl(n):
    # the points per latitude of the classical grid, read off the eckit grid
    lat, _ = eckit_geo.Grid({"grid": f"N{n}"}).to_latlons()
    return _row_sizes(lat)


def _global_pl(name):
    return _octahedral_pl(int(name[1:])) if name[0] == "O" else _normal_pl(int(name[1:]))


def _grid(name, area=None):
    spec = {"grid": name}
    if area is not None:
        spec["area"] = list(area)
    return eckit_geo.Grid(spec)


def _sorted_indices(name, area):
    """The CF index (increasing global indices) of the points of ``Grid(name, area)``.

    Independent of the module under test: the index is recovered by matching
    every point of the area grid onto the lines and columns of the global
    grid's own latitudes and longitudes.
    """
    grid_lat, grid_lon = (np.asarray(a, dtype=np.float64) for a in _grid(name, area).to_latlons())
    glat, glon = (np.asarray(a, dtype=np.float64) for a in _grid(name).to_latlons())
    lines = np.unique(glat)[::-1]  # north -> south

    starts = np.zeros(lines.size, dtype=np.int64)
    sizes = np.empty(lines.size, dtype=np.int64)
    for i, line in enumerate(lines):
        at = np.flatnonzero(glat == line)
        starts[i] = at[0]
        sizes[i] = at.size

    indices = np.empty(grid_lat.size, dtype=np.int64)
    for p in range(grid_lat.size):
        row = int(np.argmin(np.abs(lines - grid_lat[p])))
        assert abs(lines[row] - grid_lat[p]) < 1e-9
        col = int(round(np.mod(grid_lon[p], 360.0) * sizes[row] / 360.0)) % int(sizes[row])
        assert abs(np.mod(grid_lon[p] - col * 360.0 / sizes[row] + 180.0, 360.0) - 180.0) < 1e-9
        indices[p] = starts[row] + col

    return np.sort(indices)


def _point_set(grid):
    """The points of a grid, as a set of (lat, lon) with longitudes in [0, 360)."""
    lat, lon = (np.asarray(a, dtype=np.float64) for a in grid.to_latlons())
    lon = np.mod(lon, 360.0)
    lon[lon > 360.0 - 1e-6] = 0.0  # a value one rounding below 360 is the Greenwich meridian
    return set(zip(np.round(lat, 9), np.round(lon, 9)))


# --------------------------------------------------------------------------------------
# cf_reduced_gaussian_grid_name: the global grid from the mapping and the pl table
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("accumulated", [False, True])
def test_grid_name(name, accumulated):
    pl = _global_pl(name)
    subtype = SUBTYPE[name[0]]
    table = np.cumsum(pl) if accumulated else pl

    assert cf_reduced_gaussian_grid_name(_mapping(subtype), table, accumulated=accumulated) == (name, subtype)


@pytest.mark.parametrize("name", NAMES)
def test_grid_name_independent_of_points(name):
    # only the mapping and the global table come in: the number of points of
    # the file (a sub-area holds fewer than the grid) is none of this function's business
    pl = _global_pl(name)
    subtype = SUBTYPE[name[0]]
    assert cf_reduced_gaussian_grid_name(_mapping(subtype), pl) == (name, subtype)

    # the subtype spelling is what netCDF hands over: any case, padded, or bytes
    normalised = [subtype.upper(), f"  {subtype} ", subtype.encode()]
    for spelling in normalised:
        assert cf_reduced_gaussian_grid_name(_mapping(spelling), pl) == (name, subtype)

    # bytes of any spelling decode to the same subtype
    for spelling in (
        ("OCTAHEDRAL", "octahedral", "OctaHedral") if subtype == "octahedral" else ("NORMAL", "normal", "Normal")
    ):
        attrs = _mapping(spelling.encode())
        assert cf_reduced_gaussian_grid_name(attrs, pl) == (name, subtype)


def test_grid_name_greenwich_meridian():
    pl = _octahedral_pl(16)
    for value in (None, 0, 0.0, "0.0"):
        assert cf_reduced_gaussian_grid_name(_mapping("octahedral", longitude_of_first_meridian=value), pl) == (
            "O16",
            "octahedral",
        )

    for value in (10.0, -0.5, "not a number", np.array([0.0, 10.0])):
        assert cf_reduced_gaussian_grid_name(_mapping("octahedral", longitude_of_first_meridian=value), pl) is None


@pytest.mark.parametrize(
    "attrs",
    [
        None,
        "not a mapping",
        {},
        {"grid_mapping_name": "healpix"},
        _mapping("linear"),
        _mapping(42),
        _mapping(None),
    ],
)
def test_grid_name_needs_the_cf_mapping(attrs):
    assert cf_reduced_gaussian_grid_name(attrs, _octahedral_pl(16)) is None


@pytest.mark.parametrize("n", [16, 32])
def test_grid_name_rejects_pl(n):
    octahedral = _octahedral_pl(n)
    normal = _normal_pl(n)

    cases = [
        # an octahedral table which does not follow the octahedral rule
        np.r_[octahedral[:3], octahedral[3] + 4, octahedral[4:]],
        # a table which does not mirror itself
        np.concatenate((octahedral[: n // 2], octahedral[: n // 2] + 4)),
        # an odd number of lines
        octahedral[:-1],
        # no lines, a 2D table, strings, not an array
        np.array([], dtype=np.int64),
        octahedral.reshape(2, -1),
        np.array(["20"] * octahedral.size),
        [20, [24] * (octahedral.size - 1)],
        None,
    ]
    for pl in cases:
        assert cf_reduced_gaussian_grid_name(_mapping("octahedral"), pl) is None

    # each subtype is checked against its own table
    assert cf_reduced_gaussian_grid_name(_mapping("octahedral"), normal) is None
    assert cf_reduced_gaussian_grid_name(_mapping("normal"), octahedral) is None

    # a mirrored classical-shaped table of a resolution eckit does not build
    longer = np.concatenate(([normal[0]], normal, [normal[0]]))
    assert cf_reduced_gaussian_grid_name(_mapping("normal"), longer) is None


# --------------------------------------------------------------------------------------
# cf_reduced_gaussian_subset: what the index makes of the global grid
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", NAMES)
def test_subset_global(name):
    total = int(_global_pl(name).sum())

    result = cf_reduced_gaussian_subset(name, np.arange(total))

    assert result.kind == GLOBAL
    assert result.spec == {"grid": name}
    assert result.order is None
    assert result.lats is None and result.lons is None


# areas as (area, crosses_the_meridian); the wide ones hold whole lines next to partial ones
AREAS = [
    ([60, 10, 30, 40], False),
    ([60, -10, 30, 40], True),
    ([80, -100, 20, 100], True),
    ([60, 40, 30, 398], False),  # 358 degrees wide: many whole lines, some partial ones
]


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("area, crosses", AREAS)
def test_subset_rectangular(name, area, crosses):
    grid = _grid(name, area)
    index = _sorted_indices(name, area)

    result = cf_reduced_gaussian_subset(name, index)

    assert result.kind == RECTANGULAR
    assert result.grid == name
    north, west, south, east = result.spec["area"]
    assert result.spec["grid"] == name
    assert south < north
    # the area need not carry the same numbers as the one the points came from,
    # only the same points
    assert _point_set(_grid(name, result.spec["area"])) == _point_set(grid)

    # the permutation takes the file order (increasing index) to the eckit order
    if result.order is None:
        assert np.array_equal(index, np.sort(index))  # identity means the orders agree
        reordered = index
    else:
        reordered = index[result.order]

    encoding = cf_reduced_gaussian_encoding(grid)
    if crosses:
        # the eckit order of an area crossing 0 is genuinely not the increasing one
        assert result.order is not None
        assert not bool(np.all(np.diff(encoding.indices) > 0))
    assert np.array_equal(reordered, encoding.indices)


@pytest.mark.parametrize("name", NAMES)
def test_subset_rectangular_crossing_zero_permutation_is_not_identity(name):
    index = _sorted_indices(name, [60, -10, 30, 40])
    result = cf_reduced_gaussian_subset(name, index)
    assert result.kind == RECTANGULAR
    assert result.order is not None
    assert not bool(np.array_equal(result.order, np.arange(index.size)))


@pytest.mark.parametrize("name", NAMES)
def test_subset_band(name):
    area = [60, 0, 30, 360]
    index = _sorted_indices(name, area)

    result = cf_reduced_gaussian_subset(name, index)

    assert result.kind == RECTANGULAR
    north, west, south, east = result.spec["area"]
    assert west == 0.0 and east == 360.0
    assert result.order is None
    assert _point_set(_grid(name, result.spec["area"])) == _point_set(_grid(name, area))


def test_subset_wide_area_is_rectangular():
    # the reviewer's counterexample to a single largest-gap bound: O32 [60, 80, 30, 438]
    # holds 875 points, while the largest gap of every line would bound a window
    # holding all 880 of the area's lines
    grid = _grid("O32", [60, 80, 30, 438])
    assert grid.shape == (875,)
    index = _sorted_indices("O32", [60, 80, 30, 438])
    assert index.size == 875

    result = cf_reduced_gaussian_subset("O32", index)

    assert result.kind == RECTANGULAR
    assert _point_set(_grid("O32", result.spec["area"])) == _point_set(grid)


@pytest.mark.parametrize("name", NAMES)
def test_subset_single_row(name):
    pl = _global_pl(name)
    row = pl.size // 2 - 3
    start = int(pl[:row].sum())

    index = np.arange(start + 2, start + pl[row] - 2)
    result = cf_reduced_gaussian_subset(name, index)

    assert result.kind == RECTANGULAR
    north, _, south, _ = result.spec["area"]
    assert north == pytest.approx(south)

    # the points of the rebuilt grid are exactly those the index names
    rebuilt = _grid(name, result.spec["area"])
    encoding = cf_reduced_gaussian_encoding(rebuilt)
    assert np.array_equal(np.sort(encoding.indices), index)


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("size, seed", [(50, 42), (400, 7)])
def test_subset_selective_random(name, size, seed):
    pl = _global_pl(name)
    total = int(pl.sum())
    rng = np.random.default_rng(seed)
    index = np.sort(rng.choice(total, size=size, replace=False)).astype(np.int64)

    result = cf_reduced_gaussian_subset(name, index)

    assert result.kind == SELECTIVE
    assert result.spec is None
    assert result.order is None
    assert result.lats.dtype == np.float64 and result.lons.dtype == np.float64

    # the per-point coordinates are the points of the global grid at the index
    glat, glon = (np.asarray(a, dtype=np.float64) for a in _grid(name).to_latlons())
    assert np.array_equal(result.lats, glat[index])
    assert np.array_equal(result.lons, glon[index])


@pytest.mark.parametrize("name", NAMES)
def test_subset_selective_missing_row(name):
    pl = _global_pl(name)
    starts = np.concatenate(([0], np.cumsum(pl)))
    r = pl.size // 2 - 5

    # two complete lines with the line between them absent
    index = np.concatenate((np.arange(starts[r], starts[r + 1]), np.arange(starts[r + 2], starts[r + 3])))
    assert cf_reduced_gaussian_subset(name, index).kind == SELECTIVE


@pytest.mark.parametrize("name", NAMES)
def test_subset_selective_split_run(name):
    pl = _global_pl(name)
    starts = np.concatenate(([0], np.cumsum(pl)))
    r = pl.size // 2 - 5

    # one line holding two disjoint longitude runs (a point missing in their middle)
    index = np.concatenate((
        np.arange(starts[r] + 1, starts[r] + pl[r] // 2 - 1),
        np.arange(starts[r] + pl[r] // 2 + 2, starts[r] + pl[r] - 1),
        np.arange(starts[r + 1], starts[r + 2]),
    ))
    assert cf_reduced_gaussian_subset(name, index).kind == SELECTIVE


@pytest.mark.parametrize("name", NAMES)
def test_subset_rejects_index(name):
    total = int(_global_pl(name).sum())
    good = np.arange(total)

    cases = {
        "not increasing": np.array([5, 3, 1]),
        "duplicates": np.array([1, 2, 2, 3]),
        "negative": np.array([-1, 0, 5]),
        "beyond the grid": np.array([total - 1, total]),
        "not integers": np.array([1.5, 2.5]),
        "floats with a fraction": good[:10] + 0.5,
        "empty": np.array([], dtype=np.int64),
        "two dimensional": np.arange(6).reshape(2, 3),
        "strings": np.array(["1", "2"]),
        "booleans": np.array([True, False]),
        "None": None,
    }
    for label, index in cases.items():
        assert cf_reduced_gaussian_subset(name, index) is None, label


@pytest.mark.parametrize("name", ["O", "O0", "O16x", "o16", "F16", "N17", "H4", None, 42, ""])
def test_subset_rejects_grid_name(name):
    assert cf_reduced_gaussian_subset(name, np.arange(10)) is None


# --------------------------------------------------------------------------------------
# cf_reduced_gaussian_encoding: an eckit grid read back into CF metadata
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", NAMES)
def test_encoding_global(name):
    pl = _global_pl(name)
    total = int(pl.sum())

    encoding = cf_reduced_gaussian_encoding(_grid(name))

    assert encoding is not None
    assert encoding.grid == name
    assert encoding.subtype == SUBTYPE[name[0]]
    assert np.array_equal(encoding.pl, pl)
    assert encoding.pl.dtype == np.int64 and encoding.pl.size == pl.size
    assert np.array_equal(encoding.indices, np.arange(total))
    assert encoding.order is None

    glat, _ = (np.asarray(a, dtype=np.float64) for a in _grid(name).to_latlons())
    assert encoding.latitudes.dtype == np.float64
    assert encoding.latitudes.size == pl.size
    assert bool(np.all(np.diff(encoding.latitudes) < 0))  # north to south
    assert np.array_equal(encoding.latitudes, np.unique(glat)[::-1])


# the same areas as the subset tests: every point of an area grid is a global
# point, so its global index must reproduce the grid's coordinates
@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("area, crosses", AREAS + [([60, 0, 30, 360], False)])
def test_encoding_area(name, area, crosses):
    grid = _grid(name, area)
    encoding = cf_reduced_gaussian_encoding(grid)

    assert encoding is not None
    assert encoding.grid == name

    glat, glon = (np.asarray(a, dtype=np.float64) for a in _grid(name).to_latlons())
    alat, alon = (np.asarray(a, dtype=np.float64) for a in grid.to_latlons())

    # the index numbers the global points the grid holds, in the grid's own order
    assert np.array_equal(glat[encoding.indices], alat)
    dlon = np.mod(np.round(np.abs(np.mod(alon, 360.0) - glon[encoding.indices]), 12), 360.0)
    assert float(dlon.max()) < 1e-9

    increasing = bool(np.all(np.diff(encoding.indices) > 0))
    assert increasing != crosses or area[3] > 360.0
    if increasing:
        assert encoding.order is None
    else:
        assert np.array_equal(encoding.indices[encoding.order], np.sort(encoding.indices))


@pytest.mark.parametrize("name", NAMES)
def test_encoding_area_ordering_matches_the_subset_permutation(name):
    # the two permutations the helpers offer for the same crossing area are
    # inverse to each other: CF values sorted into eckit order, eckit values
    # sorted back into the increasing index order
    area = [60, -10, 30, 40]
    grid = _grid(name, area)
    encoding = cf_reduced_gaussian_encoding(grid)
    assert encoding.order is not None

    index = encoding.indices[encoding.order]
    result = cf_reduced_gaussian_subset(name, index)
    assert result.kind == RECTANGULAR and result.order is not None

    assert np.array_equal(index[result.order], encoding.indices)
    assert np.array_equal(encoding.indices[encoding.order], index)


@pytest.mark.parametrize(
    "spec",
    [
        # a rotated reduced Gaussian grid (A10)
        {"grid": "O32", "projection": {"type": "rotation", "south_pole": [18.0, -39.25]}},
        # a reduced Gaussian grid with a custom pl table
        {"type": "reduced_gg", "pl": [20] * 20 + [24] * 20 + [24] * 20 + [20] * 20},
        # not a reduced Gaussian grid at all
        {"grid": [10.0, 10.0]},
        {"grid": "H4"},
    ],
)
def test_encoding_rejects_non_plain_grids(spec):
    assert cf_reduced_gaussian_encoding(eckit_geo.Grid(spec)) is None


def test_encoding_rejects_non_grids():
    assert cf_reduced_gaussian_encoding(None) is None
    assert cf_reduced_gaussian_encoding("O32") is None


# --------------------------------------------------------------------------------------
# round trip: grid -> CF index -> sub-area spec -> grid with the same points
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize(
    "area", [[60, 10, 30, 40], [60, -10, 30, 40], [60, 0, 30, 360], [60, 80, 30, 438], [80, -100, 20, 100]]
)
def test_round_trip_area_to_index_to_area(name, area):
    grid = _grid(name, area)

    encoding = cf_reduced_gaussian_encoding(grid)
    assert encoding is not None

    index = encoding.indices if encoding.order is None else encoding.indices[encoding.order]
    assert bool(np.all(np.diff(index) > 0))

    result = cf_reduced_gaussian_subset(name, index)

    assert result.kind == RECTANGULAR
    rebuilt = _grid(name, result.spec["area"])
    # the same points, not necessarily the same area numbers
    assert _point_set(rebuilt) == _point_set(grid)
    assert rebuilt.shape == grid.shape


# --------------------------------------------------------------------------------------
# _global_geometry: the global grid every helper of the module reads its lines off
# --------------------------------------------------------------------------------------


def test_global_geometry_is_cached_and_read_only():
    from earthkit.geo.grids.xarray.reduced_gaussian import _global_geometry

    first = _global_geometry("O16")
    assert first is not None

    # a grid of millions of points is not worth materialising twice
    assert _global_geometry("O16") is first

    # names which are no plain reduced Gaussian grid are rejected without a build,
    # whether they can serve as a cache key or not
    for bad in ("F16", "", "O", None, 42, ["O16"]):
        assert _global_geometry(bad) is None

    # the cached arrays are shared between all callers: none may write into them
    for array in (first.pl, first.lats, first.starts, first.ends, first.point_lats, first.point_lons):
        assert not array.flags.writeable


def test_global_geometry_hands_out_only_writable_arrays():
    from earthkit.geo.grids.xarray.reduced_gaussian import _global_geometry

    geometry = _global_geometry("O16")

    # the tables an encoding carries are copies of the cached ones
    encoding = cf_reduced_gaussian_encoding(_grid("O16"))
    assert encoding.pl.flags.writeable and encoding.latitudes.flags.writeable
    assert encoding.indices.flags.writeable

    # and the points a selective subset locates itself with are its own array
    total = int(geometry.pl.sum())
    rng = np.random.default_rng(42)
    index = np.sort(rng.choice(total, size=40, replace=False)).astype(np.int64)
    subset = cf_reduced_gaussian_subset("O16", index)
    assert subset.kind == SELECTIVE

    subset.lats[0] += 1.0
    subset.lons[0] += 1.0

    # writing there left the cached geometry as it was
    assert geometry.point_lats.flags.writeable is False
    assert np.array_equal(_global_geometry("O16").point_lats, geometry.point_lats)
