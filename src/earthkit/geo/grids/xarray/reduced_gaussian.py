# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Pure helpers for CF encoded reduced Gaussian grids, global or sub-area.

The CF conventions (Appendix F, "Reduced Gaussian") describe a reduced Gaussian
grid with a ``reduced_gaussian`` grid mapping whose parameters always refer to
the *global* grid: the table of points per latitude and the variable of
latitudes hold all 2N lines, whichever part of the Earth a file covers. Only
the ``reduced_gaussian_index`` coordinate says which global points the file
holds: 0-based serial numbers of the global grid (west to east within a
latitude line, north to south, the first point of every line at the Greenwich
meridian), non-negative, strictly increasing and smaller than the number of
points of the global grid. A file holding every index is the global grid; any
other subset is a sub-area, regional or selective (e.g. land only).

These functions translate between that CF description and the eckit grid
specs ``{"grid": "O<n>"}``/``{"grid": "N<n>"}``, with or without an ``area``.
They hold no state and touch no dataset: they take the mapping attributes, the
``pl`` table and the index values a caller has read from a file, or an
``eckit.geo.Grid`` to encode back into CF, and return plain values. Detecting
the coordinates in an xarray object is the job of
:mod:`~earthkit.geo.grids.xarray.guesser`, and writing a CF Dataset is not
done here either.

A sub-area whose points form a longitude/latitude rectangle is given the
gridspec ``{"grid": name, "area": [north, west, south, east]}`` whose points
are exactly the file's. Because the area of such a spec is only a bound, not a
list of points, eckit keeps its own point order, which runs each latitude line
from the west edge of the area: for an area crossing the Greenwich meridian
that order is not the increasing global index order the CF file is in. The
detection therefore also returns the permutation which takes the file's values
to the eckit order, and the encoding function the permutation which takes the
eckit values back to the increasing index order of a CF file.

A subset which is not a rectangle of points (a row missing between two others,
a line holding two disjoint longitude runs, ...) is a "selective" one: this
module only computes the latitude and longitude of every point of it, from the
global grid and the index, leaving it to the caller to treat them as
unstructured points.
"""

from __future__ import annotations

import functools
import re
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

from .grid_spec import _cf_reduced_gaussian_table, _classical_pl, _integer_row

# the three forms a CF index subset of a reduced Gaussian grid can take
GLOBAL = "global"
RECTANGULAR = "rectangular"
SELECTIVE = "selective"

# the eckit grid name of a plain reduced Gaussian grid: the subtype letter and
# the number of latitude lines of one hemisphere
_NAME_PATTERN = re.compile(r"^([ON])([1-9][0-9]*)$")

# how far the latitude of a point may sit from the latitude line it is on, and
# how far its longitude from the column of that line it belongs to, both in
# degrees: the area grid and the global grid compute their lines from the same
# roots, so any real difference is floating point noise, while a point of a
# genuinely different grid is off by degrees
_LAT_TOL = 1e-7
_LON_TOL = 1e-6

# the grid spec keys a plain reduced Gaussian grid may carry: anything else in
# the spec of a grid (a projection, a rotation, a custom pl, ...) is a grid
# these functions do not encode (decision A10)
_PLAIN_SPEC_KEYS = frozenset({"grid", "area"})


@dataclass
class _Geometry:
    """The global grid an eckit name describes, read off eckit itself.

    Attributes
    ----------
    name : str
        The grid name, e.g. ``"O32"``.
    subtype : str
        ``"octahedral"`` or ``"normal"``.
    pl : np.ndarray
        Number of points on every latitude line, north to south, int64 (2N,).
    lats : np.ndarray
        Latitude of every line, north to south, float64 (2N,).
    starts : np.ndarray
        Global index of the first point of every line, int64 (2N,).
    ends : np.ndarray
        Global index of the last point of every line, int64 (2N,).
    total : int
        Number of points of the global grid.
    point_lats, point_lons : np.ndarray
        Latitude and longitude of every global point, in the global index
        order, float64 (total,).
    """

    name: str
    subtype: str
    pl: np.ndarray = field(repr=False)
    lats: np.ndarray = field(repr=False)
    starts: np.ndarray = field(repr=False)
    ends: np.ndarray = field(repr=False)
    total: int
    point_lats: np.ndarray = field(repr=False)
    point_lons: np.ndarray = field(repr=False)


@dataclass
class ReducedGaussianSubset:
    """What a CF ``reduced_gaussian_index`` subset of a global grid is."""

    grid: str
    """The name of the global grid the index numbers, e.g. ``"O32"``."""

    kind: str
    """``"global"`` when the subset is the whole grid, ``"rectangular"`` when
    its points are those of an ``area`` of it, and ``"selective"`` for any
    other subset."""

    spec: Optional[dict] = None
    """The gridspec of the subset: ``{"grid": grid}`` for a global grid,
    ``{"grid": grid, "area": [north, west, south, east]}`` for a rectangular
    one, ``None`` for a selective one. The area of a rectangular spec is only
    an equivalent to the file's points, not necessarily the area the file was
    written from: eckit keeps the points of a line between its west and east
    bound, whichever area bound produced them."""

    order: Optional[np.ndarray] = field(default=None, repr=False)
    """For a rectangular subset, the permutation which takes the points from
    the order the file holds them in (increasing global index) to the order of
    the eckit grid built from ``spec``: ``values_in_file_order[order]`` are the
    values in eckit order. ``None`` when the two orders are the same, which
    they always are for a global grid and for any area not crossing the
    Greenwich meridian."""

    lats: Optional[np.ndarray] = field(default=None, repr=False)
    """For a selective subset, the latitude of every point, in the order of the
    file's index (increasing), float64. ``None`` for a global or rectangular
    subset, whose points an eckit grid provides."""

    lons: Optional[np.ndarray] = field(default=None, repr=False)
    """For a selective subset, the longitude of every point, in the order of
    the file's index (increasing), float64. ``None`` for a global or
    rectangular subset, whose points an eckit grid provides."""


@dataclass
class ReducedGaussianEncoding:
    """The CF description of an eckit reduced Gaussian grid.

    An instance carries what a writer needs to encode the grid with the CF
    ``reduced_gaussian`` grid mapping: the mapping names the global grid, the
    latitude line table and the points per latitude are always the global ones,
    and the points of the grid are named by their global index.
    """

    grid: str
    """The name of the global grid, e.g. ``"O32"``."""

    subtype: str
    """``"octahedral"`` or ``"normal"``, the value of the CF ``grid_subtype``."""

    pl: np.ndarray = field(repr=False)
    """The global points per latitude, north to south, int64 (2N,): the CF
    ``points_per_latitude`` table."""

    latitudes: np.ndarray = field(repr=False)
    """The global latitude lines, north to south, float64 (2N,): the CF
    ``latitude_dimension`` variable."""

    indices: np.ndarray = field(repr=False)
    """The global index of every point of the grid, int64 (M,), in the order the
    grid holds it (which for an area crossing the Greenwich meridian is not the
    increasing order a CF file must be in)."""

    order: Optional[np.ndarray] = field(default=None, repr=False)
    """The permutation which sorts ``indices`` into the increasing order a CF
    ``reduced_gaussian_index`` must hold: ``indices[order]`` is increasing and
    ``values[order]`` reorders a field of the grid into that order. ``None``
    when the grid is already in it."""


def _eckit_grid(spec: dict) -> Any:
    """The ``eckit.geo.Grid`` built from ``spec``, or an exception if there is none."""
    from eckit.geo import Grid as EckitGrid

    return EckitGrid(spec)


def _parse_grid_name(name: Any) -> Optional[tuple]:
    """``(subtype, rows)`` of a plain reduced Gaussian grid name, or ``None``.

    The name of the global grid an octahedral ``O<n>`` or classical ``N<n>``
    is identified by: ``subtype`` the CF ``grid_subtype`` it maps to and
    ``rows`` the number of latitude lines of one hemisphere.
    """
    if not isinstance(name, str):
        return None

    match = _NAME_PATTERN.match(name)
    if match is None:
        return None

    subtype = "octahedral" if match.group(1) == "O" else "normal"
    return subtype, int(match.group(2))


def _global_geometry(name: Any) -> Optional[_Geometry]:
    """The global grid ``name`` names, as eckit builds it, or ``None``.

    The line table, the latitudes and the points of the grid itself are read
    off the eckit grid: the sub-area specs detected from them have to match
    eckit's points exactly, and only eckit says which lines a grid holds at
    which latitude.

    Notes
    -----
    The geometry of every global grid of the same name is built once and
    returned from a small cache: a dataset's variables and a regrid's probes
    ask for the grid of their resolution over and over, and an operational
    resolution materialises millions of points each time. The cached arrays are
    read-only, so no caller can corrupt what the next one is handed; everything
    these functions give out — a subset's per-point coordinates, an encoding's
    tables — is a fresh array or a copy, writable as it was.
    """
    # only a plain grid name is ever worth building (and only a str is hashable
    # enough to key the cache with); reject anything else up front
    if _parse_grid_name(name) is None:
        return None

    return _cached_global_geometry(name)


@functools.lru_cache(maxsize=4)
def _cached_global_geometry(name: str) -> Optional[_Geometry]:
    """The one :func:`_global_geometry` of a validated grid ``name``, built on first ask."""
    subtype, rows = _parse_grid_name(name)

    try:
        grid = _eckit_grid({"grid": name})
    except Exception:
        return None

    point_lats, point_lons = grid.to_latlons()
    point_lats = np.asarray(point_lats, dtype=np.float64).reshape(-1)
    point_lons = np.asarray(point_lons, dtype=np.float64).reshape(-1)

    if point_lats.shape != point_lons.shape or point_lats.size == 0:
        return None
    if not (bool(np.all(np.isfinite(point_lats))) and bool(np.all(np.isfinite(point_lons)))):
        return None

    # eckit runs a reduced Gaussian grid line by line, north to south, the
    # points of one line sharing its latitude
    bounds = np.flatnonzero(np.diff(point_lats) != 0) + 1
    starts = np.concatenate(([0], bounds)).astype(np.int64)
    sizes = np.diff(np.concatenate((starts, [point_lats.size]))).astype(np.int64)
    lats = point_lats[starts]

    if starts.size != 2 * rows or not bool(np.all(np.diff(lats) < 0)):
        return None

    ends = starts + sizes - 1
    geometry = _Geometry(
        name=name,
        subtype=subtype,
        pl=sizes,
        lats=lats,
        starts=starts,
        ends=ends,
        total=int(sizes.sum()),
        point_lats=point_lats,
        point_lons=point_lons,
    )

    # the one geometry is every later caller's geometry: hand it out read-only
    for array in (geometry.pl, geometry.lats, geometry.starts, geometry.ends, geometry.point_lats, geometry.point_lons):
        array.flags.writeable = False

    return geometry


def _global_indices(geo: _Geometry, lats: Any, lons: Any) -> Optional[np.ndarray]:
    """The global index of every ``(lat, lon)`` point of the grid ``geo`` names.

    Each point must sit on a line of the global grid, at one of that line's
    columns; a point which does not (a rotated grid, another resolution, a
    custom table) has no index, and the whole array is rejected. Longitudes
    are taken modulo 360, where eckit may hand a point of the line west of the
    Greenwich meridian as a negative value or as a residual of one which is
    just below 360.
    """
    lats = np.asarray(lats, dtype=np.float64).reshape(-1)
    lons = np.asarray(lons, dtype=np.float64).reshape(-1)

    if lats.shape != lons.shape or lats.size == 0:
        return None
    if not (bool(np.all(np.isfinite(lats))) and bool(np.all(np.isfinite(lons)))):
        return None

    # the nearest line to each latitude; area grids hold the global lines themselves
    ascending = geo.lats[::-1]
    pos = np.searchsorted(ascending, lats)
    left = np.clip(pos - 1, 0, ascending.size - 1)
    right = np.clip(pos, 0, ascending.size - 1)
    near = np.where(ascending[right] - lats < lats - ascending[left], right, left)
    if bool(np.any(np.abs(ascending[near] - lats) > _LAT_TOL)):
        return None
    rows = ascending.size - 1 - near

    pl = geo.pl[rows]
    cols = np.rint((lons % 360.0) * pl / 360.0).astype(np.int64) % pl
    delta = (lons - cols * (360.0 / pl) + 180.0) % 360.0 - 180.0
    if bool(np.any(np.abs(delta) > _LON_TOL)):
        return None

    return geo.starts[rows] + cols


def _verify_area(geo: _Geometry, index: np.ndarray, first: int, last: int, west: float, east: float) -> Optional[tuple]:
    """The ``(area, order)`` of the candidate window, when eckit agrees it is the subset.

    The candidate is accepted only when the points of the eckit grid built
    from it are exactly the points the index names (A1): the window itself is
    an inference, eckit's points are the fact.
    """
    area = [float(geo.lats[first]), float(west), float(geo.lats[last]), float(east)]

    try:
        grid = _eckit_grid({"grid": geo.name, "area": area})
    except Exception:
        return None

    lats, lons = grid.to_latlons()
    indices = _global_indices(geo, lats, lons)
    if indices is None or indices.size != index.size or not bool(np.array_equal(np.sort(indices), index)):
        return None

    # the file holds its points at increasing index, the grid in its own order:
    # order[p] is where in the file the p-th point of the grid is
    order = np.searchsorted(index, indices)
    if bool(np.array_equal(order, np.arange(order.size))):
        order = None

    return area, order


def _rectangular(geo: _Geometry, index: np.ndarray) -> Optional[tuple]:
    """Detect a rectangular subset: its area and the permutation to the eckit order.

    The points of a reduced Gaussian sub-area run from the line of the area's
    north bound to the line of its south bound, without a gap, and hold on
    every line one contiguous run of longitudes, cyclically. Their indices in
    the global grid say, per line, where the run starts and how long it is;
    the bounds of the area are inferred from those runs and checked against
    eckit, as :func:`_verify_area` does.
    """
    rows = np.searchsorted(geo.ends, index, side="left")
    cols = index - geo.starts[rows]

    first, last = int(rows[0]), int(rows[-1])
    splits = np.flatnonzero(np.diff(rows)) + 1
    if last - first + 1 != splits.size + 1:
        # a line between the first and the last which holds no point: no area holds that
        return None

    bounds = np.concatenate(([0], splits, [index.size]))
    partial_starts: list = []
    partial_counts: list = []
    partial_steps: list = []
    full_steps: list = []

    for i in range(splits.size + 1):
        line_cols = cols[bounds[i] : bounds[i + 1]]
        pl = int(geo.pl[first + i])
        step = 360.0 / pl

        if line_cols.size == pl:
            # the whole line: it constrains how wide the window may be, not where it starts
            full_steps.append(step)
            continue

        gaps = np.flatnonzero(np.diff(line_cols) != 1)
        if gaps.size == 0:
            start = int(line_cols[0])
        elif gaps.size == 1 and line_cols[0] == 0 and line_cols[-1] == pl - 1:
            # the run continues across the meridian: it starts after the gap
            start = int(line_cols[gaps[0] + 1])
        else:
            # two or more gaps: a line no single longitude window can hold
            return None

        partial_starts.append(start * step)
        partial_counts.append(line_cols.size)
        partial_steps.append(step)

    if not partial_starts:
        # every line whole: a band of the grid between two lines
        return _verify_area(geo, index, first, last, 0.0, 360.0)

    ps = np.asarray(partial_starts, dtype=np.float64)
    pn = np.asarray(partial_counts, dtype=np.float64)
    pdl = np.asarray(partial_steps, dtype=np.float64)
    fds = np.asarray(full_steps, dtype=np.float64) if full_steps else None

    for anchor in np.unique(ps):
        # the runs of every partial line, unwrapped into the frame that starts
        # at the anchor: the westernmost point of a run fixes where the area
        # may start, its easternmost where it may end
        sp = anchor + (ps - anchor) % 360.0
        ep = sp + (pn - 1) * pdl

        # the east bound E of the area: at or past every run end, short of the
        # first point missing from every run; the west bound W: at or before
        # every run start, past the last point missing before it — held as
        # b = W + 360, the far end of the excluded arc (E, W + 360)
        u_lo, u_hi = float(ep.max()), float((ep + pdl).min())
        b_lo, b_hi = float((sp + 360.0 - pdl).max()), float((sp + 360.0).min())

        # whole lines add that the excluded arc holds no point of theirs at
        # all, which means it has to fit into one of their gaps: moving E up
        # to the last point of a whole line at or below the arc's far end b_lo
        # moves that line's next point, the limit on b, past b_lo; the step
        # functions are monotonic, so the smallest such E is all there is to try
        u = u_lo
        if fds is not None:
            u = max(u, float(np.max(np.floor(b_lo / fds + 1e-12) * fds)))
        if not (u < u_hi):
            continue

        b = b_hi
        if fds is not None:
            k = np.floor(u / fds + 1e-12)
            b = min(b, float(np.min((k + 1) * fds)))
        if not (b_lo < b <= u + 360.0):
            continue

        result = _verify_area(geo, index, first, last, b - 360.0, u)
        if result is not None:
            return result

    return None


def cf_reduced_gaussian_grid_name(attrs: Any, pl: Any, accumulated: bool = False) -> Optional[tuple]:
    """The global grid a CF reduced Gaussian grid mapping and its pl table name.

    Parameters
    ----------
    attrs : Any
        The CF grid mapping attributes, e.g. the attributes of the
        ``grid_mapping`` variable. Only a ``reduced_gaussian`` mapping with an
        "octahedral" or a "normal" ``grid_subtype`` and the Greenwich (absent
        or 0) first meridian names a grid eckit can build.
    pl : Any
        The number of longitude points on every latitude line of the *global*
        grid, along the dimension the ``latitude_dimension`` attribute names,
        as a 1D array. Pass the ``accumulated_points_per_latitude`` array the
        same way, with ``accumulated`` set.
    accumulated : bool
        Whether ``pl`` holds the cumulative sum of the points per latitude
        rather than those points per latitude themselves.

    Returns
    -------
    tuple or None
        ``(grid, subtype)``, e.g. ``("O32", "octahedral")``, or ``None`` when
        the metadata names no global grid: not a usable mapping, a table which
        is not the table of the grid its subtype claims. A classical table is
        checked against the table eckit builds ``N<n>`` from, which means
        generating that grid once per resolution.

    Notes
    -----
    The number of points does not come in here: the mapping and the table
    always describe the global grid, sub-area or not, and it is the
    ``reduced_gaussian_index`` a file holds which says how many of its points
    are in the file (see :func:`cf_reduced_gaussian_subset`).
    :func:`~.grid_spec.cf_reduced_gaussian_grid_spec` builds the gridspec of
    the *global* case from the same checks, alongside its index.
    """
    table = _cf_reduced_gaussian_table(attrs, pl, accumulated)
    if table is None:
        return None
    subtype, counts = table

    rows = counts.size // 2
    if subtype == "octahedral":
        return f"O{rows}", "octahedral"

    reference = _classical_pl(rows)
    if reference is None or not np.array_equal(counts, reference):
        return None

    return f"N{rows}", "normal"


def cf_reduced_gaussian_subset(name: Any, index_values: Any) -> Optional[ReducedGaussianSubset]:
    """Find what a CF ``reduced_gaussian_index`` subset of a global grid is.

    The index holds the numbers the points have on the global grid, west to
    east and north to south, and says by itself which part of that grid a file
    covers: an area, or any other subset.

    Parameters
    ----------
    name : Any
        The name of the global grid the index numbers, e.g. ``"O32"``: what
        :func:`cf_reduced_gaussian_grid_name` returns for the file's mapping.
    index_values : Any
        The ``reduced_gaussian_index`` of every point of the file, as a 1D
        array. CF requires it to be integer valued, non-negative, strictly
        increasing and below the number of points of the global grid.

    Returns
    -------
    ReducedGaussianSubset or None
        The kind of subset, with its gridspec and the permutation to the eckit
        point order when it is rectangular (a global grid is one, with the
        trivial order), or with the per-point latitudes and longitudes of the
        global grid when it is selective. ``None`` when the name is no plain
        reduced Gaussian grid or the index is not a valid subset of it.
    """
    geo = _global_geometry(name)
    if geo is None:
        return None

    index = _integer_row(index_values)
    if index is None:
        return None

    if not bool(np.all(np.diff(index) > 0)):
        return None
    if int(index[0]) < 0 or int(index[-1]) >= geo.total:
        return None

    if index.size == geo.total:
        # strictly increasing, all of them below the point count: every point of the grid
        return ReducedGaussianSubset(grid=geo.name, kind=GLOBAL, spec={"grid": geo.name})

    rectangular = _rectangular(geo, index)
    if rectangular is not None:
        area, order = rectangular
        spec = {"grid": geo.name, "area": area}
        return ReducedGaussianSubset(grid=geo.name, kind=RECTANGULAR, spec=spec, order=order)

    return ReducedGaussianSubset(grid=geo.name, kind=SELECTIVE, lats=geo.point_lats[index], lons=geo.point_lons[index])


def cf_reduced_gaussian_encoding(grid: Any) -> Optional[ReducedGaussianEncoding]:
    """The CF description of an eckit reduced Gaussian grid, to write a file with.

    Parameters
    ----------
    grid : Any
        An ``eckit.geo.Grid`` on a reduced Gaussian grid: the global grid, or
        an area of it, in any point order.

    Returns
    -------
    ReducedGaussianEncoding or None
        The global grid name and subtype, the global points per latitude and
        latitude lines, the global index of every point of ``grid`` in the
        grid's order, and the permutation which sorts those indices into the
        increasing order a CF ``reduced_gaussian_index`` must hold. ``None``
        when ``grid`` is not a plain ``O<n>``/``N<n>`` grid, with or without an
        area: a rotated or projected reduced Gaussian grid, or one with a
        custom ``pl`` table, has no CF name to carry (decision A10), and
        neither has a grid whose points are not points of the global grid.
    """
    if getattr(grid, "type", None) != "reduced_gg":
        return None

    spec = getattr(grid, "spec", None)
    if not isinstance(spec, dict) or set(spec) - _PLAIN_SPEC_KEYS:
        return None

    geo = _global_geometry(spec.get("grid"))
    if geo is None:
        return None

    lats, lons = grid.to_latlons()
    indices = _global_indices(geo, lats, lons)
    if indices is None:
        return None

    order = None if bool(np.all(np.diff(indices) > 0)) else np.argsort(indices, kind="stable")

    return ReducedGaussianEncoding(
        grid=geo.name,
        subtype=geo.subtype,
        pl=geo.pl.copy(),
        latitudes=geo.lats.copy(),
        indices=indices,
        order=order,
    )
