# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Generate the CF reduced Gaussian netCDF test files for the remote test-data server.

The server (``test-data/xr/`` on sites.ecmwf.int) holds four CF encoded reduced
Gaussian files: the real 2 m temperature sample ``t2m_O1280.nc`` and the small
global files ``t2m_O32.nc``, ``t2m_N16.nc`` and ``t2m_N80.nc``, laid out like
the sample and holding the temperatures this script first wrote into them.
Rerunning it on those fixed globals gives the same files: the sources are read
as layout templates and only their values are replaced, by the same regridded
ones. This script writes fifteen files to the output directory, ready to be
uploaded:

- the three fixed globals (same names as the originals): copies of the source
  files with only the ``air_temperature`` data replaced by the MIR regridding of
  ``t2m_O1280`` onto the grid, so every file holds valid temperature;
- twelve sub-area cuts: an area not crossing the Greenwich meridian, an area
  crossing it, a latitude band of full rows and a selective (value threshold,
  non-rectangular) subset, on each of O1280 (cut from the real sample), O32 and
  N80 (cut from the fixed globals written here). A cut keeps every dimension,
  variable, attribute, fill value and compression setting of its source, with
  only the ``reduced_gaussian_index`` subsetted; the chunk size along that
  dimension shrinks to at most the new size.

The script locates a source in the ``tests/data/xr`` test-data cache when it
holds a copy of the file, reading the cache but never writing to it; otherwise
it downloads the source into the ``sources`` subfolder of the output directory.
That is the only place it writes beside the files it generates, and it never
writes anywhere outside the output directory. It verifies everything it wrote —
headers against ``ncdump -h``, encodings against netCDF4, index order, values,
grid detection through ``earthkit.geo.grids.xarray.loader`` — prints a table of
the results and saves the md5 sums of the fifteen generated files, and only
those, to ``<output>/MD5SUMS``.

Usage, from the root of the repository::

    python tools/test-data/make_cf_reduced_gaussian_files.py [--output DIR] [--no-cache]

with ``DIR`` defaulting to ``_dev/upload/xr`` and ``--no-cache`` ignoring the
cache and downloading the sources, as a fresh machine would.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys

import netCDF4
import numpy as np
import xarray as xr
from eckit.geo import Grid

from earthkit.geo.grids.array import regrid
from earthkit.geo.grids.xarray.loader import variables as grid_variables
from earthkit.geo.grids.xarray.reduced_gaussian import (
    GLOBAL,
    RECTANGULAR,
    SELECTIVE,
    cf_reduced_gaussian_encoding,
    cf_reduced_gaussian_subset,
)
from earthkit.geo.utils.testing import earthkit_test_data_path, get_test_data_path, simple_download

CUT_DIM = "reduced_gaussian_index"
"""The only dimension a cut subsets."""

DATA_VAR = "air_temperature"
"""The data variable the files hold."""

FIXED_GLOBALS = ("O32", "N16", "N80")
"""The global files whose all-zero values are replaced by a regrid of O1280."""

CUT_GRIDS = ("O1280", "O32", "N80")
"""The grids cut into sub-area files."""

# the [north, west, south, east] areas of the rectangular cases, per grid (plan B2)
EAST_AREAS = {"O1280": [50, 10, 40, 20], "O32": [60, 10, 30, 40], "N80": [60, 10, 30, 40]}
CROSS0_AREAS = {"O1280": [55, -5, 45, 5], "O32": [60, -10, 30, 40], "N80": [60, -10, 30, 40]}
# the latitude bounds of the band of full rows, per grid
BAND_LATS = {"O1280": (44.7, 45.0), "O32": (20.0, 40.0), "N80": (20.0, 40.0)}
# the value threshold of the selective case, per grid: ("below"|"above", bound);
# O1280 selects inside its cross0 area, O32 and N80 over the whole globe
SELECTIVE_RULES = {"O1280": ("below", 282.0), "O32": ("above", 300.0), "N80": ("above", 300.0)}

MAX_SIZE_KB = 200
"""The size every sub-area file must stay below (the fixed globals keep the
exact layout of the server files they replace, so their size is that one)."""

VALUE_RANGE = (200.0, 330.0)
"""The temperatures every file holds must fall within, in K."""

_INDEX_DIM_LINE = re.compile(r"^\s*reduced_gaussian_index = \d+ ;\s*$")
"""The ``ncdump -h`` line naming the size of the index dimension."""

SOURCE_SUBFOLDER = "xr"
"""The test-data subfolder the sources live in, on the server and in the cache."""

SOURCES = ("O1280",) + FIXED_GLOBALS
"""The grids of the four source files."""


def source_path(name: str, out_dir: str, use_cache: bool = True) -> str:
    """The local path of a source file, downloading it when it is not cached.

    A copy the test suite has cached under ``tests/data/xr/`` is used as it is;
    the cache is never written to. Without a cached copy the file is downloaded
    into the ``sources`` subfolder of the output directory.

    Parameters
    ----------
    name : str
        The grid name of the source, e.g. ``"O1280"``.
    out_dir : str
        The directory the script writes to.
    use_cache : bool
        When False, the cache is ignored even when it holds the file, as on a
        fresh machine.

    Returns
    -------
    str
        The path of the local copy of the source.
    """
    file = f"t2m_{name}.nc"
    if use_cache:
        cached = os.path.join(earthkit_test_data_path(SOURCE_SUBFOLDER), file)
        if os.path.exists(cached):
            return cached
    directory = os.path.join(out_dir, "sources")
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, file)
    if not os.path.exists(path):
        simple_download(url=get_test_data_path(file, subfolder=SOURCE_SUBFOLDER), target=path)
    return path


def find_ncdump() -> str:
    """The ``ncdump`` executable to compare headers with.

    Returns
    -------
    str
        The path to ``ncdump``, from the PATH or the Homebrew location on macOS.

    Raises
    ------
    RuntimeError
        When there is no ``ncdump`` to find.
    """
    exe = shutil.which("ncdump") or "/opt/homebrew/bin/ncdump"
    if not os.path.exists(exe):
        raise RuntimeError("ncdump is needed to check the file headers: install netCDF (or put it on the PATH)")
    return exe


def header_lines(ncdump: str, path: str) -> list:
    """The ``ncdump -h`` header of a file, without the line naming it.

    Parameters
    ----------
    ncdump : str
        The ``ncdump`` executable.
    path : str
        The file to dump.

    Returns
    -------
    list of str
        The header lines after the ``netcdf <name> {`` one.
    """
    out = subprocess.run([ncdump, "-h", path], capture_output=True, text=True, check=True)
    return out.stdout.splitlines()[1:]


def check_header(ncdump: str, path: str, source_path: str, label: str, subsetted: bool) -> None:
    """Check the header of a file against the header of its source.

    Parameters
    ----------
    ncdump : str
        The ``ncdump`` executable.
    path : str
        The file to check.
    source_path : str
        The file its header must match.
    label : str
        Name of the file, for the error message.
    subsetted : bool
        Whether the file is a cut: then the line naming the size of the index
        dimension may differ, on every other line the headers must agree.

    Raises
    ------
    RuntimeError
        When a header line differs.
    """
    written, source = header_lines(ncdump, path), header_lines(ncdump, source_path)
    if subsetted:
        written = [line for line in written if not _INDEX_DIM_LINE.match(line)]
        source = [line for line in source if not _INDEX_DIM_LINE.match(line)]

    if written != source:
        difference = (
            "\n".join(f"  source: {s!r}\n  file:   {w!r}" for s, w in zip(source, written) if s != w)
            or "  (the headers have a different number of lines)"
        )
        raise RuntimeError(f"The header of {label} differs from its source:\n{difference}")


def read_field(path: str) -> np.ndarray:
    """The ``air_temperature`` field of a file, as it sits in the file.

    Parameters
    ----------
    path : str
        The file to read.

    Returns
    -------
    np.ndarray
        The values of the single field, float32, of the file's index.
    """
    with netCDF4.Dataset(path) as ds:
        return np.asarray(ds[DATA_VAR][0, 0], dtype=np.float32)


def read_index(path: str) -> np.ndarray:
    """The ``reduced_gaussian_index`` coordinate of a file.

    Parameters
    ----------
    path : str
        The file to read.

    Returns
    -------
    np.ndarray
        The global index of every point of the file, int64, in the order of the file.
    """
    with netCDF4.Dataset(path) as ds:
        return np.asarray(ds[CUT_DIM][:], dtype=np.int64)


def check_values(label: str, values: np.ndarray) -> None:
    """Check that the values of a file are finite temperatures.

    Parameters
    ----------
    label : str
        Name of the file, for the error message.
    values : np.ndarray
        The values read back from the file.

    Raises
    ------
    RuntimeError
        When a value is not finite or outside the accepted range.
    """
    if not bool(np.all(np.isfinite(values))):
        raise RuntimeError(f"{label}: {int(np.count_nonzero(~np.isfinite(values)))} non-finite values")
    lo, hi = VALUE_RANGE
    if bool(values.size) and not (lo <= float(values.min()) and float(values.max()) <= hi):
        raise RuntimeError(f"{label}: values {float(values.min()):.2f}-{float(values.max()):.2f} K outside {lo}-{hi} K")


def row_geometry(name: str) -> tuple:
    """The latitude lines of a global grid, as eckit builds it.

    Parameters
    ----------
    name : str
        The grid name, e.g. ``"O32"``.

    Returns
    -------
    tuple
        ``(row_lats, starts, sizes)``: the latitude of every line north to
        south, the global index of the first point of every line and the number
        of points on every line, as float64 and int64 arrays.
    """
    lats = np.asarray(Grid({"grid": name}).to_latlons()[0], dtype=np.float64)
    bounds = np.concatenate(([0], np.flatnonzero(np.diff(lats) != 0) + 1, [lats.size]))
    return lats[bounds[:-1]], bounds[:-1], np.diff(bounds)


def check_increasing(label: str, index: np.ndarray) -> None:
    """Check that an index is a strictly increasing run of global point numbers.

    Parameters
    ----------
    label : str
        Name of the file, for the error message.
    index : np.ndarray
        The index to check.

    Raises
    ------
    RuntimeError
        When the index is empty, not strictly increasing, or negative.
    """
    if index.size == 0:
        raise RuntimeError(f"{label}: the index subset is empty")
    if bool(np.any(np.diff(index) <= 0)) or int(index[0]) < 0:
        raise RuntimeError(f"{label}: the index is not strictly increasing and non-negative")


def rectangular_indices(name: str, area: list, label: str) -> np.ndarray:
    """The CF index of the points of an eckit area grid, in increasing order.

    Parameters
    ----------
    name : str
        The global grid the area sits on.
    area : list
        The ``[north, west, south, east]`` bounds of the area.
    label : str
        Name of the file, for the error message.

    Returns
    -------
    np.ndarray
        The global index of every point of the area, int64, sorted increasing.

    Raises
    ------
    RuntimeError
        When eckit or the CF encoder rejects the area.
    """
    encoding = cf_reduced_gaussian_encoding(Grid({"grid": name, "area": list(area)}))
    if encoding is None:
        raise RuntimeError(f"{label}: eckit area {area} of {name} is not encodable as a CF reduced Gaussian grid")
    index = np.sort(np.asarray(encoding.indices, dtype=np.int64))
    check_increasing(label, index)
    subset = cf_reduced_gaussian_subset(name, index)
    if subset is None or subset.kind != RECTANGULAR:
        kind = None if subset is None else subset.kind
        raise RuntimeError(f"{label}: the points of area {area} are detected as {kind!r}, not rectangular")
    return index


def band_indices(name: str, lo: float, hi: float, label: str) -> np.ndarray:
    """The CF index of the full latitude lines of a grid between two latitudes.

    The index is taken from the eckit area grid with bounds ``[hi, 0, lo, 360]``
    (plan B6); when that grid does not hold exactly the full lines between the
    bounds, the lines are selected from the global grid itself.

    Parameters
    ----------
    name : str
        The global grid.
    lo : float
        The southern latitude bound.
    hi : float
        The northern latitude bound.
    label : str
        Name of the file, for the error message.

    Returns
    -------
    np.ndarray
        The global index of every point of the band, int64, sorted increasing.

    Raises
    ------
    RuntimeError
        When the band holds no line.
    """
    row_lats, starts, sizes = row_geometry(name)
    rows = np.flatnonzero((row_lats >= lo) & (row_lats <= hi))
    if rows.size == 0:
        raise RuntimeError(f"{label}: no latitude line of {name} between {lo} and {hi}")
    full = np.sort(np.concatenate([np.arange(starts[r], starts[r] + sizes[r], dtype=np.int64) for r in rows]))

    try:
        index = rectangular_indices(name, [hi, 0.0, lo, 360.0], label)
    except RuntimeError:
        index = np.empty(0, dtype=np.int64)
    if not np.array_equal(index, full):
        print(f"note: {label}: the area grid does not hold the full lines; selecting {full.size} points from pl")
        index = full
    check_increasing(label, index)
    return index


def selective_indices(name: str, values: np.ndarray, within: np.ndarray | None, label: str) -> np.ndarray:
    """The CF index of the points above or below a value threshold.

    Parameters
    ----------
    name : str
        The global grid the values sit on.
    values : np.ndarray
        The field the threshold is applied to, in the global index order.
    within : np.ndarray or None
        When given, the global index of the points the selection is restricted
        to; when ``None``, the whole globe is searched.
    label : str
        Name of the file, for the error message.

    Returns
    -------
    np.ndarray
        The global index of the selected points, int64, sorted increasing.

    Raises
    ------
    RuntimeError
        When nothing is selected or the selection is detected as rectangular.
    """
    rule = SELECTIVE_RULES[name]
    if within is None:
        candidates = np.arange(values.size, dtype=np.int64)
        keep = values > rule[1] if rule[0] == "above" else values < rule[1]
    else:
        candidates = within
        field = values[within]
        keep = field > rule[1] if rule[0] == "above" else field < rule[1]
    index = np.sort(candidates[keep])

    check_increasing(label, index)
    subset = cf_reduced_gaussian_subset(name, index)
    if subset is None or subset.kind != SELECTIVE:
        kind = None if subset is None else subset.kind
        raise RuntimeError(
            f"{label}: the {index.size} points {rule[0]} {rule[1]} K are detected as {kind!r}, "
            "not a selective subset: tune the threshold"
        )
    return index


def create_variable(ds: netCDF4.Dataset, source: netCDF4.Dataset, name: str, size: int) -> netCDF4.Variable:
    """Create a variable of a cut file, copying the definition of its source.

    Parameters
    ----------
    ds : netCDF4.Dataset
        The file being written.
    source : netCDF4.Dataset
        The file it is cut from.
    name : str
        The name of the variable.
    size : int
        The new size of the index dimension.

    Returns
    -------
    netCDF4.Variable
        The created variable, with the data type, dimensions, attributes and
        compression of the source; the chunk size along the index dimension
        shrinks to at most the new size, everywhere else the chunking is kept.
    """
    var = source[name]
    dims = list(var.dimensions)

    kwargs = {}
    chunks = var.chunking()
    if isinstance(chunks, list):
        kwargs["chunksizes"] = [min(c, size) if d == CUT_DIM else c for c, d in zip(chunks, dims)]
    filters = var.filters()
    if filters.get("zlib"):
        kwargs.update(zlib=True, complevel=filters["complevel"], shuffle=filters["shuffle"])

    new = ds.createVariable(name, var.dtype, dims, **kwargs)
    new.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs()})
    return new


def write_cut(source_path: str, out_path: str, index: np.ndarray) -> None:
    """Write a sub-area cut of a CF reduced Gaussian file.

    Parameters
    ----------
    source_path : str
        The file to cut.
    out_path : str
        The file to write.
    index : np.ndarray
        The global index of the points to keep, int64, strictly increasing.
    """
    with netCDF4.Dataset(source_path) as source:
        with netCDF4.Dataset(out_path, "w", format=source.file_format) as ds:
            for attr in source.ncattrs():
                ds.setncattr(attr, source.getncattr(attr))
            for dim in source.dimensions:
                if source.dimensions[dim].isunlimited():
                    ds.createDimension(dim, None)
                else:
                    size = index.size if dim == CUT_DIM else source.dimensions[dim].size
                    ds.createDimension(dim, size)

            for name in source.variables:
                var = source.variables[name]
                new = create_variable(ds, source, name, index.size)
                if var.ndim == 0:
                    new[...] = var[...]
                elif CUT_DIM in var.dimensions:
                    axis = var.dimensions.index(CUT_DIM)
                    new[...] = np.take(np.asarray(var[...]), index, axis=axis)
                else:
                    new[...] = np.asarray(var[...])


def _attr_equal(a, b) -> bool:
    """Whether two netCDF attribute values are equal, NaNs counting as equal."""
    try:
        if float(a) != float(a) and float(b) != float(b):
            return True
    except (TypeError, ValueError):
        pass
    return type(a) is type(b) and bool(a == b)


def check_encoding(label: str, cut_path: str, source_path: str) -> None:
    """Check the encoding of a cut against the encoding of its source.

    Every variable must keep the data type, dimensions, attributes, fill value
    and compression of its source; the chunk size must stay equal, except along
    the index dimension, where it must be at most the new size.

    Parameters
    ----------
    label : str
        Name of the cut, for the error message.
    cut_path : str
        The cut file to check.
    source_path : str
        The file it was cut from.

    Raises
    ------
    RuntimeError
        When any encoding detail differs.
    """
    with netCDF4.Dataset(source_path) as src, netCDF4.Dataset(cut_path) as dst:
        if list(src.variables) != list(dst.variables):
            raise RuntimeError(f"{label}: the variables differ from the source")
        for name, var in src.variables.items():
            new = dst[name]
            problems = []
            if var.dtype != new.dtype:
                problems.append(f"dtype {var.dtype} != {new.dtype}")
            if var.dimensions != new.dimensions:
                problems.append(f"dims {var.dimensions} != {new.dimensions}")
            attrs, new_attrs = var.ncattrs(), new.ncattrs()
            if attrs != new_attrs or not all(_attr_equal(var.getncattr(a), new.getncattr(a)) for a in attrs):
                problems.append(f"attributes {attrs} != {new_attrs}")
            if var.filters() != new.filters():
                problems.append(f"compression {var.filters()} != {new.filters()}")

            chunks, new_chunks = var.chunking(), new.chunking()
            if isinstance(chunks, list) != isinstance(new_chunks, list):
                problems.append(f"chunking {chunks!r} != {new_chunks!r}")
            elif isinstance(chunks, list):
                for axis, (d, c, n) in enumerate(zip(var.dimensions, chunks, new_chunks)):
                    if d == CUT_DIM:
                        if n != min(c, new.shape[axis]):
                            problems.append(f"chunk along {d}: {n} != min({c}, {new.shape[axis]})")
                    elif n != c:
                        problems.append(f"chunk along {d}: {c} != {n}")

            if problems:
                raise RuntimeError(f"{label}: variable {name}: " + "; ".join(problems))


def detect(path: str) -> dict:
    """The grid the earthkit-geo xarray reader detects in a file.

    Parameters
    ----------
    path : str
        The file to read.

    Returns
    -------
    dict
        The name of the data variable, the detected ``grid_spec``, the type of
        the eckit grid and whether a non-trivial ``point_order`` was detected.

    Raises
    ------
    RuntimeError
        When the file holds no single geographical variable.
    """
    with xr.open_dataset(path) as ds:
        variables = grid_variables(ds)
        if len(variables) != 1 or variables[0].name != DATA_VAR:
            raise RuntimeError(f"{path}: detected {len(variables)} variables, expected one named {DATA_VAR!r}")
        v = variables[0]
        return {
            "spec": v.xr_grid.grid_spec,
            "ek_type": v.ek_grid.type,
            "point_order": v.point_order is not None,
        }


def check_detection(label: str, path: str, name: str, kind: str) -> dict:
    """Check that the reader detects the grid a file is meant to hold.

    Parameters
    ----------
    label : str
        Name of the file, for the error message.
    path : str
        The file to read back.
    name : str
        The global grid the file sits on.
    kind : str
        One of ``"global"``, ``"east"``, ``"cross0"``, ``"band"``, ``"selective"``.

    Returns
    -------
    dict
        What :func:`detect` found.

    Raises
    ------
    RuntimeError
        When the detected grid is not the expected one.
    """
    found = detect(path)
    spec = found["spec"]

    if kind == "global":
        ok = spec == {"grid": name} and found["ek_type"] == "reduced_gg" and not found["point_order"]
    elif kind in ("east", "band"):
        ok = (
            isinstance(spec, dict)
            and spec.get("grid") == name
            and "area" in spec
            and found["ek_type"] == "reduced_gg"
            and not found["point_order"]
        )
    elif kind == "cross0":
        ok = (
            isinstance(spec, dict)
            and spec.get("grid") == name
            and "area" in spec
            and found["ek_type"] == "reduced_gg"
            and found["point_order"]
        )
    else:  # selective: the points are a scatter the reader gives to eckit as unstructured
        ok = spec is None and found["ek_type"].startswith("unstructured")

    if not ok:
        raise RuntimeError(f"{label}: {kind} subset of {name} detected as {found!r}")
    return found


def size_kb(path: str) -> float:
    """The size of a file in kilobytes."""
    return os.path.getsize(path) / 1024.0


def md5(path: str) -> str:
    """The md5 sum of a file, as a hex string."""
    digest = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class Report:
    """The table of the files written, with a line per file."""

    def __init__(self):
        self.rows = []

    def add(self, path: str, name: str, kind: str, points: int, found: dict) -> None:
        """Register one written file, checking the size of a sub-area file.

        Parameters
        ----------
        path : str
            The file written.
        name : str
            The global grid the file sits on.
        kind : str
            What the file holds: ``"global"`` or one of the four sub-area cases.
        points : int
            The number of points the file holds.
        found : dict
            What :func:`check_detection` found in it.
        """
        size = size_kb(path)
        if kind != "global" and size > MAX_SIZE_KB:
            raise RuntimeError(f"{os.path.basename(path)}: {size:.1f} KB, above the {MAX_SIZE_KB} KB limit")
        spec = found["spec"] if found["spec"] is not None else "unstructured"
        self.rows.append((os.path.basename(path), name, kind, f"{points:d}", f"{size:.1f}", str(spec), md5(path)))

    def print(self, out_dir: str) -> None:
        """Print the table.

        Parameters
        ----------
        out_dir : str
            The directory the files were written to, named in the header.
        """
        header = ("file", "grid", "kind", "points", "size KB", "detected spec", "md5")
        columns = list(zip(*(list(header) + list(row) for row in self.rows)))
        widths = [max(len(cell) for cell in column) for column in columns]

        print(f"\n{len(self.rows)} files written to {out_dir}:\n")
        for row in (header, *self.rows):
            print("  ".join(cell.ljust(width) for cell, width in zip(row, widths)).rstrip())


def write_fixed_globals(out_dir: str, sources: dict) -> dict:
    """Write the copies of the zero global files holding regridded temperatures.

    Parameters
    ----------
    out_dir : str
        The directory to write to.
    sources : dict
        The local paths of the source files, by name.

    Returns
    -------
    dict
        The paths written, by grid name.
    """
    o1280 = read_field(sources["O1280"])
    fixed = {}
    for name in FIXED_GLOBALS:
        source = sources[name]
        out_path = os.path.join(out_dir, os.path.basename(source))
        label = os.path.basename(out_path)

        values, spec = regrid(o1280, in_grid={"grid": "O1280"}, out_grid={"grid": name}, interpolation="linear")
        if spec != {"grid": name}:
            raise RuntimeError(f"{label}: the regrid returned grid spec {spec!r}, expected {{'grid': {name!r}}}")
        values32 = np.asarray(values, dtype=np.float32)
        if values32.shape != (int(Grid({"grid": name}).size()),):
            raise RuntimeError(f"{label}: the regrid returned {values32.shape}, expected one value per global point")

        shutil.copyfile(source, out_path)
        with netCDF4.Dataset(out_path, "a") as ds:
            ds[DATA_VAR][0, 0] = values32

        written = read_field(out_path)
        np.testing.assert_allclose(written, values, rtol=1e-6)
        check_values(label, written)
        fixed[name] = out_path
    return fixed


def write_cuts(out_dir: str, sources: dict, fixed: dict, ncdump: str, report: Report) -> None:
    """Write and check the twelve sub-area files.

    Parameters
    ----------
    out_dir : str
        The directory to write to.
    sources : dict
        The local paths of the source files, by name.
    fixed : dict
        The paths of the fixed global files written before, by name.
    ncdump : str
        The ``ncdump`` executable.
    report : Report
        Collects the line of every file written.
    """
    for name in CUT_GRIDS:
        source_path = sources[name] if name == "O1280" else fixed[name]
        field = read_field(source_path)

        rectangles = {
            "east": rectangular_indices(name, EAST_AREAS[name], f"t2m_{name}_area_east.nc"),
            "cross0": rectangular_indices(name, CROSS0_AREAS[name], f"t2m_{name}_area_cross0.nc"),
            "band": band_indices(name, *BAND_LATS[name], f"t2m_{name}_band.nc"),
        }
        within = rectangles["cross0"] if name == "O1280" else None
        rectangles["selective"] = selective_indices(name, field, within, f"t2m_{name}_selective.nc")

        for kind, index in rectangles.items():
            out_path = os.path.join(out_dir, f"t2m_{name}_{'area_' if kind in ('east', 'cross0') else ''}{kind}.nc")
            label = os.path.basename(out_path)

            write_cut(source_path, out_path, index)

            check_header(ncdump, out_path, source_path, label, subsetted=True)
            check_encoding(label, out_path, source_path)
            written_index, written_field = read_index(out_path), read_field(out_path)
            if not np.array_equal(written_index, index):
                raise RuntimeError(f"{label}: the written index differs from the requested one")
            check_increasing(label, written_index)
            if not np.array_equal(written_field, field[index]):
                raise RuntimeError(f"{label}: the written values differ from the source values at the index")
            check_values(label, written_field)
            found = check_detection(label, out_path, name, kind)
            report.add(out_path, name, kind, int(index.size), found)


def write_globals_report(sources: dict, fixed: dict, ncdump: str, report: Report) -> None:
    """Check and report the three fixed global files.

    Parameters
    ----------
    sources : dict
        The local paths of the source files, by name.
    fixed : dict
        The paths of the fixed global files, by name.
    ncdump : str
        The ``ncdump`` executable.
    report : Report
        Collects the line of every file written.
    """
    for name in FIXED_GLOBALS:
        out_path = fixed[name]
        label = os.path.basename(out_path)
        check_header(ncdump, out_path, sources[name], label, subsetted=False)
        field, index = read_field(out_path), read_index(out_path)
        subset = cf_reduced_gaussian_subset(name, index)
        if subset is None or subset.kind != GLOBAL:
            kind = None if subset is None else subset.kind
            raise RuntimeError(f"{label}: the global index is detected as {kind!r}, not global")
        check_values(label, field)
        found = check_detection(label, out_path, name, "global")
        report.add(out_path, name, "global", int(field.size), found)


def main(argv=None) -> int:
    """Run the generator.

    Parameters
    ----------
    argv : list of str, optional
        The command line arguments, by default ``sys.argv[1:]``.

    Returns
    -------
    int
        Zero, when every file was written and verified.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default=os.path.join("_dev", "upload", "xr"), help="where to write the files")
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="ignore the test-data cache and download the sources into <output>/sources",
    )
    args = parser.parse_args(argv)

    out_dir = args.output
    os.makedirs(out_dir, exist_ok=True)
    ncdump = find_ncdump()

    sources = {}
    for name in SOURCES:
        sources[name] = source_path(name, out_dir, use_cache=not args.no_cache)
    print("sources:")
    for name in SOURCES:
        print(f"  t2m_{name}.nc: {sources[name]}")

    report = Report()
    fixed = write_fixed_globals(out_dir, sources)
    write_globals_report(sources, fixed, ncdump, report)
    write_cuts(out_dir, sources, fixed, ncdump, report)

    sums = os.path.join(out_dir, "MD5SUMS")
    with open(sums, "w") as f:
        for row in sorted(report.rows):
            f.write(f"{row[6]}  {row[0]}\n")
    report.print(out_dir)
    print(f"\nmd5 sums written to {sums}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
