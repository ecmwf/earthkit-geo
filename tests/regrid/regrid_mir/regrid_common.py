# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Shared helpers for the regrid tests in ``tests/regrid/regrid_mir/``.

This module is not a test file. It holds everything the regrid tests share --
the interpolation registry, the kwarg alias list, the skip marks and the
data-access functions -- so that ``test_regrid_numpy.py``, ``test_mir_grid.py``,
``test_regrid_xarray.py``, ``test_regrid_fieldlist_grib.py``,
``test_regrid_fieldlist_array.py``, ``test_regrid_grib_message.py`` and the matrix tests
do not duplicate them.
The folder has no ``__init__.py``: pytest puts it on ``sys.path``, and the test
files import it with ``from regrid_common import ...``. The module has no side
effects at import: nothing is downloaded, and ``mir`` is only ever touched
through the ``NO_MIR`` flag of ``earthkit.geo.utils.testing``.

Public interface (fixed by plan task T1, used by T1b and T3-T6; the comparison
rules of the two bullets below are the revised ones of plan task T8)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``INTERPOLATIONS``
    ``dict`` mapping the canonical name to an :class:`Interpolation` (a frozen
    dataclass with the fields ``name``, ``aliases`` and ``compare``). The
    registry holds "nearest-neighbour" (aliases "nn", "nearest-neighbor"),
    "linear" and "grid-box-average".
    ``compare(result, reference, skip_pole_rows=False, **kwargs)`` raises
    ``AssertionError`` when the regridded ``result`` does not match the
    ``reference``; ``kwargs`` (e.g. ``rtol``, ``atol``) are forwarded to the
    underlying numpy comparison and every rule takes the ``skip_pole_rows``
    keyword (the matrix test passes it for its two ``assert_allclose`` rules).
    The rules:

    * **nearest-neighbour** is *strict*: NaN at exactly the same positions in
      both arrays, and every finite value exactly equal (``rtol=0``, ``atol=0``
      by default). No pole-row exemption and no share of values allowed to
      differ -- a systematic tie is machine dependent and can change *every*
      value, which a proportion rule cannot tell apart from a wrong result.
      Where the strict comparison fails, the caller runs the independent tie
      validation of :func:`assert_nearest_ties` on the differing positions;
      that is what makes a legitimate equidistant choice on another machine
      pass while a value from a farther source point still fails.
    * **linear** and **grid-box-average** compare the whole array with
      ``assert_allclose`` (whose ``equal_nan=True`` default already fails on a
      differing NaN mask). The pole rows of a global regular lat-lon result are
      machine dependent, so both take ``skip_pole_rows=False``; the caller sets
      it to True only for a result it knows to be a global regular lat-lon grid
      (the matrix test: target "ll-10x10" and the pair not clipped), and then
      rows 0 and -1 are left out of the comparison -- nothing else ever is.
``get_interpolation(name)``
    The :class:`Interpolation` registered under the canonical name *or* any
    alias; raises ``ValueError`` for an unknown name.
``assert_nearest_ties(result, reference, input_values, source_spec, target_spec)``
    Independent validation of the positions where a strict nearest comparison
    found a difference: each of them must hold the input value of a source
    point that is as near to the target point as the nearest source point is
    (within :data:`TIE_REL_TOL` relative and :data:`TIE_ABS_TOL` absolute
    tolerance, measured as a 3D chord distance on the unit sphere). Raises
    ``AssertionError`` listing the positions that are not ties. Either grid
    argument may also be given as an ``(lats, lons)`` pair of coordinate
    arrays, in the flat point order of the values it describes (a tuple of
    two arrays is taken as coordinates, anything else as a gridspec); that
    is the mode for a grid without a usable gridspec, e.g. the projection
    grids of a GRIB fieldlist.
``BASE_INTERPOLATION_NAMES`` / ``INTERPOLATION_NAMES``
    The parametrisation lists the old test modules declared as
    ``BASE_INTERPOLATIONS`` (["linear", "nearest-neighbour"]) and
    ``INTERPOLATIONS`` (that plus "grid-box-average"); renamed with a ``_NAMES``
    suffix because ``INTERPOLATIONS`` is the registry dict here.
``KWARG_ALIASES``
    The list of ``(kwargs, canonical_name)`` tuples of the old
    ``test_regrid_numpy_kwarg`` / ``test_regrid_fieldlist_reg_ll``
    parametrisations, e.g. ``({}, "linear")`` and
    ``({"interpolation": "nn"}, "nearest-neighbour")``.
``requires_mir`` / ``requires_ekd``
    ``pytest.mark.skipif(NO_MIR, reason="No mir available")`` and
    ``pytest.mark.skipif(NO_EKD, reason="No earthkit-data available")``; the
    flags ``NO_MIR``, ``NO_EKD`` and ``IN_GITHUB`` are re-exported for tests
    that compose their own skipif conditions.
``data_path(kind, *parts, source=..., subfolder=...)``
    Local :class:`pathlib.Path` of a test data file, in one place. New
    hierarchy: ``data_path("input", source_id)``, ``data_path("reference",
    interpolation, source_id, target_id)`` and ``data_path("grib_reference",
    interpolation, filename, target_id)`` -- the reference of the GRIB matrix
    of ``regrid_grib_cases.py``, one level deeper, at
    ``<root>/regrid/grib/<interpolation>/<filename>/<target>.npz`` -- all
    return the file of the cache :data:`CACHE_ROOT`, downloading it from the
    server when it is missing (the interpolation part is canonicalised, so an
    alias names the same file; "grib" is not an interpolation name, so a
    "reference" and a "grib_reference" path can never name the same file).
    With ``EARTHKIT_GEO_REGRID_TEST_DATA`` set, the same parts resolve under
    that local tree instead and nothing is ever downloaded (see "Where the new
    data files come from" below). Legacy files: ``data_path("legacy",
    filename, subfolder="global_0_360")`` delegates to
    ``earthkit.geo.utils.testing.get_test_data``.
``load_values(path)`` / ``save_values(path, values)``
    Read / write the ``arr_0`` key of an npz file (the layout of all the regrid
    data files, old and new); ``save_values`` creates parent folders.
``legacy_values(filename, subfolder="global_0_360")``
    Convenience: the array of a legacy data file.
``CACHE_ROOT`` / ``REGRID_TEST_DATA_ENV`` / ``GENERATOR_COMMAND`` /
``_fetch_new(relative)`` / ``_download(url, target)``
    The cache folder of the downloaded data (``tests/data``, git-ignored), the
    name of the environment variable that points the new-hierarchy files (not
    the legacy ones) at a local tree, the generator command quoted in the error message of a missing file, and the
    two functions that implement the lookup: :func:`_fetch_new` is the single
    place where a new-hierarchy file is located (override tree, cache,
    download), :func:`_download` the one request the tests fake.

Adding an interpolation
^^^^^^^^^^^^^^^^^^^^^^^

1. Write a ``_compare_<name>(result, reference, skip_pole_rows=False, **kwargs)`` function
   that raises ``AssertionError`` on mismatch (see the existing ones). Forward ``kwargs`` to
   the numpy comparison and handle the ``skip_pole_rows`` flag: an ``assert_allclose`` rule
   leaves rows 0 and -1 out of the comparison when it is True (and only then -- the flag is
   never derived from the shape of the array), a strict rule that copies values refuses it,
   because its machine-dependent differences are ties and belong in
   :func:`assert_nearest_ties`, not in an exemption.
2. Add ``Interpolation("<name>", (<aliases, or ()>), _compare_<name>)`` to the
   ``INTERPOLATIONS`` dict below; the dict key is the canonical name, which is
   also the folder name of the reference files.
3. Add the canonical name to ``INTERPOLATION_NAMES`` (or to a selection list in
   ``regrid_cases.py``) where the new interpolation should be tested, and list
   its kwargs spellings in ``KWARG_ALIASES`` if a kwarg-alias test should cover
   it. Reference files are generated by ``generate_regrid_data.py``.

Where the new data files come from
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:data:`data_path` resolves the new-hierarchy files through one function,
:func:`_fetch_new`, which turns a relative path (``regrid/input/<id>.npz``,
``regrid/<interpolation>/<source>/<target>.npz`` or, for the GRIB matrix,
``regrid/grib/<interpolation>/<file>/<target>.npz``) into a local file:

* Local override: with the environment variable named by
  :data:`REGRID_TEST_DATA_ENV` (``EARTHKIT_GEO_REGRID_TEST_DATA``) set to a
  folder that contains the ``regrid/`` tree, the file is read from that folder
  and the network is never touched; a missing file raises ``FileNotFoundError``
  naming :data:`GENERATOR_COMMAND`. This is how freshly generated data (the
  generator's output, ``_dev/test-data`` by default or its ``--output-dir``)
  are tested before they are uploaded.
* Cache and server: without the variable, the file lives under
  :data:`CACHE_ROOT` (``tests/data/regrid/``, git-ignored), and is downloaded
  from ``URL_ROOT/regrid/...`` (the ECMWF download server, which holds the
  uploaded tree) when it is not cached yet. The download goes to a temporary
  file in the same folder and is renamed atomically, so a broken download
  leaves no partial file; a failure raises ``FileNotFoundError`` naming the
  URL. Deleting the cache is always safe: the next run refetches what it needs.

No test and no path in ``regrid_cases.py`` knows the roots; everything about
where the data live is in :func:`_fetch_new` and the constants above.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pytest

from earthkit.geo.utils.testing import (
    IN_GITHUB,  # noqa: F401  # re-exported for the refactored test modules
    NO_EKD,  # noqa: F401  # re-exported for the refactored test modules
    NO_MIR,  # noqa: F401  # re-exported for the refactored test modules
    URL_ROOT,  # the root the new data tree is downloaded from (used by _fetch_new)
    get_test_data,
)
from earthkit.geo.utils.url import join_url_path

__all__ = [
    "BASE_INTERPOLATION_NAMES",
    "CACHE_ROOT",
    "GENERATOR_COMMAND",
    "INTERPOLATIONS",
    "INTERPOLATION_NAMES",
    "IN_GITHUB",
    "KWARG_ALIASES",
    "NO_EKD",
    "NO_MIR",
    "REGRID_TEST_DATA_ENV",
    "TIE_ABS_TOL",
    "TIE_MAX_NEIGHBOURS",
    "TIE_REL_TOL",
    "Interpolation",
    "assert_nearest_ties",
    "data_path",
    "get_interpolation",
    "legacy_values",
    "load_values",
    "requires_ekd",
    "requires_mir",
    "save_values",
]

# The lists the old test modules called BASE_INTERPOLATIONS and INTERPOLATIONS.
BASE_INTERPOLATION_NAMES: list[str] = ["linear", "nearest-neighbour"]
INTERPOLATION_NAMES: list[str] = BASE_INTERPOLATION_NAMES + ["grid-box-average"]

# The (kwargs, canonical interpolation name) pairs of the kwarg-alias tests,
# verifying that the default and every MIR spelling of an interpolation name
# selects that interpolation.
KWARG_ALIASES: list[tuple[dict, str]] = [
    ({}, "linear"),
    ({"interpolation": "linear"}, "linear"),
    ({"interpolation": "nearest-neighbour"}, "nearest-neighbour"),
    ({"interpolation": "nn"}, "nearest-neighbour"),
    ({"interpolation": "nearest-neighbor"}, "nearest-neighbour"),
]

# Skip marks with one wording each, replacing the per-file skipif decorators.
requires_mir = pytest.mark.skipif(NO_MIR, reason="No mir available")
requires_ekd = pytest.mark.skipif(NO_EKD, reason="No earthkit-data available")

# The new data hierarchy is the same relative tree at every root (the server, the
# override folder of EARTHKIT_GEO_REGRID_TEST_DATA and this cache):
#   <root>/regrid/input/<source>.npz
#   <root>/regrid/<interpolation>/<source>/<target>.npz
#   <root>/regrid/grib/<interpolation>/<grib file>/<target>.npz   (the GRIB matrix, data_path kind "grib_reference")
# CACHE_ROOT is the local cache the files are downloaded into (tests/data/,
# git-ignored: deleting it is always safe, the tests refetch what they need).
# The server side is URL_ROOT/regrid/..., uploaded from the generator's output
# tree by the maintainer; _fetch_new() resolves a relative path on either root.
CACHE_ROOT = Path(__file__).resolve().parents[2] / "data"

# Name of the environment variable that overrides cache and network: the folder
# of a local regrid/ tree (the generator writes _dev/test-data/regrid/ by default
# — or its --output-dir — and sets the variable to that folder itself). data_path
# then reads from that folder and never downloads; see the module docstring,
# "Where the new data files come from".
REGRID_TEST_DATA_ENV = "EARTHKIT_GEO_REGRID_TEST_DATA"

# Timeout of one data download request, in seconds (as in
# earthkit.geo.utils.testing.simple_download: without it a stalled server
# hangs the request, and the test, forever instead of failing).
DOWNLOAD_TIMEOUT = 30

# Quoted in the error message of a missing new data file; the generator is
# written by plan task T1c and the "..." stands for its --source/--target/
# --interpolation arguments.
GENERATOR_COMMAND = "python tests/regrid/regrid_mir/generate_regrid_data.py ..."

# The subfolder of the legacy test data files (the default of get_test_data).
DEFAULT_LEGACY_SUBFOLDER = "global_0_360"

# A nearest-neighbour result may legitimately differ from its reference where two
# source points are equally near to the target point: a tie candidate is accepted
# when its distance to the target point is within TIE_REL_TOL relative and
# TIE_ABS_TOL absolute tolerance of the distance of the nearest source point. The
# distances are 3D chords on the unit sphere, so the absolute tolerance is in
# chords (1e-12 is a few tens of micrometres on the Earth); the relative one
# covers the rounding of the coordinates themselves.
TIE_REL_TOL = 1e-9
TIE_ABS_TOL = 1e-12

# How many neighbours the KD-tree returns per point when the tie candidates are
# looked up. 16 is far more than a regular grid can produce at one target point;
# positions whose window is still full of equal distances (the pole rows of a
# global lat-lon source, where a whole row is at the same point) are re-queried
# without a limit, see assert_nearest_ties().
TIE_MAX_NEIGHBOURS = 16


@dataclass(frozen=True)
class Interpolation:
    """One regridding interpolation: canonical name, aliases and compare rule.

    Attributes
    ----------
    name : str
        Canonical name; the key in :data:`INTERPOLATIONS` and the folder name
        of the reference files.
    aliases : tuple of str
        Other names that select this interpolation.
    compare : callable
        ``compare(result, reference, **kwargs)`` raising ``AssertionError``
        when the regridded ``result`` does not match the ``reference``. The
        kwargs (e.g. ``rtol``, ``atol``) are forwarded to the numpy comparison;
        every rule also takes ``skip_pole_rows=False``, which leaves the pole
        rows of a global regular lat-lon result out of the ``assert_allclose``
        comparison (see :func:`_compare_linear`) -- nearest-neighbour has no
        such exemption and refuses it (see :func:`_compare_nearest`).
    """

    name: str
    aliases: tuple[str, ...]
    compare: Callable


def _prepare(result, reference):
    """Return result and reference as numpy arrays with matching shapes."""
    result = np.asarray(result)
    reference = np.asarray(reference)
    assert reference.size == result.size, f"reference has {reference.size} values, result has {result.size}"
    return result, reference.reshape(result.shape)


def _assert_same_nan_mask(result, reference):
    """Assert NaN positions of ``result`` and ``reference`` are identical; return finite mask.

    MIR NaN-fills a nearest-neighbour result outside a regional source instead of clipping
    (finding F5 of ``.ecmwf-multi-agents/notes/regrid-matrix-probe.md``), and ``np.isclose``
    is False for a pair of NaNs, so the NaN positions are checked here once and the numeric
    rules then run over the finite values only. An array without NaNs returns an all-True
    mask.
    """
    if result.size == 0:
        return np.zeros(result.shape, dtype=bool)
    result_nan = np.isnan(result)
    reference_nan = np.isnan(reference)
    assert np.array_equal(result_nan, reference_nan), (
        f"NaN positions differ: {int((result_nan ^ reference_nan).sum())} of {result.size} values "
        f"are NaN in one array but not in the other"
    )
    return ~result_nan


def _compare_nearest(result, reference, skip_pole_rows=False, **kwargs):
    """Compare a nearest-neighbour result with its reference, strictly.

    ``skip_pole_rows`` is only there to be refused: this rule has no exemption at all, its
    machine-dependent differences are ties and they are validated by
    :func:`assert_nearest_ties`, not excused.

    Result and reference must have NaN at exactly the same positions (the target points
    outside a regional source are NaN-filled, see :func:`_assert_same_nan_mask`), and then
    every finite value must be EQUAL: the default is ``rtol=0`` and ``atol=0``, and there is
    no pole-row exemption and no share of values allowed to differ.

    A nearest-neighbour tie is machine dependent: several source points at the same distance
    from a target point (a target point exactly between two cell centres of a reference-grid
    lat-lon source, or a whole source row at a pole) make MIR pick one of them, and another
    build of MIR may pick another one -- systematically, in every value at once. A proportion
    rule cannot tell that apart from a wrong result, so the difference is handled where it can
    be checked: the caller catches this ``AssertionError`` and runs :func:`assert_nearest_ties`
    over the differing positions, which accepts a value only if it is the input value of a
    source point as near to the target point as the nearest one is.
    """
    assert not skip_pole_rows, (
        "the nearest-neighbour comparison is strict and has no pole-row exemption; a tie is "
        "accepted by regrid_common.assert_nearest_ties, not by skipping rows"
    )
    result, reference = _prepare(result, reference)
    finite = _assert_same_nan_mask(result, reference)
    kwargs.setdefault("rtol", 0.0)
    kwargs.setdefault("atol", 0.0)
    actual, desired = result[finite], reference[finite]
    np.testing.assert_allclose(
        actual,
        desired,
        verbose=False,
        err_msg=(
            f"{int((actual != desired).sum())} of {actual.size} finite values differ from the reference. "
            "A nearest-neighbour difference is accepted only as a tie: see regrid_common.assert_nearest_ties"
        ),
        **kwargs,
    )


def _compare_linear(result, reference, skip_pole_rows=False, **kwargs):
    """Compare a linear result with its reference (``assert_allclose``).

    The whole array is compared -- the ``assert_allclose`` (``equal_nan=True`` by default)
    also requires the NaN positions of the two arrays to match. ``skip_pole_rows=True`` leaves
    rows 0 and -1 of a 2D result out of the comparison: the two rows of a duplicate pole are
    machine dependent in MIR, but only a grid that HAS pole rows may use the exemption, so the
    caller sets it explicitly (the matrix test: the global ``ll-10x10`` target of a pair that is
    not clipped) and it is never derived from the shape of the array.
    """
    result, reference = _prepare(result, reference)
    result, reference = _drop_pole_rows(result, reference, skip_pole_rows)
    np.testing.assert_allclose(result.flatten(), reference.flatten(), verbose=False, **kwargs)


def _compare_grid_box_average(result, reference, skip_pole_rows=False, **kwargs):
    """Compare a grid-box-average result with its reference (plain ``assert_allclose``).

    The whole array is compared; the ``assert_allclose`` (``equal_nan=True`` by default) also
    requires the NaN positions of the two arrays to match. ``skip_pole_rows=True`` is the
    explicit pole-row exemption of a global regular lat-lon result, see :func:`_compare_linear`.
    """
    result, reference = _prepare(result, reference)
    result, reference = _drop_pole_rows(result, reference, skip_pole_rows)
    np.testing.assert_allclose(result.flatten(), reference.flatten(), verbose=False, **kwargs)


def _drop_pole_rows(result, reference, skip_pole_rows):
    """Return ``result`` and ``reference`` without rows 0 and -1 when ``skip_pole_rows``.

    Without the flag the arrays come back unchanged. With it, the result must be a 2D grid
    with more than two rows, so that the exemption can never leave nothing to compare: a
    regional or clipped result of one or two rows would otherwise pass unchecked.
    """
    if not skip_pole_rows:
        return result, reference
    assert result.ndim == 2, f"skip_pole_rows needs a 2D result, this one has {result.ndim} dimension(s)"
    assert result.shape[0] > 2, (
        f"skip_pole_rows would drop rows 0 and -1 of a {result.shape[0]}-row result, "
        f"leaving {max(result.shape[0] - 2, 0)} rows to compare"
    )
    return result[1:-1], reference[1:-1]


INTERPOLATIONS: dict[str, Interpolation] = {
    interpolation.name: interpolation
    for interpolation in [
        Interpolation("linear", (), _compare_linear),
        Interpolation("nearest-neighbour", ("nn", "nearest-neighbor"), _compare_nearest),
        Interpolation("grid-box-average", (), _compare_grid_box_average),
    ]
}


def get_interpolation(name: str) -> Interpolation:
    """Return the registry entry for a canonical interpolation name or alias.

    Parameters
    ----------
    name : str
        A canonical name ("linear", "nearest-neighbour", "grid-box-average")
        or one of the aliases registered for it ("nn", "nearest-neighbor").

    Returns
    -------
    Interpolation
        The registry entry.

    Raises
    ------
    ValueError
        If ``name`` is neither a canonical name nor an alias.
    """
    interpolation = INTERPOLATIONS.get(name)
    if interpolation is not None:
        return interpolation
    for candidate in INTERPOLATIONS.values():
        if name in candidate.aliases:
            return candidate
    raise ValueError(f"Unknown interpolation {name!r}: known names are {sorted(INTERPOLATIONS)}")


def assert_nearest_ties(
    result,
    reference,
    input_values,
    source_spec,
    target_spec,
    *,
    max_neighbours: int = TIE_MAX_NEIGHBOURS,
    rel_tol: float = TIE_REL_TOL,
    abs_tol: float = TIE_ABS_TOL,
    context: str = "",
) -> None:
    """Assert that every difference between ``result`` and ``reference`` is a nearest-neighbour tie.

    The reference-free counterpart of :func:`_compare_nearest`: MIR's nearest neighbour is
    found by a search whose result depends on the build, so where several source points have
    the same distance to a target point another machine legitimately returns another one of
    their values. This function asks, for each position where ``result`` and ``reference``
    differ, whether the value in ``result`` is the input value of a source point that is as
    near to that target point as the nearest source point is -- within ``rel_tol`` relative and
    ``abs_tol`` absolute tolerance of its distance. Distances are the 3D chords of the unit
    sphere (``earthkit.geo.coord.latlon_to_xyz``) looked up in a ``scipy.spatial.cKDTree``,
    which is what a nearest-neighbour search compares anyway. A source or target given as a
    tuple of two arrays is taken as its ``(lats, lons)`` coordinates instead of a gridspec --
    anything else, as before, is a gridspec.

    It is meant to run *only* after the strict comparison has raised: the grids are built and
    their coordinates materialised here, so a pair whose arrays agree costs nothing (and, in
    fact, has nothing to explain -- calling it then fails). A differing NaN mask is never a tie
    and fails, and so does a value that is not the value of a near enough source point -- the
    independent check that "every non-NaN result value occurs in the input" cannot see the
    difference, a wrong source point's value *is* an input value.

    Parameters
    ----------
    result, reference : array-like
        The regridded array and the stored reference, same number of values.
    input_values : array-like
        The input field of the source grid, in the point order of ``source_spec`` (the array
        the regridding was called with).
    source_spec : dict or tuple of (array, array)
        Gridspec of the source the values came from, or its points as an ``(lats, lons)``
        pair of coordinate arrays, in the flat point order of ``input_values``. The
        coordinate spelling is for a source no gridspec describes, e.g. a GRIB field on a
        projection MIR does not return a spec for.
    target_spec : dict or tuple of (array, array)
        Gridspec the result is described by: the target spec, or the gridspec MIR returned for
        a clipped result -- or its points as an ``(lats, lons)`` pair, in the flat point
        order of ``result``. Its point order is the point order of ``result``; the result
        shape is only checked against the grid shape of a gridspec target.
    max_neighbours : int
        How many neighbours the KD-tree returns per differing point. Positions whose window is
        still full of equal distances (a whole source row at a pole) are re-queried without a
        limit, so this is a cost knob, not a correctness one.
    rel_tol, abs_tol : float
        Tie tolerance, in chord units: a source point is accepted as a tie candidate when its
        distance is at most ``(1 + rel_tol) * nearest_distance + abs_tol``.
    context : str
        Extra text at the head of the failure message: the matrix test names the
        (interpolation, source, target) pair there. The error of the strict
        comparison stays visible as the cause of the failure, since this runs
        inside its ``except`` block.

    Raises
    ------
    AssertionError
        For a differing NaN mask, an array that does not fit its grid, or a differing value
        that is not the input value of a source point near enough to be a tie candidate.
    """
    from scipy.spatial import cKDTree

    from earthkit.geo.grids import Grid

    result, reference = _prepare(result, reference)
    # a NaN where the other array has a value is not a tie: it is the wrong coverage
    finite = _assert_same_nan_mask(result, reference)
    differing = np.flatnonzero(finite & (result != reference))
    assert differing.size, (
        f"nothing to validate as a tie: the two arrays have the same NaN positions and {result.size} equal values"
    )

    values = np.asarray(input_values, dtype=np.float64).reshape(-1)
    # a coordinate pair has no grid metadata to check against: _xyz verifies its point count
    # against the values, and a shape it does not have cannot be compared. A gridspec is
    # built once here and the built grid is what _xyz reads, as before.
    if not _is_latlons(source_spec):
        source_spec = Grid(source_spec)
        assert values.size == source_spec.size(), (
            f"input has {values.size} values, the source grid has {source_spec.size()}"
        )
    if not _is_latlons(target_spec):
        target_spec = Grid(target_spec)
        assert result.shape == tuple(target_spec.shape), (
            f"result has shape {result.shape}, the target grid {target_spec.shape}"
        )

    # the point order of to_latlons() (or of a given coordinate pair) is the order of the
    # array's values: the input files are written as smooth_field(*grid.to_latlons()).reshape(grid.shape)
    source_xyz = _xyz(source_spec, values.size, "source")
    target_xyz = _xyz(target_spec, result.size, "target")

    tree = cKDTree(source_xyz)
    points = target_xyz[differing]
    k = max(1, min(int(max_neighbours), values.size))
    distances, indices = tree.query(points, k=k)
    distances = distances.reshape(len(points), -1)
    indices = indices.reshape(len(points), -1)

    limit = distances[:, :1] * (1.0 + rel_tol) + abs_tol
    candidates = distances <= limit
    got = result.reshape(-1)[differing]
    is_tie = np.zeros(len(points), dtype=bool)
    rows, cols = np.nonzero(candidates)
    is_tie[rows[values[indices[rows, cols]] == got[rows]]] = True

    # k neighbours may not have been enough to see the whole tie set: re-query those points
    # without a limit before calling a differing value wrong
    for row in np.flatnonzero(~is_tie & candidates[:, -1]):
        near = np.asarray(tree.query_ball_point(points[row], r=float(limit[row, 0])), dtype=np.intp)
        if near.size and np.any(values[near] == got[row]):
            is_tie[row] = True

    assert is_tie.all(), _tie_message(
        differing, got, reference.reshape(-1), result, distances, indices, values, is_tie, context, rel_tol, abs_tol
    )


def _is_latlons(x) -> bool:
    """True when an ``assert_nearest_ties`` grid argument is an ``(lats, lons)`` coordinate pair.

    A gridspec is a dict and a built ``Grid`` is not a 2-tuple, so a tuple of two
    coordinate arrays is the only spelling this accepts (see :func:`assert_nearest_ties`).
    """
    return isinstance(x, tuple) and len(x) == 2


def _xyz(grid, size, what):
    """The unit-sphere xyz of ``grid``'s points, in the order of its values' flat index.

    ``grid`` is a built ``Grid``, or an ``(lats, lons)`` pair of coordinate arrays, the
    coordinate spelling accepted by :func:`assert_nearest_ties`.
    """
    from earthkit.geo.coord import latlon_to_xyz

    if isinstance(grid, tuple):
        lats, lons = grid
    else:
        lats, lons = grid.to_latlons()
    lat = np.asarray(lats, dtype=np.float64).reshape(-1)
    lon = np.asarray(lons, dtype=np.float64).reshape(-1)
    assert lat.size == size and lon.size == size, f"{what} grid has {lat.size} coordinate pairs, expected {size}"
    return np.column_stack(latlon_to_xyz(lat, lon))


def _tie_message(differing, got, reference, result, distances, indices, values, is_tie, context, rel_tol, abs_tol):
    """The failure message of :func:`assert_nearest_ties`: up to five positions, in detail."""
    positions = np.unravel_index(differing, result.shape)
    lines = []
    for row in np.flatnonzero(~is_tie)[:5]:
        index = tuple(int(p[row]) for p in positions)
        near = indices[row][distances[row] <= distances[row, 0] * (1.0 + rel_tol) + abs_tol]
        lines.append(
            f"  position {index}: result {float(got[row])!r}, reference {float(reference[differing[row]])!r}; "
            f"the nearest source point is {distances[row, 0]:.12g} (chord) away and the tie candidates "
            f"there hold {[float(v) for v in sorted(set(values[near].tolist()))][:6]}"
        )
    failed = int((~is_tie).sum())
    more = f"\n  ... and {failed - len(lines)} more" if failed > len(lines) else ""
    prefix = f"{context}: " if context else ""
    return (
        f"{prefix}{failed} of {len(got)} values differ from the reference and are not the "
        f"value of a source point as near as the nearest one:\n" + "\n".join(lines) + more
    )


def _npz(name: str) -> str:
    """Return ``name`` with a guaranteed ``.npz`` suffix."""
    return name if name.endswith(".npz") else f"{name}.npz"


def _override_root() -> Path | None:
    """The local data tree named by ``EARTHKIT_GEO_REGRID_TEST_DATA``, or None when unset.

    Read at call time, not at import: the generator sets the variable for its
    own process and workers, and the tests set it with ``monkeypatch``.
    """
    value = os.environ.get(REGRID_TEST_DATA_ENV)
    return Path(value) if value else None


def _download(url: str, target: Path) -> None:
    """Fetch ``url`` and write its whole content to ``target`` in one go.

    The single network call of :func:`_fetch_new`, kept as its own function so
    that tests can fake or break it by monkeypatching this name (or
    ``requests.get``). The same request as
    ``earthkit.geo.utils.testing.simple_download`` -- redirects followed,
    :data:`DOWNLOAD_TIMEOUT` seconds per request -- except that the caller
    passes a temporary path which is renamed atomically afterwards, so a
    broken download cannot leave a partial file in the cache.
    """
    import requests

    response = requests.get(url, allow_redirects=True, timeout=DOWNLOAD_TIMEOUT)
    response.raise_for_status()
    Path(target).write_bytes(response.content)


def _fetch_new(relative: str | Path) -> Path:
    """Return the local path of a new-hierarchy file, downloading it when needed.

    This is the single place where the new data files are located (see the
    module docstring, "Where the new data files come from"):

    * with ``EARTHKIT_GEO_REGRID_TEST_DATA`` set, ``<that folder>/<relative>``
      is read directly and the network is never touched; a missing file is a
      ``FileNotFoundError`` naming :data:`GENERATOR_COMMAND`;
    * otherwise ``CACHE_ROOT/<relative>`` is used, fetched from
      ``URL_ROOT/<relative>`` when it is not cached yet: written to a temporary
      file in the same folder and renamed atomically, so a failed or
      interrupted download leaves no partial file behind.

    Parameters
    ----------
    relative : str or Path
        Path relative to the data root, e.g. "regrid/input/O32.npz".

    Returns
    -------
    Path
        The existing local file.

    Raises
    ------
    FileNotFoundError
        For a missing file in override mode (the message quotes
        :data:`GENERATOR_COMMAND`) and for a failed download (the message
        names the URL and the generator).
    """
    relative = Path(relative)
    override = _override_root()
    if override is not None:
        path = override / relative
        if not path.exists():
            raise FileNotFoundError(
                f"Regrid test data file {path} does not exist, although {REGRID_TEST_DATA_ENV} "
                f"points to the local tree {override}. Create it with the generator: {GENERATOR_COMMAND}"
            )
        return path

    target = CACHE_ROOT / relative
    if target.exists():
        return target

    url = join_url_path(URL_ROOT, *relative.parts)
    target.parent.mkdir(parents=True, exist_ok=True)
    # A hidden, non-.npz name in the same folder (same filesystem, so os.replace
    # is atomic) that a concurrent second process cannot clash with.
    temporary = target.parent / f".{target.name}.part.{os.getpid()}"
    try:
        _download(url, temporary)
    except Exception as error:
        temporary.unlink(missing_ok=True)
        raise FileNotFoundError(
            f"Regrid test data file {target} is not cached and downloading {url} failed "
            f"({type(error).__name__}: {error}). The data files are created by the generator "
            f"({GENERATOR_COMMAND}) and uploaded to the server by the maintainer; to work with "
            f"data that is not on the server, set {REGRID_TEST_DATA_ENV} to the folder that "
            f"contains the regrid/ tree (see tests/regrid/regrid_mir/README.md)."
        ) from error
    os.replace(temporary, target)
    return target


def data_path(
    kind: str,
    *parts: str,
    source: str | None = None,
    subfolder: str | None = None,
) -> Path:
    """Return the local path of a regrid test data file, fetching it if needed.

    Parameters
    ----------
    kind : {"input", "reference", "grib_reference", "legacy"}
        Which layout the parts follow:

        * "input": ``parts = (source_id,)`` ->
          ``<root>/regrid/input/<source_id>.npz``;
        * "reference": ``parts = (interpolation, source_id, target_id)`` ->
          ``<root>/regrid/<interpolation>/<source_id>/<target_id>.npz``,
          with the interpolation canonicalised through :func:`get_interpolation`
          (an alias names the same file as the canonical name);
        * "grib_reference": ``parts = (interpolation, filename, target_id)`` ->
          ``<root>/regrid/grib/<interpolation>/<filename>/<target_id>.npz``,
          the reference of the GRIB matrix (``regrid_grib_cases.py``), with the
          interpolation canonicalised as for "reference". The ``grib`` part is
          fixed and not an interpolation name, so this layout and the
          "reference" one can never name the same file;
        * "legacy": ``parts = (filename,)``, with ``subfolder`` (default
          "global_0_360"), delegated to
          ``earthkit.geo.utils.testing.get_test_data`` (downloads the file if
          it is not cached yet).

    For the new hierarchy ``<root>`` is the folder named by
    ``EARTHKIT_GEO_REGRID_TEST_DATA`` when that variable is set, and the cache
    :data:`CACHE_ROOT` otherwise, the missing file then downloaded from the
    server (see :func:`_fetch_new` and the module docstring, "Where the new
    data files come from").
    parts : str
        The path parts of the chosen layout, see ``kind``.
    source : {"new", "legacy"} or None
        Which data hierarchy to use; None (the default) derives it from
        ``kind`` ("legacy" for kind "legacy", "new" otherwise). A value that
        contradicts the kind raises ``ValueError``.
    subfolder : str or None
        Legacy subfolder, only for kind "legacy".

    Returns
    -------
    Path
        The local file path.

    Raises
    ------
    ValueError
        For an unknown kind, a wrong number of parts or a ``source``
        contradicting the kind.
    FileNotFoundError
        For a missing new-hierarchy file: in override mode (the message names
        :data:`GENERATOR_COMMAND`) or after a failed download from the server
        (the message names the URL).
    """
    if kind == "legacy":
        if source not in (None, "legacy"):
            raise ValueError(f"data_path('legacy', ...) needs source='legacy' or None, got {source!r}")
        if len(parts) != 1:
            raise ValueError(f"data_path('legacy', filename) takes exactly one part, got {len(parts)}")
        return Path(get_test_data(parts[0], subfolder=subfolder or DEFAULT_LEGACY_SUBFOLDER))

    if source not in (None, "new"):
        raise ValueError(f"data_path(kind={kind!r}) only supports source='new' or None, got {source!r}")

    if kind == "input":
        if len(parts) != 1:
            raise ValueError(f'data_path("input", source_id) takes exactly one part, got {len(parts)}')
        return _fetch_new(Path("regrid", "input", _npz(parts[0])))

    if kind == "reference":
        if len(parts) != 3:
            raise ValueError(
                'data_path("reference", interpolation, source_id, target_id) takes exactly three parts, '
                f"got {len(parts)}"
            )
        interpolation, source_id, target_id = parts
        interpolation = get_interpolation(interpolation).name
        return _fetch_new(Path("regrid", interpolation, source_id, _npz(target_id)))

    if kind == "grib_reference":
        if len(parts) != 3:
            raise ValueError(
                'data_path("grib_reference", interpolation, filename, target_id) takes exactly three parts, '
                f"got {len(parts)}"
            )
        interpolation, filename, target_id = parts
        interpolation = get_interpolation(interpolation).name
        return _fetch_new(Path("regrid", "grib", interpolation, filename, _npz(target_id)))

    raise ValueError(f"Unknown data_path kind {kind!r}: expected 'input', 'reference', 'grib_reference' or 'legacy'")


def load_values(path: str | Path) -> np.ndarray:
    """Read the ``arr_0`` array of an npz file.

    Parameters
    ----------
    path : str or Path
        The file to read, typically from :func:`data_path`.

    Returns
    -------
    numpy.ndarray
        The stored values.
    """
    return np.load(path)["arr_0"]


def save_values(path: str | Path, values: np.ndarray) -> Path:
    """Write ``values`` to the ``arr_0`` key of an npz file, making parent folders.

    Parameters
    ----------
    path : str or Path
        The file to write; it must end in ".npz" (as the names built by
        :func:`data_path` do), or :func:`numpy.savez_compressed` will append
        the extension itself.
    values : numpy.ndarray
        The array to store.

    Returns
    -------
    Path
        The written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, arr_0=values)
    return path


def legacy_values(filename: str, subfolder: str = DEFAULT_LEGACY_SUBFOLDER) -> np.ndarray:
    """Return the ``arr_0`` array of a legacy test data file, downloading it if needed.

    Parameters
    ----------
    filename : str
        File name inside the legacy ``subfolder``.
    subfolder : str
        Legacy subfolder, by default "global_0_360".

    Returns
    -------
    numpy.ndarray
        The stored values.
    """
    return load_values(data_path("legacy", filename, subfolder=subfolder))
