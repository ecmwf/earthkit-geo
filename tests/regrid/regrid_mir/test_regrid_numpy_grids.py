# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Matrix test: every covered source grid of the regrid table regridded onto every target.

Design
------
``regrid_cases.matrix(interpolation)`` yields the ``(source, target)`` pairs of one
interpolation; one parametrised test runs every pair of every interpolation in
``regrid_cases.SELECTION`` (nearest-neighbour: all covered sources against all four targets;
linear and grid-box-average: their small selections). Input values and MIR reference results
are npz files of the ``regrid/`` tree on the download server: ``regrid_common.data_path``
caches them under ``tests/data/regrid/`` and fetches what is missing, or -- with the
environment variable ``EARTHKIT_GEO_REGRID_TEST_DATA`` set -- reads that folder's tree
directly, which is how the files ``generate_regrid_data.py`` has just written (into
``_dev/test-data``) are tested before the upload. Data files are never written during a
pytest run. Per pair the test checks:

* the result shape equals the reference shape, and ``Grid(returned)`` of the returned
  gridspec has that shape as well;
* the returned gridspec: for an ordinary pair it describes the target
  (``Grid(returned) == Grid(target.spec)``, uid equality, so a spelling like ``{"grid":
  "H8"}`` matches the ring-order target); for a pair listed in ``regrid_cases.CLIPPED``
  (MIR clips the output to the source area, spells a too-small clip with a zero spacing or
  keeps the source rotation) it equals the gridspec recorded there, compared number-wise
  with a relative tolerance;
* the result against the saved reference through the registry compare of
  ``regrid_common.INTERPOLATIONS``: nearest-neighbour is strict (NaN at the same positions,
  every finite value equal) and, where it fails, falls back to
  ``regrid_common.assert_nearest_ties``, which accepts a differing value only if it is the
  input value of a source point as near to that target point as the nearest one -- the
  machine-dependent ties a nearest search can resolve differently elsewhere. A value that is
  not a tie candidate still fails. linear and grid-box-average compare the whole array with
  ``assert_allclose``, the pole rows left out only for a pair that really returns the whole
  global ``ll-10x10`` grid (target ``ll-10x10``, not clipped): ``skip_pole_rows`` says so, see
  ``regrid_common``;
* nearest only, independently of the reference: every non-NaN result value is exactly a
  value of the input.

Pairs that MIR cannot produce are listed in ``regrid_cases.BROKEN`` and become
``xfail(run=False)`` parameters with that reason, never silently skipped. The extra tests
check the table itself, that every non-covered source has a reason, that BROKEN and
CLIPPED name real matrix pairs, and how ``regrid_common.data_path`` finds its files: an
override tree holds exactly the files the table expects (the on-demand download cache
must hold none the table does not expect), a missing file in override mode fails with
the generator command in the message, and a cached file is fetched once with an atomic
rename while a failing download names the URL and leaves no partial file -- the three
data-access tests fake the download function with ``monkeypatch`` and never use the
network. The five named ``test_regrid_common_*`` check those comparison rules on
synthetic arrays, because the stored references cannot show what they must reject.

Extending
---------
Add a source to ``regrid_cases.SOURCES``, a target to ``regrid_cases.TARGETS``, or a new
interpolation (one ``INTERPOLATIONS`` registry entry in ``regrid_common.py`` plus one
``SELECTION`` entry -- both module docstrings explain how); this module needs no change, it
is parametrised from the tables -- a new compare rule just has to take the
``skip_pole_rows`` keyword the same way the existing ones do. Then generate the missing data
and run the tests against the local tree (the override; without it they would download the
old, upload-pending tree from the server):

    python tests/regrid/regrid_mir/generate_regrid_data.py --source <new-id>
    EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
        tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long

Finally the maintainer uploads the regenerated ``_dev/test-data/regrid/`` tree to the
download server, and plain runs fetch it into the cache again.

The module is marked ``long_test`` (skipped by the default ``-E short`` run) and needs
``mir``.
"""

import io
import math
import os
from pathlib import Path

import numpy as np
import pytest
import regrid_common
from regrid_cases import (
    BROKEN,
    CLIPPED,
    NOT_COVERED,
    SELECTION,
    SOURCES,
    check_table,
    matrix,
)
from regrid_common import (
    CACHE_ROOT,
    INTERPOLATIONS,
    REGRID_TEST_DATA_ENV,
    assert_nearest_ties,
    data_path,
    load_values,
    requires_mir,
)

from earthkit.geo.grids import Grid
from earthkit.geo.grids.array import regrid as regrid_array
from earthkit.geo.utils.testing import URL_ROOT
from earthkit.geo.utils.url import join_url_path

pytestmark = [pytest.mark.long_test, requires_mir]


def _matrix_params():
    """One ``pytest.param`` per (interpolation, source, target) pair of the whole matrix.

    The BROKEN pairs are marked ``xfail(run=False)`` with their table reason; the id is
    ``<interpolation>-<source id>-<target id>``.
    """
    params = []
    for interpolation in SELECTION:
        for source, target in matrix(interpolation):
            key = (interpolation, source.id, target.id)
            marks = []
            if key in BROKEN:
                marks.append(pytest.mark.xfail(run=False, reason=BROKEN[key]))
            id = f"{interpolation}-{source.id}-{target.id}"
            params.append(pytest.param(interpolation, source, target, id=id, marks=marks))
    return params


def _assert_gridspec(actual, expected, path="gridspec"):
    """Compare a returned gridspec with the expected one; numbers get a relative tolerance.

    Recurses through dicts (same key sets) and lists / tuples (same lengths); int and float
    values are compared with ``math.isclose`` (rel_tol 1e-9), because MIR recomputes some
    area bounds and may spell them with more digits than the table records. Anything else
    must compare equal.

    Parameters
    ----------
    actual : dict
        The gridspec returned by the regridding.
    expected : dict
        The gridspec the table says the result describes.
    path : str
        Name of ``actual`` in the failure messages, extended while recursing.
    """
    if isinstance(expected, dict):
        assert isinstance(actual, dict), f"{path}: {actual!r} is not a dict like {expected!r}"
        assert set(actual) == set(expected), f"{path}: keys {sorted(actual)} != {sorted(expected)}"
        for key, value in expected.items():
            _assert_gridspec(actual[key], value, f"{path}[{key!r}]")
    elif isinstance(expected, (list, tuple)):
        assert isinstance(actual, (list, tuple)), f"{path}: {actual!r} is not a list like {expected!r}"
        assert len(actual) == len(expected), f"{path}: {len(actual)} values, expected {len(expected)}"
        for index, (got, want) in enumerate(zip(actual, expected)):
            _assert_gridspec(got, want, f"{path}[{index}]")
    elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
        assert isinstance(actual, (int, float)) and not isinstance(actual, bool), f"{path}: {actual!r} is not a number"
        assert math.isclose(float(actual), float(expected), rel_tol=1e-9, abs_tol=0.0), (
            f"{path}: {actual!r} != {expected!r}"
        )
    else:
        assert actual == expected, f"{path}: {actual!r} != {expected!r}"


@pytest.mark.parametrize("interpolation,source,target", _matrix_params())
def test_regrid_numpy_grids(interpolation, source, target):
    """Regrid the stored input of ``source`` onto ``target`` and check shape, gridspec, values.

    Parameters
    ----------
    interpolation : str
        Canonical name of the interpolation, a key of ``regrid_cases.SELECTION``.
    source : regrid_cases.Source
        The source grid; its input values come from ``regrid/input/<source>.npz``.
    target : regrid_cases.Target
        The target grid; the reference comes from
        ``regrid/<interpolation>/<source>/<target>.npz``.
    """
    values = load_values(data_path("input", source.id))
    reference = load_values(data_path("reference", interpolation, source.id, target.id))

    result, returned = regrid_array(values, source.spec, target.spec, interpolation=interpolation)

    key = (interpolation, source.id, target.id)

    assert result.shape == reference.shape, f"{key}: result shape {result.shape}, reference {reference.shape}"

    if key in CLIPPED:
        _assert_gridspec(returned, CLIPPED[key], f"{key} gridspec")
    else:
        assert Grid(returned) == Grid(target.spec), (
            f"{key}: returned gridspec {returned!r} does not describe the target"
        )

    # the gridspec the result is described by must be the grid the values live on
    assert tuple(Grid(returned).shape) == result.shape, (
        f"{key}: Grid(returned) shape {Grid(returned).shape} != {result.shape}"
    )

    compare = INTERPOLATIONS[interpolation].compare
    if interpolation == "nearest-neighbour":
        # strict against the reference; where that fails the only accepted reason is a tie: a
        # source point as near to the target point as the nearest one, which another build of
        # MIR may legitimately have chosen instead. The validation builds both grids and their
        # coordinates, so it belongs on the failure path and nowhere else.
        try:
            compare(result, reference)
        except AssertionError:
            # the strict message stays visible as the cause of a failure of this call
            assert_nearest_ties(result, reference, values, source.spec, returned, context=str(key))

        # independent of the reference: a nearest result only moves values around, it never
        # invents any (NaN-filled positions outside the source area do not count)
        finite = result[~np.isnan(result)]
        assert np.all(np.isin(finite, values)), (
            f"{key}: result values that are not input values: {np.unique(finite[~np.isin(finite, values)])}"
        )
    else:
        # the pole rows of a global regular lat-lon result are machine dependent; the exemption
        # is asked for explicitly, and only where the pair really returns the whole target grid
        # (a clipped one may be a single row, and a single row has to be compared)
        compare(result, reference, skip_pole_rows=target.id == "ll-10x10" and key not in CLIPPED)


def test_regrid_numpy_grids_table():
    """Every source and target of the table builds and has the shape and size it claims."""
    # check_table() builds ~150 grids, which is why it runs here and not at import
    assert check_table() == []


def test_regrid_numpy_grids_coverage():
    """Every non-covered source has a reason; NOT_COVERED matches; BROKEN/CLIPPED are matrix pairs."""
    for source in SOURCES:
        if not source.covered:
            assert source.reason, f"{source.id}: not covered without a reason"
        else:
            assert not source.reason, f"{source.id}: covered but it has a reason {source.reason!r}"
    assert NOT_COVERED == {source.id: source.reason for source in SOURCES if not source.covered}

    pairs = {
        (interpolation, source.id, target.id) for interpolation in SELECTION for source, target in matrix(interpolation)
    }
    for key in BROKEN:
        assert key in pairs, f"BROKEN key {key} is not a pair of the matrix"
    for key in CLIPPED:
        assert key in pairs, f"CLIPPED key {key} is not a pair of the matrix"
    assert not set(BROKEN) & set(CLIPPED), f"BROKEN and CLIPPED share keys: {sorted(set(BROKEN) & set(CLIPPED))}"


def test_regrid_numpy_grids_broken_reasons():
    """Every BROKEN reason starts with "BROKEN: ", as the pytest xfail reason is built from it."""
    assert BROKEN, "BROKEN is empty, but the probe found failing pairs"
    for key, reason in BROKEN.items():
        assert reason.startswith("BROKEN: "), f"{key}: reason {reason!r}"


def _expected_npz_files():
    """The npz files the table expects: one input per covered source, one reference per
    matrix pair that is not BROKEN (the generator writes no file for a pair MIR cannot
    produce).

    The GRIB matrix keeps its references in the same root, under ``regrid/grib/``; they are
    entries of ``regrid_grib_cases.EXPECTED``, not of this table, and are not named here (the
    files check ignores that subtree, see below).
    """
    expected = {f"regrid/input/{source.id}.npz" for source in SOURCES if source.covered}
    expected |= {
        f"regrid/{interpolation}/{source.id}/{target.id}.npz"
        for interpolation in SELECTION
        for source, target in matrix(interpolation)
        if (interpolation, source.id, target.id) not in BROKEN
    }
    return expected


def test_regrid_numpy_grids_files():
    """The npz files on disk match the table: the exact set for a local tree, none extra for the cache.

    With ``EARTHKIT_GEO_REGRID_TEST_DATA`` set (the generator's freshly written tree,
    before it is uploaded) that folder must hold **exactly** the files the table expects:
    one input per covered source, one reference per matrix pair that is not BROKEN.
    Without the variable the files are downloaded from the server on demand, so the cache
    may hold any subset of them -- but no file the table does not expect.

    The ``regrid/grib/`` subtree of that root belongs to the GRIB matrix
    (``regrid_grib_cases.py``): it is not part of this table, so it is ignored here;
    ``test_regrid_fieldlist_grib_files`` is the files check of that subtree.
    """
    expected = _expected_npz_files()
    override = os.environ.get(REGRID_TEST_DATA_ENV)
    root = Path(override) if override else CACHE_ROOT
    folder = root / "regrid"
    found = list(folder.rglob("*.npz")) if folder.is_dir() else []
    actual = {path.relative_to(root).as_posix() for path in found if path.relative_to(folder).parts[0] != "grib"}

    if override:
        assert actual == expected, (
            f"{len(actual)} npz files under {folder} ({REGRID_TEST_DATA_ENV}), the table expects "
            f"{len(expected)}; missing: {sorted(expected - actual)[:5]}, extra: {sorted(actual - expected)[:5]}"
        )
    else:
        assert actual <= expected, (
            f"{len(actual - expected)} npz files under the download cache {folder} are not files of the "
            f"table (the cache is fetched on demand and may be deleted): {sorted(actual - expected)[:5]}"
        )


def test_regrid_numpy_grids_missing_override_file(tmp_path, monkeypatch):
    """Override mode: a missing file fails with the generator command in the message, no download."""
    monkeypatch.setenv(REGRID_TEST_DATA_ENV, str(tmp_path))

    def no_network(url, target):
        raise AssertionError(f"override mode must not download from {url}")

    monkeypatch.setattr(regrid_common, "_download", no_network)

    with pytest.raises(FileNotFoundError) as excinfo:
        data_path("input", "O32")
    message = str(excinfo.value)
    assert "generate_regrid_data.py" in message
    assert str(tmp_path / "regrid" / "input" / "O32.npz") in message


def test_regrid_numpy_grids_download_failure(tmp_path, monkeypatch):
    """Cache mode: a failing download names the URL in the error and leaves no partial file."""
    monkeypatch.delenv(REGRID_TEST_DATA_ENV, raising=False)
    monkeypatch.setattr(regrid_common, "CACHE_ROOT", tmp_path)

    def broken_download(url, target):
        Path(target).write_bytes(b"half of a file")  # like a connection that broke mid-transfer
        raise ConnectionError(f"cannot fetch {url}")

    monkeypatch.setattr(regrid_common, "_download", broken_download)

    with pytest.raises(FileNotFoundError) as excinfo:
        data_path("input", "O32")
    message = str(excinfo.value)
    assert join_url_path(URL_ROOT, "regrid", "input", "O32.npz") in message
    assert "generate_regrid_data.py" in message

    # the temp file is gone and the cache file was never created: the folder is empty again
    folder = tmp_path / "regrid" / "input"
    assert not list(folder.iterdir())


def test_regrid_numpy_grids_download_once(tmp_path, monkeypatch):
    """Cache mode: a file is fetched into the cache once and every later read uses the cache."""
    monkeypatch.delenv(REGRID_TEST_DATA_ENV, raising=False)
    monkeypatch.setattr(regrid_common, "CACHE_ROOT", tmp_path)

    calls = []

    def fake_download(url, target):
        # writes the npz bytes to exactly the path it is given, like _download does
        calls.append(url)
        buffer = io.BytesIO()
        np.savez_compressed(buffer, arr_0=np.array([1.0, 2.0, 3.0]))
        Path(target).write_bytes(buffer.getvalue())

    monkeypatch.setattr(regrid_common, "_download", fake_download)

    path = data_path("input", "O32")
    assert path == tmp_path / "regrid" / "input" / "O32.npz"
    assert load_values(path).tolist() == [1.0, 2.0, 3.0]
    assert data_path("input", "O32") == path

    assert calls == [join_url_path(URL_ROOT, "regrid", "input", "O32.npz")]
    assert not [p for p in path.parent.iterdir() if p != path]  # no temp file left behind


# ---------------------------------------------------------------------------
# the comparison rules of regrid_common, on synthetic arrays: the stored references cannot
# show what they must reject (a one-row clipped result, a wrong value in a pole row) and what
# they must accept (an equidistant source point), so those cases are built here.

# A source with two points at the same distance from a target point: (0, 5) is as far from
# (0, 0) as from (0, 10), so either of their values is a legitimate nearest result while the
# value of (10, 5) or (-10, 5) -- 10 degrees away -- is not.
_TIE_SOURCE = {
    "type": "unstructured_ll",
    "latitudes": [0.0, 0.0, 10.0, -10.0],
    "longitudes": [0.0, 10.0, 5.0, 5.0],
    "name": "EARTHKIT-GEO-T8-TIE-IN",
}
_TIE_TARGET = {
    "type": "unstructured_ll",
    "latitudes": [0.0, 0.0],
    "longitudes": [5.0, 0.0],
    "name": "EARTHKIT-GEO-T8-TIE-OUT",
}
_TIE_INPUT = np.array([1.0, 2.0, 3.0, 4.0])  # the value of each source point
_TIE_REFERENCE = np.array([1.0, 1.0])  # both target points took source point 0


def test_regrid_common_compare_nearest_is_strict():
    """The nearest rule needs the NaN positions to match and every finite value to be equal."""
    nearest = INTERPOLATIONS["nearest-neighbour"].compare
    values = np.linspace(1.0, 10.0, 20).reshape(4, 5)

    def changed(row, column, new):
        other = values.copy()
        other[row, column] = new
        return other

    assert nearest(values, values.copy()) is None
    assert nearest(np.array([]), np.array([])) is None

    # one value, one unit in the last place, is enough to fail
    with pytest.raises(AssertionError):
        nearest(changed(1, 2, np.nextafter(values[1, 2], 0.0)), values)

    # the pole rows are compared like any other row: no exemption, not even a tolerant one
    for row in (0, -1):
        with pytest.raises(AssertionError):
            nearest(changed(row, 0, values[row, 0] * 1.5), values)

    # and so is the single row of a clipped one-row result, which the old shape-based
    # exemption left completely unchecked
    with pytest.raises(AssertionError):
        nearest(np.array([[1.0, 2.0, 9.0]]), np.array([[1.0, 2.0, 3.0]]))

    # NaN: the same positions pass, a flipped value does not
    with_nan = values.copy()
    with_nan[0, 0] = np.nan
    assert nearest(with_nan, with_nan.copy()) is None
    with pytest.raises(AssertionError):
        nearest(with_nan, values)

    # the tolerance is a kwarg, its default is exact
    assert nearest(changed(2, 2, values[2, 2] * 1.0000001), values, rtol=1e-6) is None
    # a nearest rule never skips rows: the flag is refused rather than silently ignored
    with pytest.raises(AssertionError):
        nearest(values, values, skip_pole_rows=True)


@pytest.mark.parametrize("interpolation", ["linear", "grid-box-average"])
def test_regrid_common_compare_pole_rows_are_opt_in(interpolation):
    """The allclose rules compare the whole array; only the caller may ask for the pole rows."""
    compare = INTERPOLATIONS[interpolation].compare
    values = np.arange(19 * 36, dtype=np.float64).reshape(19, 36)

    def changed(rows):
        other = values.copy()
        for row in rows:
            other[row] += 1e6
        return other

    assert compare(values, values.copy()) is None
    # a difference anywhere, including the pole rows, fails by default ...
    for rows in ([9], [0], [-1], [0, -1]):
        with pytest.raises(AssertionError):
            compare(changed(rows), values)
    # ... and only the two pole rows are exempted when the flag is set
    for rows in ([0], [-1], [0, -1]):
        assert compare(changed(rows), values, skip_pole_rows=True) is None
    for rows in ([9], [1], [-2]):
        with pytest.raises(AssertionError):
            compare(changed(rows), values, skip_pole_rows=True)

    # the exemption may not leave nothing behind: one-row and two-row results, and a 1D
    # result, are refused rather than compared against nothing
    for array in (np.array([[1.0, 2.0, 3.0]]), np.zeros((2, 3)), values.ravel()):
        with pytest.raises(AssertionError):
            compare(array, array.copy(), skip_pole_rows=True)


def test_regrid_common_nearest_ties():
    """A differing value is a tie exactly when a source point that near holds it."""
    tie = assert_nearest_ties

    # the other point of the pair that is equally near: accepted, and the strict rule did fail
    # on it, which is what the fallback is for
    with pytest.raises(AssertionError):
        INTERPOLATIONS["nearest-neighbour"].compare(np.array([2.0, 1.0]), _TIE_REFERENCE)
    tie(np.array([2.0, 1.0]), _TIE_REFERENCE, _TIE_INPUT, _TIE_SOURCE, _TIE_TARGET)

    # a value of a source point 10 degrees away is not a tie -- that is the difference the
    # "every result value occurs in the input" check cannot see; a value the input does not
    # hold at all is wrong in both rules
    for wrong in ([3.0, 1.0], [4.0, 1.0], [1.5, 1.0]):
        with pytest.raises(AssertionError):
            tie(np.array(wrong), _TIE_REFERENCE, _TIE_INPUT, _TIE_SOURCE, _TIE_TARGET)

    # a NaN where the reference has a value is wrong coverage, not a tie
    with pytest.raises(AssertionError):
        tie(np.array([np.nan, 1.0]), _TIE_REFERENCE, _TIE_INPUT, _TIE_SOURCE, _TIE_TARGET)

    # nothing differing means the strict comparison had no value to explain
    with pytest.raises(AssertionError):
        tie(_TIE_REFERENCE, _TIE_REFERENCE, _TIE_INPUT, _TIE_SOURCE, _TIE_TARGET)

    # an array that does not fit the grid it claims is refused, whatever its values are
    with pytest.raises(AssertionError):
        tie(np.array([2.0, 1.0, 3.0]), np.array([1.0, 1.0, 3.0]), _TIE_INPUT, _TIE_SOURCE, _TIE_TARGET)


def test_regrid_common_nearest_ties_wider_than_the_query_window():
    """A tie set bigger than the KD-tree window -- a whole pole row of one point -- is still found."""
    source, target = {"grid": [5, 5]}, {"grid": [10, 10]}
    values = np.arange(37 * 72, dtype=np.float64).reshape(37, 72)
    reference = np.zeros((19, 36))
    reference[0, 0] = values[0, 0]  # the pole point the reference picked

    # every point of row 0 is at the north pole: any of them is a tie candidate at (0, 0),
    # and there are 72 of them, more than the 16 the tree is asked for
    other = reference.copy()
    other[0, 0] = values[0, 71]
    assert_nearest_ties(other, reference, values, source, target)

    # the point 5 degrees down the meridian is not one of them
    wrong = reference.copy()
    wrong[0, 0] = values[1, 0]
    with pytest.raises(AssertionError):
        assert_nearest_ties(wrong, reference, values, source, target)
