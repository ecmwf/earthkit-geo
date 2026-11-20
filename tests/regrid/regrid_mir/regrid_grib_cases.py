# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""The GRIB regrid matrix table: the 40 matrix GRIB files onto every target, plus the negative set.

This module is data, not a test: it is imported by the GRIB fieldlist matrix test and by
``generate_regrid_grib_data.py``, and it has no side effects at import -- no earthkit-data,
no ``mir``, no download, no grid built. It imports ``pytest`` (for :func:`cases` and
:func:`invalid_cases`) and the standard library only; exception classes named in
:class:`Broken` entries are imported lazily by :func:`cases`, and those named in
:class:`Invalid` entries are resolved by the negative test itself -- so merely reading the
tables never touches the GRIB stack.

The matrix: :data:`GRIB_FILES` (40 of the 42 GRIB files of the download server's
``test-data/grib/`` folder -- the two whose grid metadata contradicts itself are not among
them, see :data:`INVALID_GRIB_FILES` below), each read with
``earthkit.data.from_source("file", ...)`` into a fieldlist and regridded by
``earthkit.geo.regrid`` with each :data:`GRIB_COMBINATIONS` (interpolation, target) pair
onto one of :data:`GRIB_TARGETS`. MIR's result of one triple is its reference: the values of
all result fields, stacked into one npz at

    <root>/regrid/grib/<interpolation>/<filename>/<target id>.npz      (arr_0, shape (n_fields, *field shape))

read by ``regrid_common.data_path("grib_reference", interpolation, filename, target_id)``
and written by ``generate_regrid_grib_data.py``.

Contents
--------
* :data:`GRIB_FILES`: the 40 GRIB file names of the server's ``grib/`` folder that the matrix
  regrids.
* :class:`GribTarget` / :data:`GRIB_TARGETS`: the grids the files are regridded to.
* :data:`GRIB_COMBINATIONS`: the (interpolation, target id) pairs that run.
* :class:`Expected`: what MIR produced for one triple (per result field: shape, gridType,
  returned gridspec, NaN count, plus the number of fields).
* :class:`Broken`: a triple MIR cannot produce (it raises, aborts or hangs), turned into an
  xfail mark by :func:`cases`.
* :class:`KnownBug`: a triple MIR produces a result for, but a verified MIR bug makes that
  result wrong (regrid-grib-findings §C). It records the CORRECT output and becomes a strict
  xfail through :func:`cases`: the case runs, fails against the correct result while MIR is
  buggy, and an XPASS failure is the alert that MIR got fixed.
* :data:`EXPECTED`: triple -> Expected | Broken | KnownBug, the observations of the generator
  run plus the correct expectations of the known MIR bugs.
* :func:`cases`: one ``pytest.param`` per triple, with the xfail marks of the BROKEN and the
  known-bug ones.
* :data:`INVALID_GRIB_FILES` / :class:`Invalid` / :data:`INVALID` / :func:`invalid_cases`:
  the negative set -- the 2 server files that are not in the matrix because their own grid
  metadata contradicts itself (regrid-grib-findings §A), one :class:`Invalid` entry per
  combination that names the exception and the message their regridding raises, and the
  parameters of the test that asserts them with ``pytest.raises``.

Adding a file
-------------
Add its name to :data:`GRIB_FILES`; it must exist under ``URL_ROOT/grib/`` on the download
server (that is where ``data_path("legacy", filename, subfolder="grib")`` fetches it from).
Then generate its references, paste the printed entries into :data:`EXPECTED`, run the tests
against the local tree and upload once they pass (see "Generating and uploading" below):

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file <new-name.grib>

A file that cannot be regridded at all -- one whose own grid metadata contradicts itself --
is not matrix material: add it to :data:`INVALID_GRIB_FILES` with one :class:`Invalid` entry
per :data:`GRIB_COMBINATIONS` pair instead. It then has no reference and the generator does
not probe it; the negative test asserts the error its regridding raises.

Adding a target or a combination
--------------------------------
Append ``GribTarget(id, spec)`` to :data:`GRIB_TARGETS`, then add the ``(interpolation, id)``
pair to :data:`GRIB_COMBINATIONS` for every interpolation that should run onto it. The
interpolation must be registered in ``regrid_common.INTERPOLATIONS`` (canonical name; the
folder of the reference files is the canonicalised name). New pairs have no references and
no :data:`EXPECTED` entries until the generator has been run for them.

Generating and uploading
------------------------

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py                 # everything
    python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file O32.grib2 --interpolation nn

The script runs each triple in its own subprocess, writes the references under its output
root (``_dev/test-data`` by default, ``--output-dir`` to choose), prints a summary and the
pasteable :data:`EXPECTED` entries, and saves the observed expectations as JSON
(``<root>/regrid-grib-expected.json``, ``--report`` to choose). A triple that raises becomes
a ``Broken(kind="raises", exception=..., reason="BROKEN: ...")`` entry with the exception
type it observed (module-qualified when not a builtin); a timeout or abort becomes
``Broken(kind="hang", ...)``. Paste those, and the passing triples' ``Expected`` entries,
into :data:`EXPECTED` and commit -- the table is the record, the generator only observes it.
The :class:`KnownBug` entries are the exception: the generator observes the buggy MIR result of
those triples again and prints ``Expected`` lines for them, which must NOT be pasted over a
KnownBug entry -- replace such an entry (and regenerate its reference) only once MIR is fixed.
Then test the tree through the override before the upload (the tests would otherwise
download the old tree):

    EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
        tests/regrid/regrid_mir/test_regrid_fieldlist_grib.py -E long --no-cov -k grib

and upload the ``regrid/grib/`` subtree to the server -- the maintainer's step; this folder
never uploads and no test ever writes a data file.

Note on the shared ``regrid/`` root: the ``grib/`` subtree sits in the same tree as the
numpy-matrix data (``grib`` is not an interpolation name, so no ``data_path("reference",
...)`` path can ever collide with a ``grib_reference`` one), but it is NOT part of that
table: ``test_regrid_numpy_grids.py`` checks the npz files under ``<root>/regrid/`` against
``regrid_cases.py`` and ignores this ``grib/`` subtree, which ``test_regrid_fieldlist_grib.py``
checks against the table above.
"""

import builtins
import importlib
from dataclasses import dataclass

import pytest

# ---------------------------------------------------------------------------
# the files: the GRIB files of the download server's test-data/grib/ folder
# ---------------------------------------------------------------------------

# Names as they are on the server (data_path("legacy", <name>, subfolder="grib")). The
# matrix reads each with earthkit-data from_source("file", ...) and regrids the whole
# fieldlist, so files with several fields or several grids are covered as they are.
# 40 of the 42 files of that folder: the two whose grid metadata contradicts itself are in
# INVALID_GRIB_FILES below, tested as negative cases instead of regridded.
GRIB_FILES: tuple[str, ...] = (
    "5x5_multi.grib1",
    "5x5_multi.grib2",
    "alternate-scanning.grib",
    "different_grid_resolutions.grib",
    "eORCA025_T.grib",
    "healpix_H8.grib2",
    "healpix_H8_nested.grib2",
    "icon_ch1.grib2",
    "lambert_conformal.grib",
    "lambert_conformal_na.grib2",
    "mercator.grib",
    "mercator_pacific.grib2",
    "polar_stereographic_alaska.grib2",
    "polar_stereographic_north.grib1",
    "reduced_gg_N32.grib1",
    "reduced_gg_N32.grib2",
    "reduced_gg_O32.grib1",
    "reduced_gg_O32.grib2",
    "reduced_rotated_gg_N32.grib1",
    "reduced_rotated_gg_N48_global.grib1",
    "reduced_rotated_gg_O32.grib1",
    "reduced_rotated_gg_subarea_N32.grib1",
    "reduced_rotated_gg_subarea_O32.grib1",
    "regular_gg_F16.grib1",
    "regular_ll.grib1",
    "regular_ll.grib2",
    "regular_ll_single_point.grib1",
    "regular_ll_single_point_other.grib",
    "regular_ll_subarea.grib1",
    "regular_rotated_gg_F48.grib1",
    "rotated_ll.grib1",
    "rotated_ll_-180_180_5x5.grib1",
    "rotated_ll_-180_180_5x5.grib2",
    "rotated_ll_0_360_5x5.grib1",
    "rotated_ll_0_360_5x5.grib2",
    "rotated_ll_subarea.grib1",
    "rotated_ll_subarea.grib2",
    "scanning_mode_64.grib",
    "sh_t32.grib1",
    "sh_t63.grib",
)

# The two server files that are not in the matrix: each contradicts itself, so no regridding
# of them can be right, and regrid fails on the contradiction (regrid-grib-findings §A). They
# are the negative set: INVALID below records the error each combination raises, and
# invalid_cases() parametrises the test that asserts it. No reference exists for them, and
# generate_regrid_grib_data.py does not probe them: it iterates GRIB_FILES only.
INVALID_GRIB_FILES: tuple[str, ...] = (
    "regular_gg_wrong_increment.grib",
    "regular_ll_wrong_increment.grib",
)

# ---------------------------------------------------------------------------
# the targets and the (interpolation, target) combinations that run
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GribTarget:
    """One gridspec the GRIB fieldlists are regridded to.

    Parameters
    ----------
    id: str
        Short filesystem-safe name, unique among the targets and used in file names.
    spec: dict
        The gridspec passed to ``earthkit.geo.regrid`` as ``out_grid``.
    """

    id: str
    spec: dict


GRIB_TARGETS: tuple[GribTarget, ...] = (
    GribTarget("10x10", {"grid": [10, 10]}),
    GribTarget("O32", {"grid": "O32"}),
)

# The (canonical interpolation, target id) pairs of the matrix, in the order the generator
# runs them. The interpolations not named here are registered in regrid_common but not part
# of this matrix.
GRIB_COMBINATIONS: tuple[tuple[str, str], ...] = (
    ("nearest-neighbour", "10x10"),
    ("linear", "10x10"),
    ("nearest-neighbour", "O32"),
)

# ---------------------------------------------------------------------------
# what one triple does
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Expected:
    """What MIR produced for one (interpolation, file, target) triple.

    Every attribute is per result field: entry ``i`` describes field ``i`` of the regridded
    fieldlist, in order, as ``generate_regrid_grib_data.py`` observed it. The reference npz
    holds the field values stacked into one array of shape ``(n_fields, *field shape)``.

    Parameters
    ----------
    n_fields: int
        Number of fields of the regridded fieldlist; length of every per-field tuple below.
    shapes: tuple of tuple of int
        ``numpy.ndarray.shape`` of each field's values.
    grid_types: tuple
        ``metadata("gridType")`` of each field, None when the field has no gridType.
    grid_specs: tuple
        ``geography.grid_spec()`` of each field, None when MIR returns no spec for it
        (projection grids do that).
    nan_counts: tuple of int
        Number of NaN values of each field (a regional source NaN-fills the target points
        outside its area).
    """

    n_fields: int
    shapes: tuple = ()
    grid_types: tuple = ()
    grid_specs: tuple = ()
    nan_counts: tuple = ()

    def __post_init__(self):
        if self.n_fields < 1:
            raise ValueError(f"Expected needs n_fields >= 1, got {self.n_fields}")
        for name in ("shapes", "grid_types", "grid_specs", "nan_counts"):
            entries = getattr(self, name)
            if len(entries) != self.n_fields:
                raise ValueError(f"Expected.{name} has {len(entries)} entries, n_fields is {self.n_fields}")


@dataclass(frozen=True)
class Broken:
    """A triple MIR cannot produce: it raises, aborts the process or hangs.

    Parameters
    ----------
    kind: str
        "raises": the regridding fails with ``exception`` and pytest may run the case
        (``xfail(raises=..., strict=True)``); "hang": it times out or aborts the
        interpreter, so the case is never executed (``xfail(run=False)``).
    exception: str
        For kind "raises": the exception type name, module-qualified when it is not a
        builtin (e.g. "gribapi.errors.WrongGridError"). :func:`cases` imports it lazily.
    reason: str
        The xfail reason, starting with "BROKEN: " and quoting what happens.
    """

    kind: str
    exception: str = ""
    reason: str = ""

    def __post_init__(self):
        if self.kind not in ("raises", "hang"):
            raise ValueError(f'Broken.kind must be "raises" or "hang", got {self.kind!r}')
        if self.kind == "raises" and not self.exception:
            raise ValueError(f"Broken(kind='raises') needs an exception name, reason {self.reason!r}")
        if not self.reason.startswith("BROKEN: "):
            raise ValueError(f"Broken.reason must start with 'BROKEN: ', got {self.reason!r}")


@dataclass(frozen=True)
class KnownBug:
    """A triple MIR produces a result for, but the result is wrong: a verified MIR bug.

    The mirror of :class:`Broken`: MIR can produce something, it is simply not what it should
    be. Every verified bug in this table is the same one (the regrid-grib findings,
    §C): a rotated input regridded onto the bare target of its own geometry (``rotated_ll``
    10x10 -> ``{"grid": [10, 10]}``, rotated O32 -> ``{"grid": "O32"}``) comes back
    byte-for-byte unchanged, rotation kept, while a target of other geometry interpolates and
    drops the rotation. The target is unrotated, so the correct output is unrotated.

    The entry records that CORRECT output, with the same per-field meaning as in
    :class:`Expected`; a field left ``None`` is unknown -- the post-fix clip of a sub-area
    source cannot be observed while the bug hides it -- and the matrix test then asserts only
    that the result is unrotated (its gridspec holds no ``projection`` key). :func:`cases`
    marks the triple ``xfail(strict=True)``: it runs and fails against the correct result
    while MIR is buggy, and the day MIR is fixed it XPASSes -- a failure that says so, the
    alert to regenerate the triple's reference (its stored npz is the stale no-op output) and
    turn the entry into :class:`Expected`.

    Parameters
    ----------
    reason: str
        The xfail reason, starting with "KNOWN MIR BUG: " and citing regrid-grib-findings §C.
    n_fields: int
        Number of fields of the regridded fieldlist; length of every per-field tuple given.
    shapes: tuple of tuple of int, or None
        ``numpy.ndarray.shape`` of each field's values, None to not check the shapes.
    grid_types: tuple, or None
        ``metadata("gridType")`` of each field, None to not check it.
    grid_specs: tuple, or None
        ``geography.grid_spec()`` of each field, None to only check that the spec of a result
        carries no rotation projection.
    nan_counts: tuple of int, or None
        Number of NaN values of each field, None to not check it.
    """

    reason: str
    n_fields: int
    shapes: tuple | None = None
    grid_types: tuple | None = None
    grid_specs: tuple | None = None
    nan_counts: tuple | None = None

    def __post_init__(self):
        if not self.reason.startswith("KNOWN MIR BUG: "):
            raise ValueError(f"KnownBug.reason must start with 'KNOWN MIR BUG: ', got {self.reason!r}")
        if self.n_fields < 1:
            raise ValueError(f"KnownBug needs n_fields >= 1, got {self.n_fields}")
        for name in ("shapes", "grid_types", "grid_specs", "nan_counts"):
            entries = getattr(self, name)
            if entries is not None and len(entries) != self.n_fields:
                raise ValueError(f"KnownBug.{name} has {len(entries)} entries, n_fields is {self.n_fields}")


@dataclass(frozen=True)
class Invalid:
    """What one regridding of an :data:`INVALID_GRIB_FILES` file must do: fail with a given error.

    The mirror of the matrix entries: instead of a result, the file's own grid metadata
    contradicts itself, so the regridding cannot be right and must be refused (the
    regrid-grib findings, §A). Nothing is generated for such a triple and no reference exists:
    the :func:`invalid_cases` parameters go to a test that runs ``regrid`` and asserts this
    error with ``pytest.raises``.

    Parameters
    ----------
    exception: str
        The exception type the regridding raises, module-qualified when it is not a builtin
        (e.g. "gribapi.errors.WrongGridError"). The negative test imports it lazily, through
        :func:`_exception_class`, so reading this table never touches the GRIB stack.
    match: str
        Regular expression the exception message must contain, as ``pytest.raises(match=...)``
        takes it: the part of the message that says what is wrong with the file.
    reason: str
        What the file contradicts and who is at fault, starting with "INVALID INPUT: " and
        citing regrid-grib-findings §A.
    """

    exception: str
    match: str
    reason: str

    def __post_init__(self):
        if not self.exception:
            raise ValueError(f"Invalid needs an exception name, reason {self.reason!r}")
        if not self.match:
            raise ValueError(f"Invalid needs a message match, reason {self.reason!r}")
        if not self.reason.startswith("INVALID INPUT: "):
            raise ValueError(f"Invalid.reason must start with 'INVALID INPUT: ', got {self.reason!r}")


# ---------------------------------------------------------------------------
# the table: triple -> what it does
# ---------------------------------------------------------------------------

# (interpolation, filename, target id) -> Expected | Broken | KnownBug, filled from generator runs.
# The BROKEN entries are the triples MIR cannot produce: the reason quotes the exception the
# probe observed with its message and names the cause in the regrid-grib findings, whose
# section letters the comments use. The KnownBug entries are the triples MIR produces wrongly
# -- they record the CORRECT result, not what MIR answers, and the matrix test skips their
# stale reference npz (see the class docstring).
# The Expected entries below are what ``generate_regrid_grib_data.py`` observed for the rest of
# the matrix on 2026-11-20, the 111-triple table of that day (95 triples ran, 94 of them are
# entries here: the mercator linear one ran too slowly to test, see BROKEN below), plus the 6
# spectral triples the same script observed after the spectral-input fix on the same day, plus
# the 15 triples of the five files the server's grib/ folder gained that evening, each file run
# on its own with --timeout 600. 115 generator entries in all; 6 of
# them -- the rotated inputs MIR copies back unchanged onto the bare target of their own
# geometry -- were verified as a MIR bug on 2026-11-20 and became KnownBug entries (§C): 109
# Expected and 6 KnownBug. A triple that looks odd but is not a verified bug is recorded as
# observed, with a comment that says where the oddity is explained.
# That is the whole table: 120 triples, the 40 files of GRIB_FILES by the 3 combinations. The
# two wrong-increment files were the 6 BROKEN "raises" entries of §A until 2026-11-20, when
# they became the negative set instead -- the matrix went to 35 files and 105 triples, the
# BROKEN entries to the 5 that are left here, and their 6 observations moved to INVALID below;
# the five files added that evening brought it to 40 files and 120 triples. The file check
# counts exactly one npz per non-BROKEN triple: 115 references, and a §A or BROKEN triple
# never has one.
EXPECTED: dict[tuple[str, str, str], Expected | Broken | KnownBug] = {
    # --- FIXED: A. the spectral files, sh_t32.grib1 and sh_t63.grib, all three combinations of both.
    # They were BROKEN (AttributeError: 'NoneType' object has no attribute 'shape'): earthkit-geo asked
    # such a field for its latitudes/longitudes before MIR saw the message, and a spectral truncation
    # has none. Fixed in earthkit-geo on 2026-11-20: regrid_grib_field sends a field whose GRIB
    # gridType is "sh" straight to the backend's regrid_grib, which regrids spectral GRIB natively
    # (spectral -> grid through its Legendre machinery). The six Expected entries
    # they became sit at the end of the combination sections below, in table order.
    # --- BROKEN: C. the one-point grids (§B, responsible: MIR). Each field of
    # regular_ll_single_point_other.grib is a single node, and no interpolation stencil can be
    # built from one point: MIR fails inside itself instead of reporting a degenerate grid. The
    # other single-point file reaches the same assert on the way to O32; regridding it to its own
    # 10x10 frame does not interpolate at all and works, see the Expected entries.
    ("nearest-neighbour", "regular_ll_single_point.grib1", "O32"): Broken(
        kind="raises",
        exception="RuntimeError",
        reason="BROKEN: RuntimeError: Assertion failed: !first -- MIR cannot build a reduced Gaussian ring from"
        " a one-point input grid; §B, responsible: MIR",
    ),
    ("nearest-neighbour", "regular_ll_single_point_other.grib", "10x10"): Broken(
        kind="raises",
        exception="RuntimeError",
        reason="BROKEN: RuntimeError: Serious bug: KNearestNeighbours: failed to interpolate -- no nearest"
        " neighbour can be searched on a one-point input grid; §B, responsible: MIR",
    ),
    ("linear", "regular_ll_single_point_other.grib", "10x10"): Broken(
        kind="raises",
        exception="RuntimeError",
        reason="BROKEN: RuntimeError: Assertion failed: FiniteElement: non-empty interpolation (to at least one"
        " point) -- no element can be built on a one-point input grid; §B, responsible: MIR",
    ),
    ("nearest-neighbour", "regular_ll_single_point_other.grib", "O32"): Broken(
        kind="raises",
        exception="RuntimeError",
        reason="BROKEN: RuntimeError: Assertion failed: !first -- MIR cannot build a reduced Gaussian ring from"
        " a one-point input grid; §B, responsible: MIR",
    ),
    # --- BROKEN: F. mercator onto the 10x10 lat-lon grid: not a hang, it is too slow to be a
    # test. The 4.5 M scattered points of the array path (earthkit-geo keeps projected grids
    # there, §E) make MIR build its interpolation weights from scratch, which measured 233.3 s
    # on 2026-11-20, after the earthkit-data fixes; with MIR's weight cache in /tmp/cache/mir
    # the same case takes about 2 s (2.4 s re-measured the same day). A long_test case that
    # starts with a four minute cold run, and aborts the whole pytest process if it ever
    # exceeds its 300 s thread timeout, is not worth its one reference, so it never runs. The
    # probe that timed that cold run did produce a result for it (19x36, regular_ll, spec
    # {"grid": [10, 10]}, 595 NaN: the input's 9999 fill now arrives as NaN), but no reference
    # is kept for a triple that is never tested: the file check counts exactly one npz per
    # triple that runs, and this one does not run.
    ("linear", "mercator.grib", "10x10"): Broken(
        kind="hang",
        reason="BROKEN: slow: linear onto 10x10 needs about 230 s the first time, while MIR builds the"
        " interpolation weights of 4.5 M scattered points (2 s once they are cached); never run. §E",
    ),
    # ---------------------------------------------------------------------------
    # The triples MIR produced, as generate_regrid_grib_data.py observed them on 2026-11-20
    # (--timeout 120, _dev/test-data/regrid-grib-expected.json), in GRIB_COMBINATIONS order and,
    # within a combination, in GRIB_FILES order.
    # --- nearest-neighbour onto the 10x10 global lat-lon grid ---
    ("nearest-neighbour", "5x5_multi.grib1", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    ("nearest-neighbour", "5x5_multi.grib2", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §D: 2x3 clip: only 6 of the 684 target nodes are inside the source bbox
    ("nearest-neighbour", "alternate-scanning.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((2, 3),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [50, -10, 40, 10], "grid": [10, 10]},),
        nan_counts=(0,),
    ),
    ("nearest-neighbour", "different_grid_resolutions.grib", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §H: only 36 of the 684 values are NaN -- the latitude -90 row -- although the
    # input's bitmap masks 44,442 of its 1.74 M points: MIR resamples the unstructured input onto O400 first,
    # and only the south pole is forced onto a NaN node of that intermediate grid
    ("nearest-neighbour", "eORCA025_T.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(36,),
    ),
    ("nearest-neighbour", "healpix_H8.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "healpix_H8_nested.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: 1x2 clip: the two nodes of the 10x10 grid inside the ICON-CH1 area, spelled
    # grid [10, 0] like the rotated_ll_subarea.grib2 clip
    ("nearest-neighbour", "icon_ch1.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((1, 2),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [50, 0, 50, 10], "grid": [10, 0]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E: earthkit-data has no gridspec for the Lambert input, so the array path regridded it
    ("nearest-neighbour", "lambert_conformal.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E and §F: earthkit-data has no gridspec for this Lambert input either, so the array
    # path regridded it -- and nearest-neighbour filled the whole globe with values of the North-American input (0 NaN)
    ("nearest-neighbour", "lambert_conformal_na.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E: the array path regridded it, and the encoder now drops the source's scanning mode,
    # so the result spec is the plain target; the input marks 76% of its points missing through complex packing's
    # missing-value management (a valid GRIB2 alternative to a bitmap), which ecCodes decodes as missing and
    # earthkit-data reads as NaN: 591 of the 684 result values are NaN
    ("nearest-neighbour", "mercator.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(591,),
    ),
    # see regrid-grib-findings §E and §F: earthkit-data has no gridspec for this Mercator input either, so the array
    # path regridded it -- and nearest-neighbour filled the whole globe with values of the Pacific input (0 NaN)
    ("nearest-neighbour", "mercator_pacific.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E and §F: earthkit-data has no gridspec for this polar input, so the array path
    # regridded it; all 140 of its 684 NaN are exactly the nodes whose spherical-nearest input point is one of
    # the 73 its bitmap masks -- nearest neighbour, bitmap included, any distance: legitimate NaN (§F)
    ("nearest-neighbour", "polar_stereographic_alaska.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(140,),
    ),
    # see regrid-grib-findings §E: 676 of 684 NaN, and they are the input's own -- 40% of its points are missing
    ("nearest-neighbour", "polar_stereographic_north.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(676,),
    ),
    ("nearest-neighbour", "reduced_gg_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_N32.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_O32.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "reduced_rotated_gg_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "reduced_rotated_gg_N48_global.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "reduced_rotated_gg_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: 6x19 clip of the source area; its 19 NaN are the latitude-90 row, outside the input
    ("nearest-neighbour", "reduced_rotated_gg_subarea_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((6, 19), (6, 19), (6, 19), (6, 19)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
        ),
        nan_counts=(19, 19, 19, 19),
    ),
    # see regrid-grib-findings §D: the same 6x19 clip and the same all-NaN pole row
    ("nearest-neighbour", "reduced_rotated_gg_subarea_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((6, 19), (6, 19), (6, 19), (6, 19)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
        ),
        nan_counts=(19, 19, 19, 19),
    ),
    ("nearest-neighbour", "regular_gg_F16.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "regular_ll.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "regular_ll.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §B and §C: a one-node no-op copy of the input: the target spells its frame; the spec
    # reads grid [0, 0]
    ("nearest-neighbour", "regular_ll_single_point.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((1, 1), (1, 1), (1, 1), (1, 1)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C and §D: the target spells the input sub-area, so MIR echoes its 4x5 nodes instead of
    # the global grid
    ("nearest-neighbour", "regular_ll_subarea.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((4, 5), (4, 5), (4, 5), (4, 5)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "regular_rotated_gg_F48.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): MIR copies the rotated input back unchanged for this
    # same-geometry target; the entry asserts the correct UNROTATED result, so the case strictly xfails
    ("nearest-neighbour", "rotated_ll.grib1", "10x10"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated 10x10 frame and MIR answers with a"
        " no-op copy that keeps the rotation, while the target is unrotated; regrid-grib-findings §C",
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_-180_180_5x5.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_-180_180_5x5.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_0_360_5x5.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_0_360_5x5.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): no-op copy of the rotated 5x5 input for the same-geometry
    # target; the post-fix clip is unknown, so the entry asserts n_fields and an UNROTATED result only
    ("nearest-neighbour", "rotated_ll_subarea.grib1", "10x10"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated frame and MIR answers with a no-op"
        " copy of the rotated 5x5 input, while the target is unrotated; regrid-grib-findings §C",
        n_fields=4,
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
    ),
    # see regrid-grib-findings §D: 1x2 clip, the two target nodes inside the source bbox; the spec spells grid [10, 0]
    ("nearest-neighbour", "rotated_ll_subarea.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((1, 2),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [40, 10, 40, 20], "grid": [10, 0]},),
        nan_counts=(0,),
    ),
    ("nearest-neighbour", "scanning_mode_64.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # --- spectral inputs: fixed 2026-11-20, see the FIXED note above ---
    ("nearest-neighbour", "sh_t32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "sh_t63.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # --- linear onto the 10x10 global lat-lon grid ---
    ("linear", "5x5_multi.grib1", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    ("linear", "5x5_multi.grib2", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §D: the same 2x3 clip as for nearest-neighbour
    ("linear", "alternate-scanning.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((2, 3),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [50, -10, 40, 10], "grid": [10, 10]},),
        nan_counts=(0,),
    ),
    ("linear", "different_grid_resolutions.grib", "10x10"): Expected(
        n_fields=2,
        shapes=((19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §H: the same 36 NaN as for nearest-neighbour -- MIR's silent hop resamples the
    # unstructured input onto O400 with NN before any method runs, so linear shares the NN NaN set -- on the
    # same input whose bitmap masks 44,442 of its 1.74 M points
    ("linear", "eORCA025_T.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(36,),
    ),
    ("linear", "healpix_H8.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "healpix_H8_nested.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: the same 1x2 clip as for nearest-neighbour
    ("linear", "icon_ch1.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((1, 2),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [50, 0, 50, 10], "grid": [10, 0]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E: 682 of 684 NaN: MIR's linear interpolation finds no element among the scattered
    # projected points
    ("linear", "lambert_conformal.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(682,),
    ),
    # see regrid-grib-findings §E and §F: the array path again; 491 of 684 NaN -- the target nodes outside the
    # convex hull of the North-American input's scattered points, where MIR's linear interpolation finds no element
    ("linear", "lambert_conformal_na.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(491,),
    ),
    # see regrid-grib-findings §E and §F: the array path again; 341 of 684 NaN -- the target nodes outside the
    # convex hull of the Pacific input's scattered points, where MIR's linear interpolation finds no element
    ("linear", "mercator_pacific.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(341,),
    ),
    # see regrid-grib-findings §E and §F: the array path again; 522 of 684 NaN -- the target nodes outside the
    # convex hull of this polar input's scattered points, where MIR's linear interpolation finds no element
    ("linear", "polar_stereographic_alaska.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(522,),
    ),
    # see regrid-grib-findings §E: 676 NaN as for nearest-neighbour, from the input's missing points
    ("linear", "polar_stereographic_north.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(676,),
    ),
    ("linear", "reduced_gg_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "reduced_gg_N32.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "reduced_gg_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "reduced_gg_O32.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "reduced_rotated_gg_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "reduced_rotated_gg_N48_global.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "reduced_rotated_gg_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: the same 6x19 clip, with the pole row and 73 NaN
    ("linear", "reduced_rotated_gg_subarea_N32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((6, 19), (6, 19), (6, 19), (6, 19)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
        ),
        nan_counts=(73, 73, 73, 73),
    ),
    # see regrid-grib-findings §D: the same 6x19 clip, with the pole row and 73 NaN
    ("linear", "reduced_rotated_gg_subarea_O32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((6, 19), (6, 19), (6, 19), (6, 19)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
            {"area": [90, -40, 40, 140], "grid": [10, 10]},
        ),
        nan_counts=(73, 73, 73, 73),
    ),
    ("linear", "regular_gg_F16.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "regular_ll.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "regular_ll.grib2", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §B and §C: the same one-node no-op copy
    ("linear", "regular_ll_single_point.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((1, 1), (1, 1), (1, 1), (1, 1)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
            {"area": [50, 10, 50, 10], "grid": [0, 0]},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C and §D: the same 4x5 echo of the input sub-area
    ("linear", "regular_ll_subarea.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((4, 5), (4, 5), (4, 5), (4, 5)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=(
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
            {"area": [60, -10, 30, 30], "grid": [10, 10]},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "regular_rotated_gg_F48.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): the same no-op copy of the rotated input for the same-geometry
    # target; the entry asserts the correct UNROTATED result, so the case strictly xfails
    ("linear", "rotated_ll.grib1", "10x10"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated 10x10 frame and MIR answers with a"
        " no-op copy that keeps the rotation, while the target is unrotated; regrid-grib-findings §C",
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "rotated_ll_-180_180_5x5.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "rotated_ll_-180_180_5x5.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "rotated_ll_0_360_5x5.grib1", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("linear", "rotated_ll_0_360_5x5.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): the same no-op copy of the rotated 5x5 input; the entry asserts
    # n_fields and an UNROTATED result only, the post-fix clip being unknown
    ("linear", "rotated_ll_subarea.grib1", "10x10"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated frame and MIR answers with a no-op"
        " copy of the rotated 5x5 input, while the target is unrotated; regrid-grib-findings §C",
        n_fields=4,
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
    ),
    # see regrid-grib-findings §D: the same 1x2 clip
    ("linear", "rotated_ll_subarea.grib2", "10x10"): Expected(
        n_fields=1,
        shapes=((1, 2),),
        grid_types=("regular_ll",),
        grid_specs=({"area": [40, 10, 40, 20], "grid": [10, 0]},),
        nan_counts=(0,),
    ),
    ("linear", "scanning_mode_64.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # --- spectral inputs: fixed 2026-11-20, see the FIXED note above ---
    ("linear", "sh_t32.grib1", "10x10"): Expected(
        n_fields=4,
        shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
        grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
        grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("linear", "sh_t63.grib", "10x10"): Expected(
        n_fields=1,
        shapes=((19, 36),),
        grid_types=("regular_ll",),
        grid_specs=({"grid": [10, 10]},),
        nan_counts=(0,),
    ),
    # --- nearest-neighbour onto the O32 reduced Gaussian grid ---
    ("nearest-neighbour", "5x5_multi.grib1", "O32"): Expected(
        n_fields=2,
        shapes=((5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0),
    ),
    ("nearest-neighbour", "5x5_multi.grib2", "O32"): Expected(
        n_fields=2,
        shapes=((5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §D: 43 nodes inside the source bbox, spelled N32
    ("nearest-neighbour", "alternate-scanning.grib", "O32"): Expected(
        n_fields=1,
        shapes=((43,),),
        grid_types=("reduced_gg",),
        grid_specs=({"area": [48.8352409662506, -9.47368499999993, 34.8825209937735, 18.947369], "grid": "N32"},),
        nan_counts=(0,),
    ),
    ("nearest-neighbour", "different_grid_resolutions.grib", "O32"): Expected(
        n_fields=2,
        shapes=((5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0),
    ),
    # see regrid-grib-findings §H: 25 of the 5248 nodes NaN, near lats -82 to -88 -- the O32 nodes whose closest
    # node on MIR's intermediate O400 grid (the input is resampled onto it first) falls on the NaN ring left by
    # the input's masked coastal points
    ("nearest-neighbour", "eORCA025_T.grib", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(25,),
    ),
    ("nearest-neighbour", "healpix_H8.grib2", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "healpix_H8_nested.grib2", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D and §G: 13 nodes, the O32 grid clipped to the ICON-CH1 area; the gridspec of the
    # clipped reduced_gg result reads grid "N32" together with that area -- the odd part, §G
    ("nearest-neighbour", "icon_ch1.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((13,),),
        grid_types=("reduced_gg",),
        grid_specs=({"area": [48.8352409662506, 0, 43.2541946653509, 17.142858], "grid": "N32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E: earthkit-data has no gridspec for the Lambert input, so the array path regridded it
    # (its result reads back with the target spec now that the encoder resets the template's grid keys)
    ("nearest-neighbour", "lambert_conformal.grib", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E and §F: the array path regridded this Lambert input too, and nearest-neighbour
    # filled the whole globe -- all 5248 O32 nodes -- with values of the North-American input (0 NaN)
    ("nearest-neighbour", "lambert_conformal_na.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E: earthkit-geo keeps projected inputs on the array path; the input's missing points
    # (complex packing's missing-value management, no bitmap needed) arrive as NaN: 4194 of the 5248 result values
    # are NaN
    ("nearest-neighbour", "mercator.grib", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(4194,),
    ),
    # see regrid-grib-findings §E and §F: the array path regridded this Mercator input too, and nearest-neighbour
    # filled the whole globe -- all 5248 O32 nodes -- with values of the Pacific input (0 NaN)
    ("nearest-neighbour", "mercator_pacific.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §E and §F: the array path regridded this polar input; its 1596 of 5248 NaN are
    # exactly the nodes whose spherical-nearest input point is one of the 73 its bitmap masks -- legitimate NaN,
    # the same rule as on 10x10 (§F)
    ("nearest-neighbour", "polar_stereographic_alaska.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(1596,),
    ),
    # see regrid-grib-findings §E: earthkit-data has no gridspec for the polar input, so the array path regridded it;
    # 5155 of 5248 NaN, from the input's 40% missing points
    ("nearest-neighbour", "polar_stereographic_north.grib1", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(5155,),
    ),
    ("nearest-neighbour", "reduced_gg_N32.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_N32.grib2", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_O32.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "reduced_gg_O32.grib2", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated onto the geographic O32 grid, rotation gone
    ("nearest-neighbour", "reduced_rotated_gg_N32.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated onto the geographic O32 grid, rotation gone
    ("nearest-neighbour", "reduced_rotated_gg_N48_global.grib1", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): MIR copies the rotated input back unchanged for this same-geometry
    # O32 target; the entry asserts the correct UNROTATED result, so the case strictly xfails
    ("nearest-neighbour", "reduced_rotated_gg_O32.grib1", "O32"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated O32 frame and MIR answers with a"
        " no-op copy that keeps the rotation, while the target is unrotated; regrid-grib-findings §C",
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: 488 of the 5248 nodes, clipped to the source area and spelled N32
    ("nearest-neighbour", "reduced_rotated_gg_subarea_N32.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((488,), (488,), (488,), (488,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=(
            {"area": [90, -40, 40.463648178115, 140], "grid": "N32"},
            {"area": [90, -40, 40.463648178115, 140], "grid": "N32"},
            {"area": [90, -40, 40.463648178115, 140], "grid": "N32"},
            {"area": [90, -40, 40.463648178115, 140], "grid": "N32"},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # KNOWN MIR BUG (regrid-grib-findings §C): no-op copy of the rotated input for its own O32 frame; the post-fix
    # clip (§D) is unknown, so the entry asserts n_fields and an UNROTATED result only
    ("nearest-neighbour", "reduced_rotated_gg_subarea_O32.grib1", "O32"): KnownBug(
        reason="KNOWN MIR BUG: the target spells the input's own rotated frame and MIR answers with a no-op"
        " copy of the rotated input, clipped to its area and still rotated, while the target is unrotated;"
        " regrid-grib-findings §C",
        n_fields=4,
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
    ),
    ("nearest-neighbour", "regular_gg_F16.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "regular_ll.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "regular_ll.grib2", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: 100 nodes, clipped to the source area and spelled N32
    ("nearest-neighbour", "regular_ll_subarea.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((100,), (100,), (100,), (100,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=(
            {"area": [59.9970201084913, -10, 32.091943881744, 30], "grid": "N32"},
            {"area": [59.9970201084913, -10, 32.091943881744, 30], "grid": "N32"},
            {"area": [59.9970201084913, -10, 32.091943881744, 30], "grid": "N32"},
            {"area": [59.9970201084913, -10, 32.091943881744, 30], "grid": "N32"},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated onto the geographic O32 grid, rotation gone
    ("nearest-neighbour", "regular_rotated_gg_F48.grib1", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated onto the geographic O32 grid, rotation gone
    ("nearest-neighbour", "rotated_ll.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_-180_180_5x5.grib1", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_-180_180_5x5.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_0_360_5x5.grib1", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §C: interpolated, and the input's rotation is gone
    ("nearest-neighbour", "rotated_ll_0_360_5x5.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # see regrid-grib-findings §D: 972 nodes: the O32 grid clipped to the source bbox and spelled N32
    ("nearest-neighbour", "rotated_ll_subarea.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((972,), (972,), (972,), (972,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=(
            {"area": [90, 0, 40.463648178115, 360], "grid": "N32"},
            {"area": [90, 0, 40.463648178115, 360], "grid": "N32"},
            {"area": [90, 0, 40.463648178115, 360], "grid": "N32"},
            {"area": [90, 0, 40.463648178115, 360], "grid": "N32"},
        ),
        nan_counts=(0, 0, 0, 0),
    ),
    # see regrid-grib-findings §D: 20 nodes inside the source bbox
    ("nearest-neighbour", "rotated_ll_subarea.grib2", "O32"): Expected(
        n_fields=1,
        shapes=((20,),),
        grid_types=("reduced_gg",),
        grid_specs=({"area": [46.0447266311017, 4.5, 34.8825209937735, 21.428572], "grid": "N32"},),
        nan_counts=(0,),
    ),
    ("nearest-neighbour", "scanning_mode_64.grib", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
    # --- spectral inputs: fixed 2026-11-20, see the FIXED note above ---
    ("nearest-neighbour", "sh_t32.grib1", "O32"): Expected(
        n_fields=4,
        shapes=((5248,), (5248,), (5248,), (5248,)),
        grid_types=("reduced_gg", "reduced_gg", "reduced_gg", "reduced_gg"),
        grid_specs=({"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}, {"grid": "O32"}),
        nan_counts=(0, 0, 0, 0),
    ),
    ("nearest-neighbour", "sh_t63.grib", "O32"): Expected(
        n_fields=1,
        shapes=((5248,),),
        grid_types=("reduced_gg",),
        grid_specs=({"grid": "O32"},),
        nan_counts=(0,),
    ),
}

# The form of the lines ``generate_regrid_grib_data.py`` prints for pasting into the table above
# (one line per triple it runs, wrapped here for readability; never read by anything):
#
#     ("nearest-neighbour", "regular_ll.grib2", "10x10"): Expected(
#         n_fields=4,
#         shapes=((19, 36), (19, 36), (19, 36), (19, 36)),
#         grid_types=("regular_ll", "regular_ll", "regular_ll", "regular_ll"),
#         grid_specs=({"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}, {"grid": [10, 10]}),
#         nan_counts=(0, 0, 0, 0),
#     ),

# ---------------------------------------------------------------------------
# the negative set: the files whose regridding must fail, and with what
# ---------------------------------------------------------------------------

# (interpolation, filename, target id) -> Invalid, one entry per GRIB_COMBINATIONS pair of
# every INVALID_GRIB_FILES file: the errors the matrix probe observed on 2026-11-20 for the
# six triples that were BROKEN "raises" entries here until the two files left the matrix, and
# re-observed the same day by hand for this set. Both files contradict themselves, so the
# refusal is the correct behaviour and the message is the finding (regrid-grib-findings §A):
# regular_gg holds Ni*Nj = 192*64 = 12288 for 18432 data points, which ecCodes itself refuses
# to interpret; regular_ll encodes 2 degree increments on a 5 degree grid, which MIR's own
# GRIB reader asserts on. Responsible: the data files, in both cases.
INVALID: dict[tuple[str, str, str], Invalid] = {
    ("nearest-neighbour", "regular_gg_wrong_increment.grib", "10x10"): Invalid(
        exception="gribapi.errors.WrongGridError",
        match="Grid description is wrong or inconsistent",
        reason="INVALID INPUT: gribapi.errors.WrongGridError: Grid description is wrong or inconsistent -- the"
        " file holds Ni*Nj = 12288 for 18432 points, ecCodes refuses the grid; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
    ("nearest-neighbour", "regular_ll_wrong_increment.grib", "10x10"): Invalid(
        exception="RuntimeError",
        match="Assertion failed: ni == ni_",
        reason="INVALID INPUT: RuntimeError: Assertion failed: ni == ni_ -- the file encodes 2 degree"
        " increments on a 5 degree grid, MIR's reader refuses the contradiction; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
    ("linear", "regular_gg_wrong_increment.grib", "10x10"): Invalid(
        exception="gribapi.errors.WrongGridError",
        match="Grid description is wrong or inconsistent",
        reason="INVALID INPUT: gribapi.errors.WrongGridError: Grid description is wrong or inconsistent -- the"
        " file holds Ni*Nj = 12288 for 18432 points, ecCodes refuses the grid; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
    ("linear", "regular_ll_wrong_increment.grib", "10x10"): Invalid(
        exception="RuntimeError",
        match="Assertion failed: ni == ni_",
        reason="INVALID INPUT: RuntimeError: Assertion failed: ni == ni_ -- the file encodes 2 degree"
        " increments on a 5 degree grid, MIR's reader refuses the contradiction; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
    ("nearest-neighbour", "regular_gg_wrong_increment.grib", "O32"): Invalid(
        exception="gribapi.errors.WrongGridError",
        match="Grid description is wrong or inconsistent",
        reason="INVALID INPUT: gribapi.errors.WrongGridError: Grid description is wrong or inconsistent -- the"
        " file holds Ni*Nj = 12288 for 18432 points, ecCodes refuses the grid; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
    ("nearest-neighbour", "regular_ll_wrong_increment.grib", "O32"): Invalid(
        exception="RuntimeError",
        match="Assertion failed: ni == ni_",
        reason="INVALID INPUT: RuntimeError: Assertion failed: ni == ni_ -- the file encodes 2 degree"
        " increments on a 5 degree grid, MIR's reader refuses the contradiction; regrid-grib-findings §A,"
        " responsible: the data file",
    ),
}


def _exception_class(name: str) -> type:
    """Resolve an exception name to its class: a builtin, or module-qualified, imported now.

    Called by :func:`cases` and by the negative test only, so importing this module needs
    neither earthkit-data nor the GRIB stack; collecting the cases does need the module that
    defines a named exception on ``sys.path`` -- which is the case wherever the GRIB tests
    run at all.
    """
    builtin = getattr(builtins, name, None)
    if isinstance(builtin, type) and issubclass(builtin, BaseException):
        return builtin
    module_name, _, class_name = name.rpartition(".")
    if not module_name:
        raise ValueError(f"cannot resolve exception {name!r}: neither a builtin nor a dotted name")
    class_ = getattr(importlib.import_module(module_name), class_name)
    if not (isinstance(class_, type) and issubclass(class_, BaseException)):
        raise ValueError(f"{name} is not an exception class")
    return class_


def cases():
    """Yield one ``pytest.param`` per (interpolation, file, target) triple of the matrix.

    The value is the triple itself, in :data:`GRIB_COMBINATIONS` order and, within a
    combination, in :data:`GRIB_FILES` order. The id is ``<interpolation>-<filename>-
    <target id>``; file names contain ``-`` and ``.``, which pytest ids take as they are.
    A triple listed in :data:`EXPECTED` as :class:`Broken` carries its xfail mark:
    ``xfail(raises=<exception class>, strict=True)`` for kind "raises" (the case runs and
    must raise exactly that), ``xfail(run=False)`` for kind "hang" (never executed, because
    it hangs or aborts the interpreter). A :class:`KnownBug` triple carries
    ``xfail(strict=True)`` with no ``raises=`` and no ``run=False``: the case runs, fails
    against the correct result the entry records while MIR is buggy, and XPASSes as a failure
    once MIR is fixed. An ``Expected`` or missing entry runs unmarked.
    """
    keys = {
        (interpolation, filename, target_id)
        for interpolation, target_id in GRIB_COMBINATIONS
        for filename in GRIB_FILES
    }
    for key in sorted(set(EXPECTED) - keys):
        raise ValueError(f"EXPECTED names {key}, which is not a (combination, file) triple of this table")

    for interpolation, target_id in GRIB_COMBINATIONS:
        for filename in GRIB_FILES:
            key = (interpolation, filename, target_id)
            marks = []
            entry = EXPECTED.get(key)
            if isinstance(entry, Broken):
                if entry.kind == "raises":
                    marks.append(
                        pytest.mark.xfail(raises=_exception_class(entry.exception), strict=True, reason=entry.reason)
                    )
                else:
                    marks.append(pytest.mark.xfail(run=False, reason=entry.reason))
            elif isinstance(entry, KnownBug):
                marks.append(pytest.mark.xfail(strict=True, reason=entry.reason))
            yield pytest.param(key, id=f"{interpolation}-{filename}-{target_id}", marks=marks)


def invalid_cases():
    """Yield one ``pytest.param`` per (interpolation, file, target) triple of the negative set.

    The values are the three parts of the triple -- interpolation, filename and target id --
    in :data:`GRIB_COMBINATIONS` order and, within a combination, in
    :data:`INVALID_GRIB_FILES` order; the id is ``<interpolation>-<filename>-<target id>``,
    the form of :func:`cases`. Every triple carries its :class:`Invalid` entry in
    :data:`INVALID`, and no mark: the test itself asserts the error with ``pytest.raises``.
    Unlike :func:`cases` this resolves nothing -- the exception class named by an entry is
    imported by the test, through :func:`_exception_class`, when it runs -- so collecting the
    negative cases works without earthkit-data, on a machine that cannot run them.
    """
    keys = {
        (interpolation, filename, target_id)
        for interpolation, target_id in GRIB_COMBINATIONS
        for filename in INVALID_GRIB_FILES
    }
    for key in sorted(set(INVALID) - keys):
        raise ValueError(f"INVALID names {key}, which is not a (combination, invalid file) triple of this table")
    for key in sorted(keys - set(INVALID)):
        raise ValueError(f"INVALID has no entry for {key}: every negative file needs one per combination")

    for interpolation, target_id in GRIB_COMBINATIONS:
        for filename in INVALID_GRIB_FILES:
            key = (interpolation, filename, target_id)
            yield pytest.param(
                interpolation, filename, target_id, id=f"{interpolation}-{filename}-{target_id}", marks=[]
            )
