# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Spherical harmonics: the ``sh`` "grid" and its ``T<T>`` names.

``sh`` is not a real grid: a ``T<T>`` spec describes a triangular spectral
truncation, i.e. the coefficients of a spherical harmonics expansion, and has no
points (``docs/concepts/grid/types/spectral.rst``). Metadata still works, point
operations do not; both are asserted here exactly as observed (probed in
``.ecmwf-multi-agents/notes/grid-probe.md``, re-verified for these tests). The
errors raised by the point operations are EXPECTED behaviour for a grid without
points, not bugs.

The one real defect is :attr:`Grid.order`: reading it lets an uncaught C++
``eckit::NotImplemented`` escape the bindings and aborts the interpreter (SIGABRT,
uncatchable from Python). It can therefore never be read inside a test run; see
:class:`TestSHOrder`.
"""

import re
import signal
import subprocess
import sys
from math import prod

import pytest

from earthkit.geo.grids import Grid

pytestmark = pytest.mark.long_test

# Truncations that build and give correct metadata (all probed; T1 works too, T0 is
# rejected because the truncation must be positive).
TRUNCATIONS = [10, 21, 42, 63, 85, 159, 319, 639]

# The truncation the point operations were probed on (notes/grid-probe.md, sh table).
PROBED_TRUNCATION = 21

# Whole-earth placeholder. ``bounding_box()`` does not raise for ``sh`` and does not
# describe an extent either: it is the nominal global area (probe finding F5).
WHOLE_EARTH = (90.0, 0.0, -90.0, 360.0)

# Exception texts, as observed.
NOT_GRIDDABLE = "GridError: [SphericalHarmonics does not support gridded functionality]"
NO_LAT = "Not implemented: Grid: lat() is not implemented for type 'sh'"
NO_LON = "Not implemented: Grid: lon() is not implemented for type 'sh'"
NO_BOXES = "Serious bug: Representation::gridBoxes() not implemented for SphericalHarmonics"
NO_REORDER = "Not implemented: Grid: reorder() is not implemented for type 'sh'"

# The ways of building the same T<T> grid, all equivalent.
WAYS = ["name", "truncation", "positional", "spec"]


def build(way: str, truncation: int) -> Grid:
    """Build ``T<truncation>`` in one of the equivalent ways."""
    name = f"T{truncation}"
    if way == "name":
        return Grid(grid=name)
    if way == "truncation":
        return Grid(type="sh", truncation=truncation)
    if way == "positional":
        return Grid(name)
    if way == "spec":
        return Grid({"grid": name})
    raise AssertionError(f"unknown build way: {way}")


def expected_size(truncation: int) -> int:
    """Number of real spherical harmonics coefficients of a triangular truncation.

    Degrees ``n = 0..T`` carry ``n + 1`` orders each, so ``(T + 1)(T + 2) / 2``
    complex coefficients, i.e. twice as many real numbers. That is the
    ``(T + 1)(T + 2)`` of ``docs/concepts/grid/types/spectral.rst``.
    """
    return 2 * sum(n + 1 for n in range(truncation + 1))


class TestSHSpec:
    """The spec of an ``sh`` grid, and its canonical ``T<T>`` name."""

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_canonical_spec_is_the_name(self, truncation):
        for way in WAYS:
            assert build(way, truncation).spec == {"grid": f"T{truncation}"}

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_type(self, truncation):
        for way in WAYS:
            assert build(way, truncation).type == "sh"

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_all_ways_are_equal(self, truncation):
        grids = [build(way, truncation) for way in WAYS]
        assert all(grid == grids[0] for grid in grids)

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_round_trip_through_the_spec(self, truncation):
        grid = build("name", truncation)
        assert Grid(grid.spec) == grid

    def test_truncations_differ(self):
        a, b = Grid(grid="T42"), Grid(grid="T63")
        assert a != b
        assert a.spec != b.spec
        assert a.uid != b.uid
        assert a.size() != b.size()


class TestSHMetadata:
    """Metadata that works for ``sh``, checked against the truncation formula."""

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_shape_is_the_flat_coefficient_count(self, truncation):
        assert build("name", truncation).shape == (expected_size(truncation),)

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_size_matches_the_formula(self, truncation):
        assert expected_size(truncation) == (truncation + 1) * (truncation + 2)  # docs' closed form
        grid = build("name", truncation)
        assert grid.size() == expected_size(truncation)
        assert grid.size() == prod(grid.shape)

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_bounding_box_is_the_whole_earth_placeholder(self, truncation):
        assert build("name", truncation).bounding_box() == WHOLE_EARTH

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_spec_str(self, truncation):
        assert build("name", truncation).spec_str == f'{{"grid":"T{truncation}"}}'

    @pytest.mark.parametrize("truncation", TRUNCATIONS)
    def test_uid_is_stable(self, truncation):
        # Same uid for every build way, and for a repeated build in the same process.
        # (It is also stable across processes: these uids match the ones written down
        # in the probe notes, which ran in a separate interpreter.)
        uids = {build(way, truncation).uid for way in WAYS} | {build("name", truncation).uid}
        assert len(uids) == 1, uids
        (uid,) = uids
        assert re.fullmatch(r"[0-9a-f]{32}", uid), uid

    def test_uids_differ_between_truncations(self):
        uids = {build("name", truncation).uid for truncation in TRUNCATIONS}
        assert len(uids) == len(TRUNCATIONS)

    def test_projection_and_figure_are_reachable(self):
        # Not meaningful for sh (there is nothing to project): eckit hands back its
        # default projection and figure. Checked to exist, not for their values.
        grid = Grid(grid="T63")
        assert grid.projection is not None
        assert grid.figure is not None


class TestSHUnsupportedPointOperations:
    """Point operations: ``sh`` has no points, so these raise.

    Expected behaviour, not a bug: there is nothing to answer them with. Recorded
    here with the exact exception class (``RuntimeError`` carrying the eckit error
    text) and message. Every member below works on a point grid such as
    ``Grid(grid=[1, 1])`` (control-checked while these tests were written).

    The probe notes cover ``to_latlons``, ``to_points``, ``lat``, ``lon``,
    ``distinct_latitudes``, ``distinct_longitudes`` and ``grid_box_areas``; the rest
    (``to_xyz``, ``point_coordinates``, ``to_unstructured_ll``, ``search``, ``crop``,
    ``cropped``, ``x``, ``y``, ``reorder``) are recorded by these tests.
    """

    # fmt: off
    @pytest.mark.parametrize(
        "operation, message",
        [
            pytest.param(lambda g: g.to_latlons(), NOT_GRIDDABLE, id="to_latlons"),
            pytest.param(lambda g: g.to_points(), NOT_GRIDDABLE, id="to_points"),
            pytest.param(lambda g: g.to_xyz(), NOT_GRIDDABLE, id="to_xyz"),
            pytest.param(lambda g: g.point_coordinates, NOT_GRIDDABLE, id="point_coordinates"),
            pytest.param(lambda g: g.to_unstructured_ll(), NOT_GRIDDABLE, id="to_unstructured_ll"),
            pytest.param(lambda g: g.search(), NOT_GRIDDABLE, id="search"),
            pytest.param(lambda g: g.crop({"area": [70, -20, 30, 20]}), NOT_GRIDDABLE, id="crop"),
            pytest.param(lambda g: g.cropped({"area": [70, -20, 30, 20]}), NOT_GRIDDABLE, id="cropped"),
            pytest.param(lambda g: g.lat(), NO_LAT, id="lat"),
            pytest.param(lambda g: g.distinct_latitudes(), NO_LAT, id="distinct_latitudes"),
            pytest.param(lambda g: g.y(), NO_LAT, id="y"),
            pytest.param(lambda g: g.lon(), NO_LON, id="lon"),
            pytest.param(lambda g: g.distinct_longitudes(), NO_LON, id="distinct_longitudes"),
            pytest.param(lambda g: g.x(), NO_LON, id="x"),
            pytest.param(lambda g: g.reorder("i+j+"), NO_REORDER, id="reorder"),
            pytest.param(lambda g: g.grid_box_areas(), NO_BOXES, id="grid_box_areas"),
        ],
    )
    # fmt: on
    def test_operation_is_not_supported(self, operation, message):
        # T21, the truncation the extra members were probed on; the messages name the
        # type, not the truncation, so they do not depend on it (except for
        # grid_box_areas(), checked separately below).
        grid = Grid(grid=f"T{PROBED_TRUNCATION}")
        with pytest.raises(RuntimeError, match=re.escape(message)):
            operation(grid)

    def test_grid_box_areas_message_names_the_truncation(self):
        for truncation in (10, PROBED_TRUNCATION):
            message = f"{NO_BOXES}[truncation={truncation}]"
            with pytest.raises(RuntimeError, match=re.escape(message)):
                Grid(grid=f"T{truncation}").grid_box_areas()


class TestSHSpecErrors:
    """``truncation`` is required and must be positive (reference/grids.rst)."""

    @pytest.mark.parametrize(
        "kwargs, message",
        [
            ({"type": "sh"}, "SpecError: [truncation]"),
            ({"type": "sh", "truncation": 0}, "SpecError: [SphericalHarmonics: truncation must be positive]"),
            ({"type": "sh", "truncation": -1}, "SpecError: [Bad cast: Negative value cannot be cast to unsigned type]"),
        ],
    )
    def test_invalid_spec(self, kwargs, message):
        with pytest.raises(RuntimeError, match=re.escape(message)):
            Grid(**kwargs)

    # The name form is deliberately not checked: `Grid(grid="T0")` is rejected with
    # the misleading "Grid: cannot build grid without 'type'" (probe finding F6),
    # after eckit has dumped a list of all factories to stderr.


class TestSHOrder:
    """:attr:`Grid.order` aborts the interpreter for ``sh`` (the one real defect)."""

    @pytest.mark.xfail(
        run=False,
        reason="BROKEN: reading .order aborts the interpreter (SIGABRT, uncaught C++ eckit::NotImplemented)",
    )
    def test_order_is_not_implemented(self):
        # Kept as a marker so the gap is listed by `-E long` without ever risking the
        # run: run=False means this body is never executed. It would abort Python.
        assert Grid(grid="T63").order == "i+j-"

    def test_order_abort_confined_to_a_child_process(self):
        # What actually documents the behaviour: reading .order kills the process.
        # Deliberately done in a child that is allowed to die. If eckit starts
        # reporting this as a catchable Python exception this test fails, which is the
        # signal to turn the xfail above into a normal test.
        code = (
            "from earthkit.geo.grids import Grid\n"
            "grid = Grid(grid='T10')\n"
            "print('BEFORE-ORDER', flush=True)\n"
            "print(grid.order)\n"
            "print('AFTER-ORDER', flush=True)\n"
        )
        result = subprocess.run(
            [sys.executable, "-u", "-c", code],
            capture_output=True,
            text=True,
            timeout=120,
        )
        reached = [line for line in result.stdout.splitlines() if line in ("BEFORE-ORDER", "AFTER-ORDER")]
        assert reached == ["BEFORE-ORDER"], (
            f"expected the child to die while reading .order, it printed {reached} and exited "
            f"with {result.returncode}; stderr: {result.stderr[-300:]}"
        )
        assert result.returncode == -signal.SIGABRT, (
            f"expected the child to abort (SIGABRT), exit status was {result.returncode}; "
            f"stderr: {result.stderr[-300:]}"
        )
