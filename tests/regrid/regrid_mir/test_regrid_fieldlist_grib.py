# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""earthkit-data GRIB fieldlist regridding tests (``earthkit.geo.regrid``) against MIR references.

The GRIB half of the old ``test_regrid_fieldlist.py``; its non-GRIB (array) half is
``test_regrid_fieldlist_array.py``. What the two files share -- the fieldlist loader, the
10x10 reference check and the ``..._to_ll`` case table and test body -- is
``regrid_fieldlist_common.py``, while the shared pieces they have with the rest of the
folder, the kwarg alias parametrisation, the ``mir`` / ``earthkit-data`` skip marks and
the loading of the legacy npz references, still come from ``regrid_common``.

The old use cases are unchanged: a regular lat-lon and an octahedral Gaussian GRIB input
regridded to a 10x10 global lat-lon grid, as a fieldlist and as a single field (the
``field_type="grib"`` cases of the parametrised ``..._to_ll`` test), the deprecated
``grid=`` kwarg, the grib1 and grib2 lat-lon -> points -> lat-lon round trips (one
parametrised test) and the check that regridding to points does not make ecCodes print a
GridUnknownError. Their reference comparisons deliberately keep the strict
``compare_global_ll_results`` helper of ``earthkit.geo.utils.testing`` instead of the
registry compare of ``regrid_common.get_interpolation(...).compare`` (used by the matrix
tests): these checks keep the tolerance and the pole-row handling they always had. The GRIB
inputs are read from the test-data URL with ``get_test_data_path`` (a URL for
``from_source``, not a local file, so it is not a ``regrid_common.data_path`` case).

Two speeds
----------
*fast* -- everything above, the checks of the GRIB table and of the comparison rules that
need no MIR result (``test_regrid_fieldlist_grib_table``,
``test_regrid_fieldlist_grib_case_marks``,
``test_regrid_fieldlist_grib_nearest_ties_in_coordinates``,
``test_regrid_fieldlist_grib_nearest_values_are_input_values`` and
``test_regrid_fieldlist_grib_nearest_packing_allowance``), one smoke test
(``test_regrid_fieldlist_grib_smoke_nearest_10x10``) that regrids ``regular_ll.grib2`` and
``reduced_gg_O32.grib1`` from the server with property checks only: no reference file, so
it says something even before the GRIB references below are on the server; the two
regressions of the spectral-input fix -- a spectral field with no explicit ``in_grid``
goes straight to the backend's GRIB path --
``test_regrid_fieldlist_grib_spectral_nearest_10x10`` -- ``sh_t32.grib1`` onto the 10x10
grid, property checks only -- and ``test_regrid_fieldlist_grib_spectral_needs_grib_backend``
-- a stub backend without ``regrid_grib`` must refuse spectral input with a clear ValueError;
and the §E regression ``test_regrid_fieldlist_grib_projected_array_path`` -- a projected
``lambert_conformal.grib`` whose field is given a grid spec must stay on the array path,
shown with a spy backend and checked with property checks only.
Of these, all but the marks, the smoke test and the three regressions are dependency-free --
the table check and the three rule checks hold without ``mir`` and without earthkit-data,
which is what keeps the default suite saying something about the matrix on a machine that
cannot run it; only ``test_regrid_fieldlist_grib_case_marks`` (``cases()`` resolves the
exception classes of the BROKEN entries), the smoke test and the three regressions need the
GRIB stack and skip when it is missing (``test_regrid_fieldlist_grib_spectral_nearest_10x10``
and the §E regression also need ``mir``, the backend test does not).

*long* (``-E long``, marked ``long_test``) -- the matrix ``test_regrid_fieldlist_grib_matrix``,
the file check ``test_regrid_fieldlist_grib_files`` and the negative set
``test_regrid_fieldlist_grib_invalid_file``.

The matrix
----------
40 of the 42 GRIB files of the download server's ``test-data/grib/`` folder (``GRIB_FILES``;
the two whose grid metadata contradicts itself are the negative set below), regridded
whole-fieldlist as it is read by ``from_source("file", ...)``, with each combination of
``GRIB_COMBINATIONS`` (nearest-neighbour or linear onto the 10x10 global
lat-lon grid, nearest-neighbour onto ``{"grid": "O32"}``): 120 cases, one per
(interpolation, file, target) triple, parametrised by ``regrid_grib_cases.cases()`` with the
id ``<interpolation>-<file>-<target>`` and the xfail mark of a BROKEN or KnownBug triple.

What one triple does is recorded in ``regrid_grib_cases.EXPECTED`` -- MIR's own result, as
observed by ``generate_regrid_grib_data.py``, including the results that look odd (a regional
source that clips the target to its area, a result that is NaN over most of the target because
its input is mostly missing values). A triple MIR cannot produce is a ``Broken`` entry:
``xfail(raises=<exception>, strict=True)`` when it raises, ``xfail(run=False)`` when it hangs
or aborts; five are, the four one-point triples that fail inside MIR and the mercator linear
case too slow to run. Six triples produce a result that is wrong: the verified MIR bug of
regrid-grib-findings §C (a rotated input regridded onto the
bare target of its own geometry comes back unchanged, rotation kept) is a ``KnownBug`` entry
recording the correct, unrotated result and marked ``xfail(strict=True)``; the case runs and
fails against that correct result while MIR is buggy, and its stored reference npz is the
stale no-op output, so the reference comparison is skipped for a KnownBug triple. Per case
the test checks, against the table and the reference file: the number of
fields; per field the value shape, the ``gridType``, ``geography.grid_spec()``, the NaN
count and the values; that ``param``, ``level``, ``date`` and ``time`` are the input's; and
for nearest-neighbour that every non-NaN value of a field is one of that input field's values
(MIR copies, it never invents -- exact membership, with room only for what the two GRIB
messages' own packing can explain, see :func:`_assert_values_from_input`; not checked for a
spectral input, whose array holds coefficients, not grid values).
Nothing of the input is read before the regridding: unpacking a file's values can fail with an
exception of its own, and the table records the regridding's -- which is why the negative test
below reads nothing either.

The negative set
----------------
The two files the matrix does not take -- ``regular_gg_wrong_increment.grib``, whose Ni*Nj
contradicts its number of data points, and ``regular_ll_wrong_increment.grib``, 2 degree
increments on a 5 degree grid (regrid-grib-findings §A, responsible: the data files) -- are
not regridded but asserted to fail: ``test_regrid_fieldlist_grib_invalid_file`` runs the same
3 combinations on each of them (6 cases, ids as in the matrix) and checks the error with
``pytest.raises`` against the exception type and message its ``regrid_grib_cases.INVALID``
entry records. They have no reference npz, and the generator, which iterates ``GRIB_FILES``,
never probes them.

The data
--------
Two file trees. ``data_path("legacy", <file>, subfolder="grib")`` finds the GRIB **inputs**
(``URL_ROOT/grib/`` on the download server, cached in ``tests/data/grib/``, git-ignored); that
kind always comes from the server, the override below does not touch it.
``data_path("grib_reference", interpolation, <file>, target_id)``, i.e.
``<root>/regrid/grib/<interpolation>/<file>/<target>.npz`` with the values of all result
fields stacked into ``arr_0``, holds MIR's reference (cached in ``tests/data/regrid/grib/``,
downloaded from ``URL_ROOT/regrid/grib/`` when missing).
``EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data`` redirects only that second, new-hierarchy
kind: the references are then read from the local tree and never downloaded, which is how
freshly generated ones are tested before the upload (``tests/regrid/regrid_mir/README.md``,
"The GRIB matrix"); the inputs still come from the server cache. Without the override and
without the upload the matrix cases fail with a ``FileNotFoundError`` naming the generator;
the fast tests do not read references, so they are unaffected.

A nearest result is compared strictly (``get_interpolation("nearest-neighbour").compare``:
NaN at the same positions, every finite value equal) and, where that fails, only values that
belong to an equally near source point are accepted, by
``regrid_common.assert_nearest_ties`` -- the tie a nearest search resolves differently on
another machine. Here both grids are passed as ``(lats, lons)`` coordinates taken from the
fields themselves (``field.geography.latlons()``): a tie is a question about the points, and the
coordinate mode asks nothing else -- so it works for every grid of the matrix, the projected
inputs included, whether or not their message answers a gridspec. A linear result skips its
pole rows only when the table
records the whole global 10x10 grid, never a clipped one. Every long case also carries
``pytest.mark.timeout(300, method="thread")``, so a MIR hang ends the run with a message
instead of stalling it.
"""

import datetime
import math
import os
from pathlib import Path

import numpy as np
import pytest
from regrid_common import (
    CACHE_ROOT,
    KWARG_ALIASES,
    REGRID_TEST_DATA_ENV,
    assert_nearest_ties,
    data_path,
    get_interpolation,
    legacy_values,
    load_values,
    requires_ekd,
    requires_mir,
)
from regrid_fieldlist_common import (
    METADATA_KEYS,
    METADATA_KEYS_GRID_TYPE,
    TO_LL_CASES,
    check_regrid_to_ll,
    create_test_fieldlist,
    run_to_ll,
)
from regrid_grib_cases import (
    EXPECTED,
    GRIB_COMBINATIONS,
    GRIB_FILES,
    GRIB_TARGETS,
    INVALID,
    INVALID_GRIB_FILES,
    Broken,
    Expected,
    Invalid,
    KnownBug,
    _exception_class,
    cases,
    invalid_cases,
)

from earthkit.geo import regrid
from earthkit.geo.grids._regrid.data.fieldlist.grib import PROJECTED_GRID_TYPES
from earthkit.geo.utils.testing import NO_EKD

if not NO_EKD:
    from earthkit.data import create_fieldlist, from_source

# The six target points of the lat-lon -> points regriddings.
_POINTS_LATS = [40.0, 40.0, 40.0, 30.0, 30.0, 30.0]
_POINTS_LONS = [10.0, 20.0, 30.0, 5.0, 25.0, 35.0]
_POINTS_NUM = 6

# Metadata of the two fields of the 5x5_multi GRIB messages, in earthkit-data names.
_ROUND_TRIP_METADATA_KEYS = ["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"]
_ROUND_TRIP_METADATA_REF = [
    ["2t", 2, datetime.datetime(2024, 3, 23, 12, 0), datetime.timedelta(0)],
    ["2t", 2, datetime.datetime(2024, 3, 24, 0, 0), datetime.timedelta(seconds=43200)],
]

# ---------------------------------------------------------------------------
# the GRIB matrix: the helpers its one long test and the fast checks share
# ---------------------------------------------------------------------------

# target id -> GribTarget, and the 120 (interpolation, file, target) triples of the matrix in
# table order; the 6 of the negative set are parametrised straight from invalid_cases().
_TARGET_BY_ID = {target.id: target for target in GRIB_TARGETS}
_TRIPLES = [
    (interpolation, filename, target_id) for interpolation, target_id in GRIB_COMBINATIONS for filename in GRIB_FILES
]

# The shape of the global 10x10 lat-lon grid, and the gridspec that is exactly that grid.
_GLOBAL_10X10_SHAPE = (19, 36)
_GLOBAL_10X10_SPEC = {"grid": [10, 10]}

# Relative size of one float32 ulp: the precision the values of these messages are handled in
# (``_assert_values_from_input``). The packing allowance is on top of it, not instead of it.
_VALUE_ULP = float(np.finfo(np.float32).eps)

# The packings whose ``binaryScaleFactor`` / ``decimalScaleFactor`` really are the step between
# two values the message can hold: simple, complex and complex-with-spatial-differencing packing
# all round onto ``referenceValue + X * 2**binaryScaleFactor * 10**-decimalScaleFactor``, the two
# complex ones group at a time but on that same lattice (see :func:`_packing_quantum`). The lossy
# compressions use the two keys differently and are not here: the eight grid_ccsds fields of the
# matrix deviate from the value they copied by less than half a float32 ulp and the two grid_jpeg
# ones copy exactly, so neither needs a packing allowance at all.
_QUANTUM_PACKINGS = ("grid_simple", "grid_complex", "grid_complex_spatial_differencing")

# The two files of the fast smoke test: one regular lat-lon, one reduced Gaussian, both with a
# clean nearest-neighbour -> 10x10 triple in regrid_grib_cases.EXPECTED.
_SMOKE_FILES = ["regular_ll.grib2", "reduced_gg_O32.grib1"]


def _grib_fieldlist(filename):
    """The server GRIB file *filename*, read the way the matrix and the generator read it.

    ``data_path("legacy", ..., subfolder="grib")`` fetches it into ``tests/data/grib/``
    (git-ignored) and ``from_source("file", ...)`` gives the GRIB fieldlist, so the input is
    the same file the references were generated from.
    """
    path = data_path("legacy", filename, subfolder="grib")
    return from_source("file", str(path)).to_fieldlist()


def _matrix_params():
    """The parameters of the matrix test: ``regrid_grib_cases.cases()``, or one placeholder.

    ``cases()`` resolves the exception class each BROKEN "raises" triple names, and such a
    class may come from ``gribapi``, ie need earthkit-data to be importable (since the two
    §A files left the matrix, today's BROKEN entries name builtins only). A machine without
    earthkit-data cannot run the matrix at all -- it must still be able to import this module
    for its fast tests, so the collection failure becomes a single parameter, which
    ``requires_ekd`` skips with the usual reason. The negative test needs no such guard:
    ``invalid_cases()`` resolves nothing.
    """
    if NO_EKD:
        return [pytest.param(_TRIPLES[0], id="no-earthkit-data")]
    return list(cases())


def _assert_gridspec(actual, expected, path="gridspec"):
    """Compare a field's returned gridspec with what the table recorded; numbers get a tolerance.

    Recurses through dicts (same key set) and lists / tuples (same length); int and float
    values are compared with ``math.isclose`` (rel_tol 1e-9), because MIR recomputes some area
    bounds and may spell them with more digits than the table records. Anything else must
    compare equal, and ``None`` -- a field that answers no gridspec at all -- matches ``None``
    only. No entry of today's table records one: the array path's re-encoded message carries the
    target spec whatever its input grid was, and a triple whose correct post-fix spec nothing can
    observe leaves its ``grid_specs`` out of the table entirely, which the matrix test handles
    without asking this helper. The branch is here in case a result stops answering a spec.
    """
    if expected is None:
        assert actual is None, f"{path}: {actual!r} is not None, as the table records"
    elif isinstance(expected, dict):
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


def _is_global_10x10(spec) -> bool:
    """True when a returned gridspec is exactly the whole global 10x10 lat-lon grid.

    The pole-row exemption of the linear comparison is for that grid alone: a result clipped
    to the source area spells something else (an ``area``, fewer points, a kept rotation), any
    other spec carries keys besides ``grid`` -- a projected one spells a ``projection`` with its
    x/y extents -- and a field that answers no spec at all says ``None``, so nothing is ever
    exempted for them.
    """
    if not isinstance(spec, dict) or set(spec) != set(_GLOBAL_10X10_SPEC):
        return False
    grid = spec["grid"]
    return isinstance(grid, (list, tuple)) and [float(x) for x in grid] == [10.0, 10.0]


def _field_points(field, where):
    """The ``(lats, lons)`` of a field's own points, the coordinate mode of the tie check.

    A regridded GRIB field is a message, so its coordinates come from ecCodes through
    ``geography.latlons()``. That works for every result of today's matrix: the array path
    re-encodes its target grid into a message ecCodes reads back, projected inputs included.
    Should a result ever come back
    unreadable, the tie check has no target grid to work with, so the difference from the
    reference stays unexplained: say which case that is instead of reporting an ecCodes
    traceback, and let the strict comparison's own error remain the reason the case fails.
    """
    try:
        lats, lons = field.geography.latlons()
    except Exception as error:
        raise AssertionError(
            f"{where}: its values differ from the reference, and the grid of this result cannot be read back"
            " (its latitudes and longitudes are not there), so the difference cannot be checked for ties"
        ) from error
    assert lats is not None and lons is not None, f"{where}: the field has no latitudes/longitudes to tie-check"
    return (np.asarray(lats), np.asarray(lons))


def _packing_quantum(field):
    """The step between two values a GRIB field's packing can hold, None when it has none.

    Simple, complex and complex-with-spatial-differencing packing store ``value = referenceValue
    + X * 2**binaryScaleFactor * 10**-decimalScaleFactor`` with an integer ``X``, so that factor
    is the step between two values the message can hold; the two complex ones write their values
    group by group, but every group starts on the same lattice. Measured over every message of
    the matrix that reports a step -- simple, complex or differenced -- all of its values do sit
    on it (largest residual 4e-15 of a step for simple and complex, 3.5e-13 for the two
    grid_complex_spatial_differencing files added on 2026-11-20, float64 arithmetic).
    ``decimalScaleFactor`` is 0 in most of them, 1 in ``mercator.grib``, ``mercator_pacific.grib2``
    and ``polar_stereographic_alaska.grib2`` and 2 in two fields of a one-point file, and it
    scales the same way. ``None`` for a field whose packing is none of those three -- the lossy
    compressions use the keys differently, see :data:`_QUANTUM_PACKINGS` -- and for a field that
    does not answer at all, which is the case of an array written outside a GRIB message.
    """
    try:
        if field.metadata("packingType") not in _QUANTUM_PACKINGS:
            return None
        binary = field.metadata("binaryScaleFactor")
        decimal = field.metadata("decimalScaleFactor")
        return float(2.0 ** float(binary) * 10.0 ** (-float(decimal)))
    except Exception:  # noqa: BLE001
        # any answer but the numbers above means "this field has no value step to report"
        return None


def _packing_allowance(*fields):
    """Half the coarsest value step among ``fields``: what a re-packing may have moved a value.

    The encoder rounds onto its own lattice, so a value copied into a new message is at most half
    a step away from the one it came from; the step of the field it came from is taken into the
    same bound, because that value was itself rounded once when its file was written. ``0.0`` when
    no field reports a step, which leaves the float32 ulp margin of
    :func:`_assert_values_from_input` as the only allowance.
    """
    quanta = [quantum for quantum in (_packing_quantum(field) for field in fields) if quantum is not None]
    return 0.5 * max(quanta) if quanta else 0.0


def _assert_values_from_input(values, input_values, where, allowance=0.0):
    """Fail unless every non-NaN value of ``values`` is one of ``input_values``.

    The reference-free rule of a nearest-neighbour result: MIR copies source values around, it
    never computes new ones (the positions it NaN-fills outside a regional source do not
    count). **Exact membership is the rule**; the only deviation accepted is the one the two
    GRIB messages can account for, ``allowance`` -- half the packing value step of the coarser
    of them, :func:`_packing_allowance` -- plus one float32 ulp of the value, the precision the
    numbers are handled in. When the coarser step is the input's own, the allowance cannot reach a
    second input value, because the distinct values of a field are that step apart or more. A result
    packed more coarsely than its input could in principle reach one, so this holds for the matrix
    as measured, not in general -- and it is the situation the two
    grid_complex_spatial_differencing files are in, see the fourth bullet below. Measured over the
    164 nearest result fields of the matrix whose two messages report a packing step and whose
    input holds more than one value, half the coarser step is never more than half the spacing of
    the input values (max and median ratio 0.5, reached by every field whose input values are
    exactly one step apart, ``polar_stereographic_north.grib1``, ``mercator.grib``,
    ``eORCA025_T.grib`` and ``icon_ch1.grib2`` among them; the last digit of those ratios is
    float32 noise) except for those four Pacific and Alaska fields, whose result message packs
    onto a lattice of 0.25 / 0.125 over the inputs' 0.1 spacing (ratios 1.25 and 0.625). The 8
    fields of the two grid_ccsds files report no step and get the ulp alone, so do the 2 fields of
    the grid_jpeg file, which copy exactly, and the 4 of the one-point input are out of that
    ratio because a single value has no spacing.

    Measured on this tree (2026-11-20, all 178 result fields of the 73 non-spectral nearest
    triples of the matrix -- a spectral input is exempt from this rule, its array holds
    spherical-harmonic coefficients, not grid values): 168 of them hold
    values that are *exactly* input values, and every deviation of the ten that do not is either
    the float32 precision of a value or the lattice of the message that holds it --

    * ``mercator.grib`` onto ``10x10`` and onto ``O32``, 4.30312502 for the input's 4.3 and 8.1875
      for its 8.2: both messages are grid_complex, so a copied value sits on the lattice of the
      result, 0.4 of its own step (``binaryScaleFactor`` -7 onto 10x10, -5 onto O32) -- a sixteenth
      and a quarter of the allowance, whose half step comes from the input file's coarser complex
      packing (``binaryScaleFactor`` 0 with ``decimalScaleFactor`` 1: a step of 0.1);
    * ``reduced_rotated_gg_subarea_N32.grib1`` / ``_O32.grib1`` onto ``10x10``, 5.722e-06 and
      7.629e-06 at values of 82.43 and 71.36: 0.58 and 0.90 ulp; their grid_simple step is
      0.03125, so the packing allowance is 0.0156 -- 2000 times more than the deviation;
    * ``reduced_gg_N32.grib2`` / ``_O32.grib2`` onto ``10x10``, 7.276e-12 at a value of
      4.859e-04: 0.13 and 0.12 ulp; grid_ccsds, so again no packing allowance.
    * ``mercator_pacific.grib2`` and ``polar_stereographic_alaska.grib2``, both onto ``10x10`` and
      onto ``O32``: every value sits 0.049994-0.049997 off the input value it copied. Both inputs
      are grid_complex_spatial_differencing with ``binaryScaleFactor`` 0 and ``decimalScaleFactor``
      1 (a 0.1 lattice), and their result message re-packs differenced with ``binaryScaleFactor``
      -2 (step 0.25, the Pacific file) or -3 (step 0.125, the Alaska file) -- 0.2 and 0.4 of the
      result's own step, inside the 0.125 / 0.0625 half step of the coarser of the two messages
      that the allowance is;

    So the ulp margin covers the two grid_ccsds fields, the complex lattice covers the two
    mercator ones (at 0.4 of a step), the differenced lattice covers the four Pacific and Alaska
    ones (at 0.2 and 0.4 of a step) and the simple one covers the two rotated sub-area ones (at a
    deviation two thousand times smaller than the half step their packing allows); the half-step
    bound is room for a re-packing onto a lattice as coarse as the coarser of the two messages --
    which is exactly what the Pacific and Alaska results do. A nearest
    result that needs more than that is a finding, not something to widen the bound for.
    """
    values = np.asarray(values).reshape(-1)
    input_values = np.asarray(input_values).reshape(-1)
    finite = values[~np.isnan(values)]
    if finite.size == 0:
        return
    sources = np.unique(input_values[~np.isnan(input_values)])
    if sources.size == 0:
        raise AssertionError(f"{where}: the input holds no value a result could come from")
    sources.sort()
    # the closest input value to each result value, by binary search in the sorted uniques
    position = np.searchsorted(sources, finite)
    left = sources[np.clip(position - 1, 0, sources.size - 1)]
    right = sources[np.clip(position, 0, sources.size - 1)]
    take_left = np.abs(finite - left) <= np.abs(finite - right)
    matched = np.where(take_left, left, right)
    offset = np.where(take_left, np.abs(finite - left), np.abs(finite - right))
    bound = allowance + _VALUE_ULP * np.maximum(np.abs(finite), np.abs(matched))
    copies = offset <= bound
    bad = int((~copies).sum())
    assert bad == 0, (
        f"{where}: {bad} result values that are not values of the input, e.g. "
        f"{finite[~copies][:5]} (nearest input value {matched[~copies][:5]}, further apart than "
        f"the {allowance:g} packing allowance and one float32 ulp)"
    )


# ---------------------------------------------------------------------------
# the old use cases, unchanged
# ---------------------------------------------------------------------------


@requires_mir
@requires_ekd
@pytest.mark.parametrize("_kwarg,interpolation", KWARG_ALIASES)
@pytest.mark.parametrize("case", TO_LL_CASES)
def test_regrid_fieldlist_grib_to_ll(case, _kwarg, interpolation):
    run_to_ll(case, _kwarg, interpolation, field_type="grib")


@requires_mir
@requires_ekd
def test_regrid_fieldlist_grib_deprecated_grid_kwarg():
    ds = create_test_fieldlist("5x5.grib", field_type="grib")

    interpolation = "linear"

    v_ref = legacy_values(f"out_5x5_10x10_{interpolation}.npz")
    metadata_ref = ds.metadata(list(METADATA_KEYS_GRID_TYPE))

    r = regrid(ds, grid={"grid": [10, 10]}, interpolation=interpolation)

    check_regrid_to_ll(
        r, v_ref, interpolation, METADATA_KEYS_GRID_TYPE, metadata_ref, is_field=False, check_grid_type=False
    )


@requires_mir
@requires_ekd
# The old grib2 test had a commented-out parametrisation over ["grib", "array"]; the
# array field type is not covered by these round trips, so it stays "grib" only.
@pytest.mark.parametrize(
    "filename",
    [pytest.param("5x5_multi.grib1", id="grib1"), pytest.param("5x5_multi.grib2", id="grib2")],
)
def test_regrid_fieldlist_grib_ll_to_points_round_trip(filename):
    ds = create_test_fieldlist(filename, subfolder="grib", field_type="grib")

    # Regrid from lat-lon to points
    lats = list(_POINTS_LATS)
    lons = list(_POINTS_LONS)
    out_grid_points = {"latitudes": lats, "longitudes": lons}

    r = regrid(ds, out_grid=out_grid_points, interpolation="nn")

    # mir turns grib1 to grib2 internally for unstructured input grids
    ref_vals = np.array([
        [288.44410706, 286.67457581, 281.93043518, 303.16481018, 292.76637268, 293.47926331],
        [289.20581055, 278.41088867, 271.97143555, 294.8894043, 285.50463867, 283.16674805],
    ])

    metadata_keys = list(_ROUND_TRIP_METADATA_KEYS)

    assert len(r) == 2

    for i, f in enumerate(r):
        assert f.shape == (_POINTS_NUM,)
        lats_res, lons_res = f.geography.latlons()
        assert np.allclose(lats_res, lats)
        assert np.allclose(lons_res, lons)
        assert np.allclose(f.values, ref_vals[i])

    assert r.get(metadata_keys) == _ROUND_TRIP_METADATA_REF

    # Encode each regridded field to a new GRIB message and create a new fieldlist from them
    fields = []
    for f in r:
        f_grib = from_source("memory", f.sync().message()).to_fieldlist()[0]
        fields.append(f_grib)

    ds_r = create_fieldlist(fields)

    assert len(ds_r) == 2
    for i, f in enumerate(ds_r):
        assert f.shape == (_POINTS_NUM,)
        # A GRIB message cannot hold the points of an unstructured grid, so the field that comes
        # back from it has no coordinates at all and latlons() answers (None, None).
        lats_res, lons_res = f.geography.latlons()
        assert lats_res is None
        assert lons_res is None
        assert np.allclose(f.values, ref_vals[i])

    assert ds_r.get(metadata_keys) == _ROUND_TRIP_METADATA_REF
    assert ds_r.metadata("edition") == [2, 2]

    # Regrid from points back to lat-lon
    out_grid_ll = {"grid": [5, 5]}

    r_ll = regrid(ds_r, in_grid=out_grid_points, out_grid=out_grid_ll, interpolation="nn")

    assert len(r_ll) == 2
    for f in r_ll:
        assert f.shape == (37, 72)
        assert f.geography.grid_spec().items() >= out_grid_ll.items()
        lats_res, lons_res = f.geography.latlons()
        assert lats_res.shape == (37, 72)
        assert lons_res.shape == (37, 72)


def _ekd_requests_grid_spec_for_zero_uid(message_bytes):
    """Return True if the installed earthkit-data asks ecCodes for the "gridSpec" of a zero-uid grid.

    Such a request makes ecCodes/eckit print "GridUnknownError" to stderr. The check is
    behaviour-based, so it also works for local development versions of earthkit-data.
    """

    class _Spy:
        def __init__(self, handle):
            self._handle = handle
            self.keys = []

        def get(self, name, *args, **kwargs):
            self.keys.append(name)
            return self._handle.get(name, *args, **kwargs)

        def __getattr__(self, name):
            return getattr(self._handle, name)

    try:
        from earthkit.data.field.grib.geography import GribGeographyBuilder
        from earthkit.data.readers.grib.handle import MemoryGribHandle

        handle = MemoryGribHandle.from_message(message_bytes)

        # the input must be an unstructured grid with an all-zero grid uid
        assert handle.get("gridType") == "unstructured_grid"
        uid = handle.get("uuidOfHGrid")
        uid = uid.hex() if isinstance(uid, bytes) else str(uid).strip()
        assert uid != "" and set(uid) == {"0"}

        spy = _Spy(handle)
        GribGeographyBuilder.build(spy)
    except (ImportError, AttributeError):
        return True

    return "gridSpec" in spy.keys


@requires_mir
@requires_ekd
@pytest.mark.parametrize("filename", ["5x5_multi.grib1", "5x5_multi.grib2"])
def test_regrid_fieldlist_grib_ll_to_points_no_grid_unknown_error(filename, capfd):
    ds = create_test_fieldlist(filename, subfolder="grib")

    lats = list(_POINTS_LATS)
    lons = list(_POINTS_LONS)

    # Regrid to points, encode and read back the result. The stderr output is captured throughout.
    r = regrid(ds, out_grid={"latitudes": lats, "longitudes": lons}, interpolation="nn")

    assert len(r) == 2
    messages = []
    for f in r:
        assert f.values.shape == (_POINTS_NUM,)
        lats_res, lons_res = f.geography.latlons()
        assert np.allclose(lats_res, lats)
        assert np.allclose(lons_res, lons)
        messages.append(f.sync().message())

    for msg in messages:
        f_grib = from_source("memory", msg).to_fieldlist()[0]
        assert f_grib.geography.shape() == (_POINTS_NUM,)

    stderr = capfd.readouterr().err

    # Only run the probe after the output has been read so that its own stderr cannot change the result.
    forbidden = ["GridUnknownError", "GridSpec: geo::Exception"]
    if any(s in stderr for s in forbidden) and _ekd_requests_grid_spec_for_zero_uid(messages[0]):
        pytest.xfail("earthkit-data without the zero grid-uid fix prints GridUnknownError")

    for s in forbidden:
        assert s not in stderr, stderr


# ---------------------------------------------------------------------------
# the fast end of the matrix: the table itself, and the two rules the stored references
# cannot demonstrate on their own
# ---------------------------------------------------------------------------


def test_regrid_fieldlist_grib_table():
    """Every matrix triple has an entry and every negative triple a reason; both tables are complete.

    Fast on purpose: the table is what the long cases are planned from, and a hole in it must
    show up in the default run, not only when the references are there. Dependency-free on
    purpose too -- it reads :data:`regrid_grib_cases.EXPECTED`, :data:`regrid_grib_cases.INVALID`
    and nothing else, because without
    earthkit-data :func:`regrid_grib_cases.cases` may not be able to resolve the exception classes
    of the BROKEN entries and the matrix runs on one placeholder parameter
    (``_matrix_params``). The
    marks themselves are the subject of ``test_regrid_fieldlist_grib_case_marks``, which skips
    when that dependency is missing.
    """
    assert len(GRIB_FILES) == len(set(GRIB_FILES)), "GRIB_FILES lists a file twice"
    assert len(GRIB_TARGETS) == len({target.id for target in GRIB_TARGETS}), "GRIB_TARGETS lists an id twice"
    targets = {target.id for target in GRIB_TARGETS}
    for interpolation, target_id in GRIB_COMBINATIONS:
        assert get_interpolation(interpolation).name == interpolation, (
            f"{interpolation}: not the canonical name of a registered interpolation"
        )
        assert target_id in targets, f"{target_id}: not a target of GRIB_TARGETS"
    assert len(set(GRIB_COMBINATIONS)) == len(GRIB_COMBINATIONS), "GRIB_COMBINATIONS lists a pair twice"

    assert set(_TRIPLES) == set(EXPECTED), (
        f"{len(set(EXPECTED))} entries for {len(_TRIPLES)} triples: "
        f"missing {sorted(set(_TRIPLES) - set(EXPECTED))[:5]}, unknown {sorted(set(EXPECTED) - set(_TRIPLES))[:5]}"
    )

    for key, entry in EXPECTED.items():
        if isinstance(entry, Broken):
            assert entry.reason.startswith("BROKEN: "), f"{key}: reason {entry.reason!r}"
            if entry.kind == "raises":
                assert entry.exception, f"{key}: a BROKEN raises entry without an exception name"
        elif isinstance(entry, KnownBug):
            assert entry.reason.startswith("KNOWN MIR BUG: "), f"{key}: reason {entry.reason!r}"
            assert "regrid-grib-findings" in entry.reason, f"{key}: a KnownBug reason must cite the report"
            assert entry.n_fields >= 1, f"{key}: {entry!r}"
        else:
            assert isinstance(entry, Expected), f"{key}: {entry!r} is neither Expected, Broken nor KnownBug"
            assert entry.n_fields >= 1, f"{key}: {entry!r}"

    # the negative set: the files it tests are not in the matrix, and each of them has one
    # entry per combination, naming the error and the message the regridding must fail with
    assert len(INVALID_GRIB_FILES) == len(set(INVALID_GRIB_FILES)), "INVALID_GRIB_FILES lists a file twice"
    shared = set(GRIB_FILES) & set(INVALID_GRIB_FILES)
    assert not shared, f"{sorted(shared)}: in GRIB_FILES and INVALID_GRIB_FILES, a file is either matrix or negative"
    negative = {
        (interpolation, filename, target_id)
        for interpolation, target_id in GRIB_COMBINATIONS
        for filename in INVALID_GRIB_FILES
    }
    assert set(INVALID) == negative, (
        f"{len(INVALID)} entries for {len(negative)} negative triples: "
        f"missing {sorted(negative - set(INVALID))[:5]}, unknown {sorted(set(INVALID) - negative)[:5]}"
    )
    for key, entry in INVALID.items():
        assert isinstance(entry, Invalid), f"{key}: {entry!r} is not an Invalid entry"
        assert entry.exception, f"{key}: a negative entry without an exception name"
        assert entry.match, f"{key}: a negative entry without a message match"
        assert entry.reason.startswith("INVALID INPUT: "), f"{key}: reason {entry.reason!r}"
        assert "regrid-grib-findings" in entry.reason, f"{key}: a negative reason must cite the report"

    # the ids pytest shows, spelled here rather than through cases() so that this holds without
    # earthkit-data; cases() is checked against them by the test below
    ids = [f"{interpolation}-{filename}-{target_id}" for interpolation, filename, target_id in _TRIPLES]
    ids += [f"{interpolation}-{filename}-{target_id}" for interpolation, filename, target_id in sorted(INVALID)]
    assert len(set(ids)) == len(ids), f"duplicate case ids: {sorted({i for i in ids if ids.count(i) > 1})}"


@requires_ekd
def test_regrid_fieldlist_grib_case_marks():
    """``cases()`` walks the table in order and carries the xfail mark of every BROKEN and known-bug triple.

    Marked ``requires_ekd`` because ``cases()`` resolves the exception class each BROKEN "raises"
    entry names, and such a class may come from ``gribapi`` (earthkit-data's GRIB stack); the
    BROKEN entries left today name builtins only. No other test looks at the ids and marks the
    matrix is parametrised with.
    """
    params = list(cases())
    assert [param.values[0] for param in params] == _TRIPLES, "cases() does not walk the table in order"
    ids = [param.id for param in params]
    assert ids == [f"{interpolation}-{filename}-{target_id}" for interpolation, filename, target_id in _TRIPLES], (
        "cases() does not build the <interpolation>-<file>-<target> id of every triple"
    )
    for param in params:
        entry = EXPECTED[param.values[0]]
        xfail = [mark for mark in param.marks if mark.name == "xfail"]
        if isinstance(entry, Broken):
            assert len(xfail) == 1, f"{param.values[0]}: {len(xfail)} xfail marks, one expected"
            if entry.kind == "raises":
                assert xfail[0].kwargs.get("raises") is not None, f"{param.values[0]}: xfail without raises="
                assert xfail[0].kwargs.get("strict") is True, f"{param.values[0]}: a raises= xfail must be strict"
            else:
                assert xfail[0].kwargs.get("run") is False, f"{param.values[0]}: a hang must not be executed"
        elif isinstance(entry, KnownBug):
            assert len(xfail) == 1, f"{param.values[0]}: {len(xfail)} xfail marks, one expected"
            assert xfail[0].kwargs.get("strict") is True, f"{param.values[0]}: a known-bug xfail must be strict"
            assert "raises" not in xfail[0].kwargs, f"{param.values[0]}: a known-bug xfail must not name an exception"
            assert "run" not in xfail[0].kwargs, f"{param.values[0]}: a known-bug case must really run"
        else:
            assert not xfail, f"{param.values[0]}: a triple the table marks neither BROKEN nor KnownBug is xfailed"


# A source with two points equally near a target point, given as coordinates because that is
# the mode the GRIB matrix needs: (0, 0) is as far from (0, 5) as (0, 10) is, so either of
# their values is a legitimate nearest result, while the value of (10, 5) -- 10 degrees away
# -- is not. No gridspec is involved anywhere: the points are the only thing the check sees,
# exactly as for a GRIB field on a projection grid.
_TIE_LATS = (np.array([0.0, 0.0, 10.0, -10.0]), np.array([0.0, 10.0, 5.0, 5.0]))
_TIE_INPUT = np.array([1.0, 2.0, 3.0, 4.0])  # the value of each source point
_TIE_LATS_TARGET = (np.array([0.0, 0.0]), np.array([5.0, 0.0]))
_TIE_REFERENCE = np.array([1.0, 1.0])  # both target points took source point 0


def test_regrid_fieldlist_grib_nearest_ties_in_coordinates():
    """``assert_nearest_ties`` in coordinate mode: a tie passes, a far point and a NaN do not."""

    def tie(result):
        return assert_nearest_ties(result, _TIE_REFERENCE, _TIE_INPUT, _TIE_LATS, _TIE_LATS_TARGET, context="ties")

    # the strict rule really does fail on the other point of the tie pair, which is what the
    # fallback is for, and the tie check then accepts it
    other = np.array([2.0, 1.0])
    with pytest.raises(AssertionError):
        get_interpolation("nearest-neighbour").compare(other, _TIE_REFERENCE)
    tie(other)

    # the value of a source point 10 degrees away is not a tie -- the "every result value
    # occurs in the input" check cannot see that difference, a wrong point's value *is* an
    # input value -- and neither is a value the input does not hold at all
    for wrong in ([3.0, 1.0], [4.0, 1.0], [1.5, 1.0]):
        with pytest.raises(AssertionError):
            tie(np.array(wrong))

    # a changed NaN mask is never a tie: it is the wrong coverage
    with pytest.raises(AssertionError):
        tie(np.array([np.nan, 1.0]))
    with pytest.raises(AssertionError):
        assert_nearest_ties(_TIE_REFERENCE, np.array([1.0, np.nan]), _TIE_INPUT, _TIE_LATS, _TIE_LATS_TARGET)

    # and coordinates that do not hold one point per value are refused, not silently truncated
    with pytest.raises(AssertionError):
        assert_nearest_ties(other, _TIE_REFERENCE, _TIE_INPUT, _TIE_LATS[:1], _TIE_LATS_TARGET)


# A GRIB field that only answers the three packing keys _packing_quantum asks for, so the
# tolerance of the membership rule can be checked on numbers of the test's own choosing. The two
# values of the "input" are 10 apart: everything the rule accepts in between it accepts because a
# packing says so, never because the gap is wide.
_PACKING_INPUT = np.array([0.0, 10.0])


class _PackedField:
    """A stand-in for a GRIB field: ``metadata()`` answers what it was given and raises otherwise."""

    def __init__(self, keys):
        self._keys = keys

    def metadata(self, key):
        return self._keys[key]


def _simple_packing(binary, decimal=0):
    """A field with simple packing: value step 2**binary * 10**-decimal, as ecCodes reports it."""
    return _PackedField({"packingType": "grid_simple", "binaryScaleFactor": binary, "decimalScaleFactor": decimal})


def _complex_packing(binary, decimal=0):
    """The same field, complex packed: the same value step, stored one group at a time."""
    return _PackedField({"packingType": "grid_complex", "binaryScaleFactor": binary, "decimalScaleFactor": decimal})


def _differenced_packing(binary, decimal=0):
    """The same field once more, complex packed with spatial differencing.

    The values are coded group by group as differences, but they sit on the very same lattice --
    measured over the ``mercator_pacific.grib2`` and ``polar_stereographic_alaska.grib2``
    messages of the matrix, all of their values are within 3.5e-13 of a step of it. It is also
    what the regridding's encoder keeps for those two inputs, with a coarser
    ``binaryScaleFactor`` than the input's own.
    """
    return _PackedField({
        "packingType": "grid_complex_spatial_differencing",
        "binaryScaleFactor": binary,
        "decimalScaleFactor": decimal,
    })


# The two values of one of the complex-packed mercator.grib messages, and the value the array
# path's result holds where the input had 4.3. That result is packed with referenceValue 0.6 and a
# step of 2**-7, and 4.3 is 473.6 steps from its reference, so the copy is stored as
# 0.6 + 474 * 2**-7 = 4.303125: 0.4 of a step away from the value it copied. _assert_values_from_input
# has to accept that, and 4.29 -- the copy of nothing -- still has to fail under the result's own
# allowance (the two-message allowance, 0.05 from the input's 0.1 step, would accept it).
_COMPLEX_INPUT = np.array([4.3, 4.4])
_COMPLEX_COPY = 4.30312502


def test_regrid_fieldlist_grib_nearest_values_are_input_values():
    """``_assert_values_from_input`` is exact membership, with one float32 ulp of margin only.

    The rule without a file and without MIR: what a nearest result may hold. Half the distance
    between two input values is deliberately one of the things it may not accept -- an invented
    value in the middle of the range is not a copy of anything, however near an input value it is.
    """
    # a copy passes, including where the source is NaN (that is coverage, not a value)
    _assert_values_from_input(np.array([0.0, 10.0, np.nan]), _PACKING_INPUT, "copies")
    _assert_values_from_input(np.array([np.nan, np.nan]), _PACKING_INPUT, "all NaN")

    # and everything else is an invented value
    for wrong in (0.1, 4.0, 6.0, 5.0, -0.1, 10.1):
        with pytest.raises(AssertionError, match="not values of the input"):
            _assert_values_from_input(np.array([wrong]), _PACKING_INPUT, f"invented {wrong}")
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([0.0, 4.0]), np.array([10.0]), "one input value")

    # an input that holds nothing but NaN cannot be the source of anything
    with pytest.raises(AssertionError, match="no value"):
        _assert_values_from_input(np.array([1.0]), np.array([np.nan]), "no input value")

    # the margin really is one float32 eps: a deviation below it is the precision the numbers are
    # handled in (the two grid_ccsds fields of the matrix deviate by 0.13 and 0.12 eps at their
    # largest difference, and that is all the room they get), eight times it is a different value
    _assert_values_from_input(np.array([1.0 + 5e-8]), np.array([1.0]), "half an eps")
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([1.0 + 1e-6]), np.array([1.0]), "eight eps")


def test_regrid_fieldlist_grib_nearest_packing_allowance():
    """The only extra allowance is half the value step of the message's lattice; no step, no allowance."""
    # binaryScaleFactor -4: the message can hold steps of 0.0625, so a copy may sit 0.03125 away
    fine = _simple_packing(-4)
    assert _packing_quantum(fine) == 0.0625
    assert _packing_allowance(fine) == 0.03125
    assert _packing_allowance(fine, _PackedField({})) == 0.03125
    _assert_values_from_input(np.array([0.03]), _PACKING_INPUT, "packed 0.03", _packing_allowance(fine))
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([0.05]), _PACKING_INPUT, "packed 0.05", _packing_allowance(fine))

    # the coarser of the two fields decides: a step of 0.5 allows a quarter, not a thirty-second
    coarse = _simple_packing(-1)
    assert _packing_allowance(fine, coarse) == 0.25
    _assert_values_from_input(np.array([0.2]), _PACKING_INPUT, "packed 0.2", _packing_allowance(coarse))
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([0.4]), _PACKING_INPUT, "packed 0.4", _packing_allowance(coarse))

    # complex packing rounds onto the same lattice, and it is what the array path inherits when it
    # re-encodes mercator.grib: binaryScaleFactor -7 onto 10x10 and -5 onto O32, so the copy of the
    # input's 4.3 is 4.30312502 -- 0.4 of a step away, inside the half step that -7 allows
    complex_result = _complex_packing(-7)
    assert _packing_quantum(complex_result) == 0.0078125
    assert _packing_allowance(complex_result) == 0.00390625
    _assert_values_from_input(
        np.array([_COMPLEX_COPY]), _COMPLEX_INPUT, "complex copy", _packing_allowance(complex_result)
    )
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([4.29]), _COMPLEX_INPUT, "invented 4.29", _packing_allowance(complex_result))
    # and the coarser message decides here too: the mercator.grib input is complex packed with
    # binaryScaleFactor 0 and decimalScaleFactor 1, a step of 0.1 whose half still stops short of
    # the next input value, 0.1 away
    assert _packing_allowance(complex_result, _complex_packing(0, decimal=1)) == 0.05

    # spatial differencing is complex packing coded group by group as differences, on the same
    # lattice, and the regridding's encoder keeps it: the 0.1 lattice of the inputs of
    # mercator_pacific.grib2 and polar_stereographic_alaska.grib2 (bsf 0, dsf 1) comes back packed
    # bsf -2 (a 0.25 lattice) and bsf -3 (0.125), so a copied value sits half an input step --
    # 0.05 -- off, inside the half step of the coarser of the two messages; 0.15 is the copy of nothing
    differenced = _differenced_packing(-2)
    assert _packing_quantum(differenced) == 0.25
    assert _packing_allowance(differenced) == 0.125
    _assert_values_from_input(np.array([0.05]), _PACKING_INPUT, "differenced copy", _packing_allowance(differenced))
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(np.array([0.15]), _PACKING_INPUT, "invented 0.15", _packing_allowance(differenced))
    # the coarser of the two decides here too, and it can be either of them
    assert _packing_allowance(differenced, _differenced_packing(0, decimal=1)) == 0.125
    assert _packing_allowance(_differenced_packing(-3), _differenced_packing(0, decimal=1)) == 0.0625

    # a packing whose keys are not a value step gives no allowance: the lossy compressions use the
    # two scale factors their own way, the eight grid_ccsds fields of the matrix are within half a
    # float32 ulp of the values they copied and the two grid_jpeg ones copy exactly, so none of
    # them needs an allowance at all
    for keys in (
        {"packingType": "grid_ccsds", "binaryScaleFactor": -17, "decimalScaleFactor": 0},
        {"packingType": "grid_jpeg", "binaryScaleFactor": -5, "decimalScaleFactor": 0},
        {"packingType": "grid_simple", "binaryScaleFactor": None, "decimalScaleFactor": 0},
        {"packingType": "grid_simple", "binaryScaleFactor": -4},
        {},
    ):
        assert _packing_quantum(_PackedField(keys)) is None, keys
    assert _packing_allowance(_PackedField({}), _simple_packing(0)) == 0.5
    assert _packing_allowance() == 0.0
    with pytest.raises(AssertionError, match="not values of the input"):
        _assert_values_from_input(
            np.array([0.1]), _PACKING_INPUT, "no step to point at", _packing_allowance(_PackedField({}))
        )


@requires_mir
@requires_ekd
@pytest.mark.parametrize("filename", _SMOKE_FILES)
def test_regrid_fieldlist_grib_smoke_nearest_10x10(filename):
    """The fast end of the matrix on two server files: properties only, no reference file.

    A regular lat-lon and a reduced Gaussian GRIB input, regridded as a fieldlist to the 10x10
    global lat-lon grid with nearest-neighbour -- both triples clean in
    ``regrid_grib_cases.EXPECTED``. The checks are the ones of the long case with the numbers
    the table records written out here instead: field count, shape, gridType, gridspec, no
    NaN, the input metadata kept, and every result value copied from the input. Nothing is
    read from the reference tree, so this test passes before the GRIB references are uploaded.
    """
    ds = _grib_fieldlist(filename)
    metadata_ref = ds.metadata(METADATA_KEYS)

    out = regrid(ds, out_grid=_GLOBAL_10X10_SPEC, interpolation="nearest-neighbour")

    assert len(out) == len(ds), f"{filename}: {len(out)} result fields for {len(ds)} input fields"
    for index, field in enumerate(out):
        where = f"{filename} field {index}"
        values = np.asarray(field.to_numpy())
        assert values.shape == _GLOBAL_10X10_SHAPE, f"{where}: shape {values.shape}"
        assert field.metadata("gridType") == "regular_ll", f"{where}: gridType {field.metadata('gridType')!r}"
        _assert_gridspec(field.geography.grid_spec(), _GLOBAL_10X10_SPEC, f"{where} gridspec")
        assert not np.isnan(values).any(), f"{where}: {int(np.isnan(values).sum())} NaN values"
        _assert_values_from_input(values, np.asarray(ds[index].to_numpy()), where, _packing_allowance(field, ds[index]))

    assert out.metadata(METADATA_KEYS) == metadata_ref, f"{filename}: the regridding changed param, level, date or time"


@requires_mir
@requires_ekd
def test_regrid_fieldlist_grib_spectral_nearest_10x10():
    """Fast regression: a spectral GRIB fieldlist regrids onto 10x10.

    ``sh_t32.grib1`` holds spherical-harmonic coefficients with no grid points, and
    earthkit-geo used to ask such a field for its latitudes/longitudes before handing the
    message to MIR -- which aborted inside earthkit-data on the None it had no answer for
    (``latlons()`` answers ``(None, None)`` for a grid without coordinates now). The message goes to
    MIR's GRIB path now. The checks are those of ``test_regrid_fieldlist_grib_smoke_nearest_10x10``
    with the numbers the table records for this triple, minus the values-from-input rule: a
    spectral field's array holds coefficients, not grid values, so the result values cannot
    be members of it.
    """
    filename = "sh_t32.grib1"
    ds = _grib_fieldlist(filename)
    metadata_ref = ds.metadata(METADATA_KEYS)

    out = regrid(ds, out_grid=_GLOBAL_10X10_SPEC, interpolation="nearest-neighbour")

    assert len(out) == len(ds) == 4, f"{filename}: {len(out)} result fields for {len(ds)} input fields"
    for index, field in enumerate(out):
        where = f"{filename} field {index}"
        values = np.asarray(field.to_numpy())
        assert values.shape == _GLOBAL_10X10_SHAPE, f"{where}: shape {values.shape}"
        assert field.metadata("gridType") == "regular_ll", f"{where}: gridType {field.metadata('gridType')!r}"
        _assert_gridspec(field.geography.grid_spec(), _GLOBAL_10X10_SPEC, f"{where} gridspec")
        assert not np.isnan(values).any(), f"{where}: {int(np.isnan(values).sum())} NaN values"

    assert out.metadata(METADATA_KEYS) == metadata_ref, f"{filename}: the regridding changed param, level, date or time"


@requires_ekd
def test_regrid_fieldlist_grib_spectral_needs_grib_backend():
    """A backend without ``regrid_grib`` refuses spectral input with a clear ValueError.

    The array path cannot regrid a spectral field: its values are coefficients, not grid
    points. The stub backend stands in for any backend that only regrids arrays, and the
    error is raised before any regridding, so this needs no MIR.
    """
    from earthkit.geo.grids._regrid.backends import Backend

    class _ArrayOnlyBackend(Backend):
        def regrid(self, values, in_grid, out_grid, **kwargs):
            raise AssertionError("the array path must not be reached for a spectral field")

    ds = _grib_fieldlist("sh_t32.grib1")
    with pytest.raises(ValueError, match="[Ss]pectral.*regrid_grib"):
        regrid(ds, out_grid=_GLOBAL_10X10_SPEC, interpolation="nearest-neighbour", backend=_ArrayOnlyBackend())


# The GRIB gridSpec ecCodes computes for a projected message whose projection names the
# figure of the earth, spelled like the one earthkit-data hands out for mercator.grib (see
# regrid-grib-findings §E), but with the 475x475 extent and the Lambert conformal projection
# of the file the test below regrids.
_PROJECTED_SPEC = {
    "projection": {
        "type": "lcc",
        "figure": "wgs84_sphere",
        "lat_0": 39.0,
        "lon_0": -98.0,
        "lat_1": 33.0,
        "lat_2": 45.0,
    },
    "grid": [475, 475],
    "x": [-2031639.2, 2031639.2, 8583.0],
    "y": [-1471943.1, 2605000.0, 8583.0],
    "false": [0.0, 0.0],
    "order": "i+-j+",
}


@requires_mir
@requires_ekd
def test_regrid_fieldlist_grib_projected_array_path(monkeypatch):
    """Fast regression of regrid-grib-findings §E: a projected field keeps the array path.

    earthkit-data can answer a grid spec like :data:`_PROJECTED_SPEC` for a projected
    GRIB field. With it the field would look structured and go to MIR's GRIB path, which
    1.30.2 cannot run for projected messages; ``regrid_grib_field`` therefore builds the
    input grid of a projected field (GRIB gridType lambert, polar_stereographic,
    mercator or lambert_azimuthal_equal_area) from its latitudes/longitudes when no
    ``in_grid`` is given, and the field stays on the array path. The spec is faked here
    so the test says the same whatever the installed earthkit-data answers for a
    projected grid, and the spy backend records which of the two paths the regridding
    took. The checks are those of ``test_regrid_fieldlist_grib_smoke_nearest_10x10``
    with the numbers the table records for ``lambert_conformal.grib`` onto 10x10.
    """
    from earthkit.data.field.grib.geography import GribGeography

    from earthkit.geo.grids._regrid.backends.mir import MirBackend

    class _PathSpyBackend(MirBackend):
        """The MIR backend, counting the calls of its two regridding paths."""

        def __init__(self):
            super().__init__()
            self.grib_calls = 0
            self.array_calls = 0

        def regrid_grib(self, *args, **kwargs):
            self.grib_calls += 1
            return super().regrid_grib(*args, **kwargs)

        def regrid(self, *args, **kwargs):
            self.array_calls += 1
            return super().regrid(*args, **kwargs)

    # the fake answers the projected spec for a projected field only, so the regridded
    # regular_ll result keeps the spec it really has
    real_grid_spec = GribGeography.grid_spec

    def fake_grid_spec(self):
        if self.handle.get("gridType", default=None) in PROJECTED_GRID_TYPES:
            return dict(_PROJECTED_SPEC)
        return real_grid_spec(self)

    monkeypatch.setattr(GribGeography, "grid_spec", fake_grid_spec)

    filename = "lambert_conformal.grib"
    ds = _grib_fieldlist(filename)
    assert ds[0].geography.grid_spec() is not None, "the projected field must be given a grid spec to route on"
    metadata_ref = ds.metadata(METADATA_KEYS)

    backend = _PathSpyBackend()
    out = regrid(ds, out_grid=_GLOBAL_10X10_SPEC, interpolation="nearest-neighbour", backend=backend)

    assert backend.grib_calls == 0, f"{filename}: the projected field was sent to the GRIB path"
    assert backend.array_calls == len(ds), (
        f"{filename}: {backend.array_calls} array-path calls for {len(ds)} projected fields"
    )

    assert len(out) == len(ds) == 1, f"{filename}: {len(out)} result fields for {len(ds)} input fields"
    for index, field in enumerate(out):
        where = f"{filename} field {index}"
        values = np.asarray(field.to_numpy())
        assert values.shape == _GLOBAL_10X10_SHAPE, f"{where}: shape {values.shape}"
        assert field.metadata("gridType") == "regular_ll", f"{where}: gridType {field.metadata('gridType')!r}"
        _assert_gridspec(field.geography.grid_spec(), _GLOBAL_10X10_SPEC, f"{where} gridspec")
        assert not np.isnan(values).any(), f"{where}: {int(np.isnan(values).sum())} NaN values"
        _assert_values_from_input(values, np.asarray(ds[index].to_numpy()), where, _packing_allowance(field, ds[index]))

    assert out.metadata(METADATA_KEYS) == metadata_ref, f"{filename}: the regridding changed param, level, date or time"


# ---------------------------------------------------------------------------
# the long matrix: one case per (interpolation, file, target) triple
# ---------------------------------------------------------------------------


def _expected_grib_npz_files():
    """The npz files the GRIB table expects: one per triple that is not BROKEN.

    The generator writes no reference for a triple MIR cannot produce, and its report JSON is
    kept outside the ``regrid/`` tree, so this is the complete list of what
    ``<root>/regrid/grib/`` holds. The KnownBug triples are in it: their stored npz is the
    stale no-op output of the §C MIR bug, kept until MIR is fixed, and the matrix test never
    reads it for those triples.
    """
    return {
        f"regrid/grib/{interpolation}/{filename}/{target_id}.npz"
        for interpolation, filename, target_id in _TRIPLES
        if not isinstance(EXPECTED.get((interpolation, filename, target_id)), Broken)
    }


@pytest.mark.long_test
def test_regrid_fieldlist_grib_files():
    """The GRIB npz files on disk match the table: the exact set for a local tree, none extra for the cache.

    The GRIB half of what ``test_regrid_numpy_grids.py`` checks for the numpy matrix (that
    check ignores this subtree): with ``EARTHKIT_GEO_REGRID_TEST_DATA`` set the folder must
    hold **exactly** one reference per non-BROKEN triple of ``regrid_grib_cases``, and without
    it the download cache may hold any subset -- but no file the table does not expect.
    """
    expected = _expected_grib_npz_files()
    override = os.environ.get(REGRID_TEST_DATA_ENV)
    root = Path(override) if override else CACHE_ROOT
    folder = root / "regrid" / "grib"
    actual = {path.relative_to(root).as_posix() for path in folder.rglob("*.npz")} if folder.is_dir() else set()

    if override:
        assert actual == expected, (
            f"{len(actual)} npz files under {folder} ({REGRID_TEST_DATA_ENV}), the GRIB table expects "
            f"{len(expected)}; missing: {sorted(expected - actual)[:5]}, extra: {sorted(actual - expected)[:5]}"
        )
    else:
        assert actual <= expected, (
            f"{len(actual - expected)} GRIB npz files under the download cache {folder} are not files of the "
            f"table (the cache is fetched on demand and may be deleted): {sorted(actual - expected)[:5]}"
        )


@requires_mir
@requires_ekd
@pytest.mark.long_test
@pytest.mark.timeout(300, method="thread")
@pytest.mark.parametrize("key", _matrix_params())
def test_regrid_fieldlist_grib_matrix(key):
    """Regrid one whole GRIB fieldlist onto one target and check it against the table and MIR.

    Parameters
    ----------
    key : tuple
        The (canonical interpolation, filename, target id) triple of the case, from
        ``regrid_grib_cases.cases()``. A BROKEN triple brings its xfail mark with it: a
        "raises" case runs and must fail with the exception the table names (``strict``, so a
        MIR that stops failing here is a failure, not a silence), a "hang" case is never
        executed. A KnownBug triple brings ``xfail(strict=True)``: it runs and is checked
        against the CORRECT result its entry records (regrid-grib-findings §C); its stored npz
        is the stale no-op output of the bug, so the reference comparison is skipped for it,
        while the property checks, the reference-free nearest rule and the metadata check run.
    """
    interpolation, filename, target_id = key
    entry = EXPECTED[key]
    target = _TARGET_BY_ID[target_id]

    ds = _grib_fieldlist(filename)

    # a "raises" case ends here, with the exception the table records. Nothing of the input is
    # read before this point on purpose: a file whose own grid metadata contradicts itself (§A,
    # the two INVALID_GRIB_FILES of the negative test below) fails while its values are
    # unpacked, with a different exception than the one regrid raises, so reading the input
    # first would hide that finding.
    out = regrid(ds, out_grid=target.spec, interpolation=interpolation)

    in_fields = list(ds)
    in_values = [np.asarray(field.to_numpy()) for field in in_fields]
    metadata_ref = ds.metadata(METADATA_KEYS)

    if not isinstance(entry, (Expected, KnownBug)):
        raise AssertionError(
            f"{key}: the table marks this BROKEN ({entry.reason}), but MIR produced a result. "
            f"Reprobe it with generate_regrid_grib_data.py --file {filename}"
        )

    assert len(out) == entry.n_fields, f"{key}: {len(out)} result fields, the table has {entry.n_fields}"
    assert len(in_fields) == entry.n_fields, (
        f"{key}: the input has {len(in_fields)} fields, the result {entry.n_fields}"
    )

    # a KnownBug triple is checked against the CORRECT result its entry records; its stored npz
    # is the stale no-op output of the §C MIR bug, and comparing values to it would keep the
    # case failing after a MIR fix, killing the XPASS alert. Expected entries have every field
    # below, so all the is-not-None guards only ever take effect for a KnownBug entry.
    known_bug = isinstance(entry, KnownBug)
    reference = None
    if not known_bug:
        reference = load_values(data_path("grib_reference", interpolation, filename, target_id))
        assert reference.ndim == entry.shapes[0].__len__() + 1, (
            f"{key}: reference {reference.shape} stacks {entry.shapes}"
        )
    compare = get_interpolation(interpolation).compare

    for index, field in enumerate(out):
        where = f"{key} field {index}"
        values = np.asarray(field.to_numpy())
        shape = tuple(values.shape)

        if entry.shapes is not None:
            assert shape == tuple(entry.shapes[index]), f"{where}: shape {shape}, the table has {entry.shapes[index]}"
        if entry.grid_types is not None:
            assert field.metadata("gridType") == entry.grid_types[index], (
                f"{where}: gridType {field.metadata('gridType')!r}, the table has {entry.grid_types[index]!r}"
            )
        spec = field.geography.grid_spec()
        if entry.grid_specs is not None:
            _assert_gridspec(spec, entry.grid_specs[index], f"{where} gridspec")
        else:
            # the KnownBug triples with an unknown post-fix spec (a sub-area clip nothing can
            # observe while the bug hides it) assert only that the result is unrotated
            still_rotated = isinstance(spec, dict) and "projection" in spec
            assert not still_rotated, (
                f"{where}: gridspec {spec!r} keeps the rotation of the input, so the result is not the"
                " unrotated grid the target asks for (regrid-grib-findings §C)"
            )
        if entry.nan_counts is not None:
            nan = int(np.isnan(values).sum())
            assert nan == entry.nan_counts[index], f"{where}: {nan} NaN values, the table has {entry.nan_counts[index]}"

        if known_bug:
            # no reference for the values; the reference-free nearest rule stays -- the correct
            # result must hold input values just as much as the buggy one
            if interpolation == "nearest-neighbour" and in_fields[index].metadata("gridType") != "sh":
                _assert_values_from_input(values, in_values[index], where, _packing_allowance(field, in_fields[index]))
            continue

        ref = reference[index]
        assert ref.shape == shape, f"{where}: the reference has {ref.shape}, the result {shape}"

        if interpolation == "nearest-neighbour":
            try:
                compare(values, ref)
            except AssertionError:
                # the only accepted reason for a difference is a tie: a source point as near to
                # this target point as the nearest one. Both grids are the points themselves: the
                # tie is a question about points and the coordinate mode asks nothing else, so it
                # covers every grid of the matrix, the projected inputs among them.
                assert_nearest_ties(
                    values,
                    ref,
                    in_values[index],
                    _field_points(in_fields[index], f"{where} input"),
                    _field_points(field, where),
                    context=where,
                )
            # independent of the reference: a nearest result only moves values around, and the
            # only deviation allowed for it is the one the two GRIB packings can explain -- not
            # checked for a spectral input (GRIB gridType "sh"): the field's array holds
            # spherical-harmonic coefficients, not grid values, so the result values cannot be
            # members of it
            if in_fields[index].metadata("gridType") != "sh":
                _assert_values_from_input(values, in_values[index], where, _packing_allowance(field, in_fields[index]))
        else:
            # the pole rows of a global regular lat-lon result are machine dependent, and only
            # a result the table records as the whole unclipped 10x10 grid may skip them
            skip = target_id == "10x10" and shape == _GLOBAL_10X10_SHAPE and _is_global_10x10(entry.grid_specs[index])
            compare(values, ref, skip_pole_rows=skip)

    assert out.metadata(METADATA_KEYS) == metadata_ref, f"{key}: the regridding changed param, level, date or time"


# ---------------------------------------------------------------------------
# the negative set: the files that are not in the matrix because they cannot regrid correctly
# ---------------------------------------------------------------------------


@requires_mir
@requires_ekd
@pytest.mark.long_test
@pytest.mark.timeout(300, method="thread")
@pytest.mark.parametrize("interpolation, filename, target_id", list(invalid_cases()))
def test_regrid_fieldlist_grib_invalid_file(interpolation, filename, target_id):
    """Regrid one file of ``INVALID_GRIB_FILES`` and assert the error ``INVALID`` records.

    The two files contradict their own grid metadata (regrid-grib-findings §A, responsible: the
    data files): ``regular_gg_wrong_increment.grib`` holds Ni*Nj = 192*64 = 12288 for 18432
    data points, which ecCodes refuses to interpret, and ``regular_ll_wrong_increment.grib``
    encodes 2 degree increments on a 5 degree grid, which MIR's own GRIB reader asserts on.
    Neither is a result MIR could be compared against, so they are not matrix triples and have
    no reference npz: what their regridding must do is fail, and this test says how, with
    ``pytest.raises`` on the exception type and message of its :data:`regrid_grib_cases.INVALID`
    entry. The three combinations of the matrix run on each of them, so a reader sees that the
    refusal comes from the grid, not from the interpolation or the target.

    Parameters
    ----------
    interpolation : str
        The canonical name of the interpolation the case regrids with, from
        ``regrid_grib_cases.invalid_cases()``.
    filename : str
        The GRIB file, one of ``regrid_grib_cases.INVALID_GRIB_FILES``.
    target_id : str
        The id of the ``regrid_grib_cases.GRIB_TARGETS`` entry the case regrids onto.
    """
    key = (interpolation, filename, target_id)
    entry = INVALID[key]
    target = _TARGET_BY_ID[target_id]
    # resolved here and not at collection: gribapi comes with earthkit-data, and invalid_cases()
    # resolves nothing, which is what keeps this module importable, and collectable, without it
    exception_class = _exception_class(entry.exception)

    ds = _grib_fieldlist(filename)

    # as in the matrix test, nothing of the input is read before the regridding: unpacking the
    # values of regular_gg_wrong_increment.grib fails with a ValueError about the reshape, an
    # error of reading, not of regridding, and not the one the table records for this triple
    with pytest.raises(exception_class, match=entry.match):
        regrid(ds, out_grid=target.spec, interpolation=interpolation)
