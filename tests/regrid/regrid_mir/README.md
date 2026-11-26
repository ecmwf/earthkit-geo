# Regridding tests (MIR)

These tests check the regridding of `earthkit.geo` against **MIR**, ECMWF's interpolation
library. MIR is an optional dependency: without it these tests are skipped with the reason "No
mir available"; a few of them also need earthkit-data and skip with "No earthkit-data
available".

Words used below:

- **source grid**: the grid the values live on before regridding.
- **target grid**: the grid the values are moved onto.
- **gridspec**: the dictionary that describes a grid, e.g. `{"grid": [10, 10]}` (global
  regular lat-lon, 10 degree spacing) or `{"grid": "O32"}` (a reduced Gaussian grid).
- **interpolation**: how the values are moved. Three are registered: `nearest-neighbour`
  (aliases `nn`, `nearest-neighbor`), `linear` and `grid-box-average`.
- **reference**: a regridding result that MIR produced once and that was saved to a npz file
  (npz = a zip of numpy arrays, here always with the single key `arr_0`). Tests compare a new
  result against it. See [The data files](#the-data-files).

## What is in the folder

Thirteen python files: six modules that are not tests (the two shared helper modules, the grid
table, the GRIB table and the two data generators), one table-driven test module for the numpy
matrix, and six modules of older test cases — one of which, `test_regrid_fieldlist_grib.py`, now
also holds the table-driven GRIB matrix.

- `regrid_common.py` — not a test: the interpolation registry (canonical name, aliases and how
  a result is compared with its reference, plus `assert_nearest_ties()`, the reference-free
  check of a nearest-neighbour tie), the `interpolation=` kwarg spellings, the `requires_mir` /
  `requires_ekd` skip marks, and the one data access function `data_path()`.
- `regrid_fieldlist_common.py` — not a test: what the two fieldlist test modules share,
  the fieldlist loader, the 10x10 reference check and the `..._to_ll` case table with its
  test body, so that each of them is only the GRIB or the array half of the same use
  cases. It has no side effects at import and imports without earthkit-data installed.
- `regrid_cases.py` — not a test: the grid table (sources, targets, per-interpolation
  selection, `MAX_POINTS`) plus the `BROKEN` and `CLIPPED` findings of the probe.
- `regrid_grib_cases.py` — not a test: the GRIB matrix table (`GRIB_FILES`, `GRIB_TARGETS`,
  `GRIB_COMBINATIONS`) and, per (interpolation, file, target) triple, what MIR produced
  (`Expected`), why it cannot (`Broken`), or what it *should* produce where a verified MIR bug
  answers wrong (`KnownBug`, §C), plus the negative set (`INVALID_GRIB_FILES`, `INVALID`,
  `invalid_cases()`) — the 2 server files whose own grid metadata contradicts itself, kept out
  of the matrix and tested for the errors their regridding raises. See
  [The GRIB matrix](#the-grib-matrix).
- `generate_regrid_data.py` — not a test: the script that writes the numpy matrix's npz files,
  by hand only, never during a test run.
- `generate_regrid_grib_data.py` — not a test: the same for the GRIB matrix: it regrids every
  matrix GRIB file onto every target, writes the references and prints the table entries (it
  iterates `GRIB_FILES`, so it never probes the files of the negative set).
- `test_regrid_numpy_grids.py` — the numpy matrix: every covered source of the table onto every
  target, for the interpolations of its selection. Marked `long_test`.
- `test_regrid_numpy.py` — the older numpy array cases: a handful of grids onto a 10x10
  global lat-lon with reference files, the `interpolation=` kwarg spellings, shape-only
  any-to-any and any-to-points cases, gridspec spellings.
- `test_mir_grid.py` — smoke test: any gridspec MIR accepts regrids onto any other, checked by
  shape only, no data files.
- `test_regrid_xarray.py` — regridding `xarray.DataArray` / `xarray.Dataset` objects.
- `test_regrid_fieldlist_grib.py` — regridding earthkit-data GRIB fieldlists and single
  fields, plus the GRIB round trips and the ecCodes grid-spec check (30 tests that always run,
  one smoke test, two fast spectral regressions and the fast array-path regression among them;
  four of the 30 — the table, the two comparison rules and the packing allowance — need no
  earthkit-data at all, so they still check the matrix where the GRIB stack cannot be
  imported), the GRIB matrix: 120 `long_test` cases over `regrid_grib_cases.cases()` plus the
  files check of that matrix, and the negative set: 6 `long_test` cases over
  `regrid_grib_cases.invalid_cases()`.
- `test_regrid_fieldlist_array.py` — the same inputs as non-GRIB (array) fieldlists and
  single fields (16 tests).
- `test_regrid_grib_message.py` — regridding a GRIB message (a buffer, not a fieldlist).

Three families:

- **numpy arrays**, through `earthkit.geo.grids.array.regrid`: the table-driven matrix of
  `test_regrid_numpy_grids.py` plus the older use cases of `test_regrid_numpy.py` and
  `test_mir_grid.py`. This is the family that got the new structure.
- **xarray**, through `earthkit.geo.regrid`: `test_regrid_xarray.py`.
- **fieldlists and GRIB messages**, through `earthkit.geo.regrid`:
  `test_regrid_fieldlist_grib.py`, `test_regrid_fieldlist_array.py` and
  `test_regrid_grib_message.py`.

The xarray, fieldlist and GRIB message files are **refactored only**: same use cases, same
data files, same assertions, with the duplicated lists, marks and test bodies replaced by
helpers from `regrid_common.py`. The fieldlist file is also **split** in two, the GRIB
half and the non-GRIB (array) half, with what both need in `regrid_fieldlist_common.py`;
its 36 old tests are the same 36 (20 + 16). The GRIB half has since grown the matrix of
[The GRIB matrix](#the-grib-matrix), which is why the folder totals below are larger than the
ones of 2026-11-19.

## Running the tests

Activate an environment with MIR installed (for example `source ~/venv/dev/bin/activate`) and
work from the repository root. `pytest.ini` turns on coverage for every run, so add `--no-cov`
for a plain one.

`tests/conftest.py` adds the `-E` option: the default `-E short` skips every test marked
`long_test` — as in `tests/grid/`, and in this folder that is the two matrices, their files
checks and the negative GRIB set, and nothing else — while `-E long` runs it. The counts in
the comments are from 2026-11-20, with MIR and
earthkit-data installed, and are the **server/cache** numbers: the 115 GRIB matrix references
are not on the download server yet, so a plain long run has 109 failures (the six `KnownBug`
cases of the §C MIR bug xfail on their own checks without ever reading a reference), and the
same run against the locally generated tree (`EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data`,
see [The GRIB matrix](#the-grib-matrix)) is green:

```
python -m pytest tests/regrid/regrid_mir --no-cov -q          # short: 204 passed, 548 skipped, 8 xfailed
python -m pytest tests/regrid/regrid_mir -E long --no-cov -q  # everything: 628 passed, 109 failed, 5 skipped, 18 xfailed
EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
    tests/regrid/regrid_mir -E long --no-cov -q               # the same, on the local tree: 737 passed, 5 skipped, 18 xfailed
```

The 18 xfailed are the 7 BROKEN pairs of the numpy matrix plus the 5 BROKEN and the 6
`KnownBug` (strict-xfail) triples of the GRIB matrix; the 6 negative cases pass — they assert
their errors with `pytest.raises`. The short run counts the one GRIB triple that is never
executed (`xfail(run=False)`) among its xfailed rather than its skipped; the six
strict `xfail(strict=True)` marks of the `KnownBug` triples only speak in a long run, where
their cases really execute.

The long run fetches the npz files it needs from the download server into the cache
`tests/data/regrid/` (499 files, 74 MB for the numpy matrix, 115 more for the GRIB one — the
first run is slow, every later one reads from the cache, and deleting the folder is safe). The 5
skipped in the long run have nothing to do with `-E`: they are older use cases that skip for
their own reason, 3 in `test_regrid_numpy.py` (a reference that needs updating, one case marked
too slow) and 2 in `test_regrid_xarray.py` (currently failing). Without `-E long` the numpy
matrix module runs nothing at all:

```
python -m pytest tests/regrid/regrid_mir/test_regrid_numpy_grids.py --no-cov -q  # 417 skipped, 7 xfailed
```

The matrix has one test per (interpolation, source, target) pair, with the id
`<interpolation>-<source id>-<target id>`; select with `-k`:

```
python -m pytest tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -q -k "nearest-neighbour-O32-ll-10x10"
# 1 passed, 420 deselected

python -m pytest tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -q -k O32
# 121 passed, 2 xfailed
```

### Combinations MIR cannot do (BROKEN)

`regrid_cases.BROKEN` holds one entry per (interpolation, source, target) triple that MIR
cannot produce — it raises, aborts or hangs — with the MIR error as its reason. Each becomes a
test parameter marked `pytest.mark.xfail(run=False, reason="BROKEN: ...")`, which means:
**never executed**, listed in the report as `xfailed` with that reason. The run stays green
while MIR still fails, and the failure is on the list instead of being silently skipped. There
are 7 of them, the "7 xfailed" of the reports above (`-k O32` above shows the two that involve
O32).

**Because the test never runs, pytest can never notice that one of them works again**: an
`xfail(run=False)` parameter cannot `xpass`, and the generator writes no reference file for it
either, so a MIR upgrade that fixes a BROKEN combination stays invisible. Reprobe them by hand
after an upgrade — `--include-broken` makes the generator run those pairs too and report what
still fails, which is the only way a BROKEN combination can say anything. One command per
combination, from the repository root:

```
python tests/regrid/regrid_mir/generate_regrid_data.py --include-broken --force \
    --interpolation nearest-neighbour --source swisslv95 --target ll-10x10
```

Run it for each of the 7 entries and look at the status line: `ok` means MIR produced a result,
so the reference file is there — delete that entry from `BROKEN`, and the pair comes back as a
normal test parameter on the next run. A `failed`, `aborted` or `timeout` line quoting the same
error means the entry stays. Two notes on the command: `--force` is only needed when an earlier
reprobe has already written the reference (a BROKEN pair has no file, so it would run without
it), and it then also rewrites the **input** file of that source — the same deterministic field
from the same coordinates, so the references computed from it stay valid (the warning about
`--force` on inputs, in [Regenerating](#regenerating), is about changing `smooth_field`).
Appending `--dry-run` to the command prints what would run and writes nothing; so does
`--include-broken --dry-run` without filters, which then lists all 412 pairs, the 7 BROKEN ones
among them.

The GRIB matrix marks its BROKEN triples the same way, with one difference worth knowing: a
GRIB triple that *raises* is marked `xfail(raises=<exception>, strict=True)` and really runs, so
a MIR or earthkit-data fix that makes it work turns into an `XPASS` and a **failure** that says
so, instead of staying invisible. Only the one triple that is too slow to run
(`mercator.grib` linear onto `10x10`) is `xfail(run=False)`. The six `KnownBug` triples of the
§C MIR bug are marked in the same strict spirit — `xfail(strict=True)`, no `raises=`, and they
really run and fail against the *correct* result their table entry records — see
[The GRIB matrix](#the-grib-matrix).

## The grid table (`regrid_cases.py`)

**Sources.** All 149 grids of the grid tests (`tests/grid/test_grid_types.py`,
`test_grid_named.py`, `test_grid_sh.py`, `test_grid.py`), copied from those files as gridspec
dicts so that a change in a grid test never silently changes a regrid test, and de-duplicated
only where two spellings are the same grid in the same order. 94 of them are covered, i.e.
they are regrid sources here; the rest stay in the table with a reason. A source id becomes a
file name, so it uses letters, digits, `-`, `_` and `.` only.

**Targets.** The four grids sources are regridded onto:

| id            | gridspec                                      | result shape | what it is                         |
| ------------- | --------------------------------------------- | ------------ | ---------------------------------- |
| `ll-10x10`    | `{"grid": [10, 10]}`                          | 19 x 36      | global regular lat-lon, coarse     |
| `O32`         | `{"grid": "O32"}`                             | 5248         | global octahedral reduced Gaussian |
| `H8-ring`     | `{"grid": "H8", "order": "ring"}`             | 768          | HEALPix, Nside 8, no rows          |
| `ll-1x1-area` | `{"grid": [1, 1], "area": [72, -25, 25, 45]}` | 48 x 71      | regional 1 degree grid over Europe |

**`MAX_POINTS = 1_000_000`.** A source of a million points or more gets no data files: its
coordinates would be gigabytes, and regridding them can take minutes. It stays in the table,
flagged `too large: <n> points`, listed but not covered; raising the constant is the only
change needed to cover it. 43 of the 149 sources are above the gate (the smallest has
1,036,800 points).

**`NOT_COVERED`.** Derived from the table: source id to reason. Apart from the 43 too-large
grids it holds 8 `sh` grids (a spherical-harmonic truncation has no grid points), 3 Cartesian
`regular_xy` grids without a projection (no lat/lon at all), and `HN10` (whose `to_latlons()`
raises: "Nside must be a power of 2"). Not being a source is a statement about the coordinates,
not about the grid: `check_table()` builds all of them and compares shape and size, so a spec
that does not build is reported wherever it is in the table.

**`SELECTION`.** What runs against what, per interpolation:

- `nearest-neighbour`: every covered source against every target — 94 x 4 = 376 pairs.
- `linear` and `grid-box-average`: a hand-picked subset each, the same 9 small sources
  (`ll-1x1`, `ll-area-string`, `F32`, `O32`, `N32`, `H8`, `H8-nested`, `ORCA2_T`, `pi_C` —
  global and regional lat-lon, regular and reduced Gaussian, HEALPix in both orders, ORCA and
  FESOM) against 2 global targets (`ll-10x10`, `O32`) — 18 pairs each.

That is 412 pairs; 7 are BROKEN, so 405 have data files. `test_regrid_numpy_grids.py` collects
424 tests: one per pair, four that check the table itself — that every spec builds with the
shape and size it claims (`check_table()`, which builds **all** 149 sources including the
non-covered ones, plus the 4 targets), that every non-covered source has a reason and that
BROKEN and CLIPPED name real pairs of the matrix, that every BROKEN reason starts with
`BROKEN: `, and that the npz files on disk match the table (the exact set for a local tree
named by `EARTHKIT_GEO_REGRID_TEST_DATA`, no unexpected file in the download cache) — three
that check how `data_path` finds and fetches the files (missing override file, failing
download, fetched-once cache; all three fake the download function, no network), and five
that check the comparison rules of `regrid_common.py` on synthetic arrays (see
[How a result is compared](#how-a-result-is-compared)).

## Adding a grid, a target or an interpolation

### A source grid

1. Copy its gridspec from the grid test that uses it and add one `_source(...)` line to the
   type section it belongs to in `regrid_cases.py`, with the id, the spec, the shape (from
   `Grid(spec).shape`) and, for a grid that declares an area, `region="regional"`. The size
   and the coverage come from the shape, so a large grid needs no extra work. The module
   docstring lists the details (naming, spellings to merge, a grid that cannot be a source at
   all needs an explicit `reason=`).

1. Check the table: it prints `[]` when the table is right. It builds every grid, which is why
   nothing like it runs at import; run it from `tests/regrid/regrid_mir/`, since the folder is
   not a package:

   ```
   python -c "import regrid_cases as c; print(c.check_table())"
   ```

1. From the repository root, see what the new grid would run, then write its files. `--source`
   takes ids, and nearest-neighbour picks up a new covered source by itself:

   ```
   python tests/regrid/regrid_mir/generate_regrid_data.py --source <new-id> --dry-run
   python tests/regrid/regrid_mir/generate_regrid_data.py --source <new-id>
   ```

1. Run the matrix on it — on the local tree, since the new files are not on the server
   until they are uploaded ([The data files](#the-data-files)):

   ```
   EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
       tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -k <new-id>
   ```

1. When the tests pass, upload the `regrid/` tree (or just the new files) to the download
   server — the maintainer's step — after which plain runs fetch the new data too.

For `linear` or `grid-box-average` to cover it as well, add the id to `SMALL_SOURCES` (or to
the selection of that interpolation) in `regrid_cases.py`, then generate with
`--source <new-id> --interpolation linear` and go on as in steps 3 and 4. Two worked examples of
filtered runs, which print the pairs without writing anything:

```
python tests/regrid/regrid_mir/generate_regrid_data.py --source O32 --dry-run
python tests/regrid/regrid_mir/generate_regrid_data.py --source pi_C --target O32 --interpolation linear --dry-run
```

### A target

Append `Target(id, spec, shape)` to `TARGETS` in `regrid_cases.py`. `nearest-neighbour` takes
it for every covered source; the other two selections name their targets and need editing.
Then check the table (step 2 above), generate, run against the local tree, and upload the
new files once it passes (as for a source grid, steps 3 to 5):

```
python tests/regrid/regrid_mir/generate_regrid_data.py --target O32 --dry-run
EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
    tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -k O32
```

### An interpolation

1. In `regrid_common.py`: write a `_compare_<name>(result, reference, skip_pole_rows=False, **kwargs)` function that raises `AssertionError` on mismatch, and register
   `Interpolation("<name>", (<aliases, or ()>), _compare_<name>)` in `INTERPOLATIONS`. The
   canonical name is also the folder of its reference files. The matrix test passes
   `skip_pole_rows=True` for a pair that really returns the whole global `ll-10x10` grid, so an
   `assert_allclose` rule takes that keyword; a strict rule that copies values, like
   nearest-neighbour's, refuses it and leaves the machine-dependent differences to the tie
   validation of the next section.

1. In `regrid_cases.py`: add one key to `SELECTION` with the source ids and target ids it
   should cover — every covered source and every target for an interpolation that must run on
   all of them, a subset otherwise. Nothing else changes: `test_regrid_numpy_grids.py` is
   parametrised from `SELECTION`, and `BROKEN` / `CLIPPED` are keyed by that same name.

1. Generate its references, run them against the local tree (the override — the new files
   are not on the server yet), and upload the `regrid/` tree once the tests pass:

   ```
   python tests/regrid/regrid_mir/generate_regrid_data.py --interpolation <new-name> --dry-run
   python tests/regrid/regrid_mir/generate_regrid_data.py --interpolation <new-name>
   EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
       tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -k <new-name>
   ```

   `--interpolation` takes aliases as well as canonical names, so this prints the 376
   nearest-neighbour pairs that a run for `nn` would (re-)generate:

   ```
   python tests/regrid/regrid_mir/generate_regrid_data.py --interpolation nn --dry-run
   ```

## The data files

One npz file per array, always with the single key `arr_0`, in a tree whose layout is the
same at every root:

```
<root>/regrid/input/<source>.npz                        the input values of a source grid
<root>/regrid/<interpolation>/<source>/<target>.npz     MIR's result: the reference
<root>/regrid/grib/<interpolation>/<grib file>/<target>.npz   the same, for the GRIB matrix
```

The tree (94 inputs + 405 references, 499 files, 74 MB) lives on the ECMWF download server
under `https://sites.ecmwf.int/repository/earthkit-geo/test-data/regrid/` —
`URL_ROOT + "/regrid/"` of `earthkit.geo.utils.testing`. `regrid_common.data_path` resolves
new-hierarchy files through `_fetch_new()`, the one place that knows the roots: a file is
read from the local cache `tests/data/regrid/` (git-ignored — deleting it is always safe,
the next run fetches what it needs again) and, when it is not cached yet, downloaded from
the server into a temporary file in the same folder and renamed atomically, so a broken
download never leaves a half file behind. A failed download **fails** the test with an
error that names the URL, with a hint at the generator; it is never silently skipped.
Uploading a regenerated tree to the server is the maintainer's manual step — the generator
only writes locally.

Data that is not on the server yet — freshly generated files — is tested through the
override: set the environment variable `EARTHKIT_GEO_REGRID_TEST_DATA` to a folder that
contains the `regrid/` tree — the generator's output folder, `_dev/test-data` by default
or wherever its `--output-dir` (or the `EARTHKIT_GEO_REGRID_GENERATE_DIR` it defaults to)
pointed. (The generator sets the variable for its own process and workers anyway, so
generating never touches the network.) `data_path` then reads the files of the **new
hierarchy** — the `input`, `reference` and `grib_reference` kinds — from that folder and never
downloads them; the `legacy` kind is not affected, so the GRIB inputs of the fieldlist tests
and the older shared npz files keep coming from the server cache (see [Other data
files](#other-data-files)). A missing override file **fails** with a message that quotes the
generator command.
The tests check the two modes separately: an override tree must hold exactly the files the
table expects; the cache must hold no file the table does not expect (it need not hold all
of the expected ones — they are fetched on demand). Each matrix checks its own part of the
tree: `test_regrid_numpy_grids_files` looks at `regrid/` and ignores `regrid/grib/`, and
`test_regrid_fieldlist_grib_files` looks at `regrid/grib/` alone.

That is also why the run reports the generators write sit outside `regrid/`, at
`<output root>/regrid-run-report.json` and `<output root>/regrid-grib-expected.json`
(`_dev/test-data/...` with the default output folder): the tree is uploaded as it is, no
sidecar files.

The inputs are not noise: every source grid gets the same deterministic smooth function of its
latitude and longitude (`smooth_field` in `generate_regrid_data.py`), so a regridding has a
known shape and a nearest-neighbour result can be checked without a reference at all. The
**references are a regression record made by MIR itself, not an independent truth**: a green
matrix says that MIR still answers the way it did when the files were written, not that the
answer is correct. Two families of checks do not use the reference: the result shape and the
returned gridspec, and, for `nearest-neighbour`, that every non-NaN output value is exactly one
of the input values (nearest only moves values around, it never invents any).

### How a result is compared

The rules live in `INTERPOLATIONS` in `regrid_common.py`; `test_regrid_numpy_grids.py` runs the
rule of the interpolation for each pair, and keeps the two reference-free checks of the previous
section (shape and returned gridspec, and the values of a nearest result being input values).

- **`nearest-neighbour` is strict**: `NaN` at exactly the same positions in result and
  reference, and every finite value **equal** (`rtol=0`, `atol=0`; the tolerance is a kwarg of
  the rule, but nothing passes it). No pole-row exemption, no share of values allowed to
  differ.
- **Ties are handled where they can be checked, not excused.** MIR looks its nearest neighbour
  up with a search, and where several source points have the *same* distance to a target point
  another build of MIR may answer with another one of their values — systematically, in every
  value at once. `ll-reference` onto `ll-1x1-area` is the plain case: all 3,408 target points
  sit exactly between four source cell centres, so 3,405 of the values may legitimately differ
  elsewhere. A "99% of the values must match" rule cannot tell that from a broken result, so
  when the strict comparison fails the test runs `regrid_common.assert_nearest_ties` on the
  positions that differ: each of them must hold the input value of a source point whose
  distance to the target point is within `1e-9` relative and `1e-12` absolute tolerance (3D
  chords on the unit sphere, looked up in a `scipy` KD-tree) of the distance of the nearest
  source point; anything else fails. That is an independent check, and it only runs on failure,
  so a passing pair costs nothing. A value from a *farther* source point fails, which the
  "every output value occurs in the input" check cannot catch on its own.
- **`linear` and `grid-box-average`** compare the whole array with `assert_allclose` (default
  tolerance, and its `equal_nan=True` already requires the NaN positions to match). The two
  rows of a duplicate pole are machine dependent in MIR, so both rules take a
  `skip_pole_rows=` keyword: the matrix test sets it only where the result really is the whole
  global regular lat-lon target — target `ll-10x10` and the pair not listed in `CLIPPED` — and
  then rows 0 and -1 are left out, nothing else ever is. Deriving the exemption from the shape
  instead, as an earlier version did, left every regional 2D result unchecked: a clipped
  one-row result *is* its pole rows, so `xy-laea-grid-bbox-reverse` onto `ll-10x10` (1 x 1)
  accepted any value at all, and replacing the reference value with the input value of a source
  point on the other side of the map passed.
- **A `BROKEN` pair is never compared**: it is an `xfail(run=False)` parameter, see
  [Combinations MIR cannot do (BROKEN)](#combinations-mir-cannot-do-broken).

### Regenerating

```
python tests/regrid/regrid_mir/generate_regrid_data.py --help
```

The generator writes its tree under an **output root** — the folder that will contain
`regrid/` — and never uploads or downloads anything. Without options that is
`_dev/test-data`, the default output folder; `--output-dir PATH` picks another one
(relative paths count from the current working directory, the folder is created if
missing), and the environment variable `EARTHKIT_GEO_REGRID_GENERATE_DIR` is the default
when the option is not given (the option wins over the variable, the variable over
`_dev/test-data`). The generator sets `EARTHKIT_GEO_REGRID_TEST_DATA` to the output root
for itself and its workers, so its own reads stay local whatever the root is. After a
run, test the tree through that override — with `EARTHKIT_GEO_REGRID_TEST_DATA` set to
**that folder** — and, once it is green, upload that folder's `regrid/` tree to the
download server — the maintainer's step, after which plain runs fetch the new files into
`tests/data/regrid/`:

```
EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest tests/regrid/regrid_mir -E long --no-cov -q

python tests/regrid/regrid_mir/generate_regrid_data.py --output-dir /tmp/regrid-data --source O32
EARTHKIT_GEO_REGRID_TEST_DATA=/tmp/regrid-data python -m pytest \
    tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -k O32
```

The whole table (94 inputs + 376 + 18 + 18 pairs), and inputs only, four workers, overwriting
what exists. Both here with `--dry-run`, which prints one line per file and writes nothing;
drop it to actually write. (The probe of the whole matrix took about 4 minutes of wall clock
with `--jobs 3`.)

```
python tests/regrid/regrid_mir/generate_regrid_data.py --dry-run
python tests/regrid/regrid_mir/generate_regrid_data.py --inputs-only --jobs 4 --force --dry-run
```

An existing file is never overwritten without `--force`, and the filters of the sections above
(`--source`, `--target` and `--interpolation`, all three repeatable) keep a run cheap.
`--output-dir PATH` chooses the output root the run writes to (default
`EARTHKIT_GEO_REGRID_GENERATE_DIR`, or `_dev/test-data`), `--jobs N` runs N pair
subprocesses at a time, `--timeout SECONDS` (default 120) is per pair, and
`--include-broken` also runs the pairs listed in BROKEN — which is how the probe covered
all 412 of them. Every pair is regridded in its own subprocess because MIR can abort the
interpreter or hang on some grid combinations. The exit code is 0 even when pairs fail: a
failure is a finding for `BROKEN`, and the run prints them.

Regenerate when a new MIR version changes results: run it with `--force` over what it touched
(or everything), re-run the tests and look at what fails. Each failure is either a change to
accept — the reference is then simply written again — or a MIR bug to report. Beware that
`--force` on the inputs invalidates the references computed from them: regenerate those too,
with the same command without `--inputs-only` (and with filters to keep it cheap).

## The GRIB matrix

The numpy matrix regrids arrays the tests build themselves. The **GRIB matrix** regrids the real
GRIB files of the download server, the way a user reads them: `regrid_grib_cases.GRIB_FILES` is
**40 of the 42 files** of the server's `test-data/grib/` folder — the two whose grid metadata
contradicts itself are the negative set below —, each one read whole with
`earthkit.data.from_source("file", ...)` and passed to `earthkit.geo.regrid` with a gridspec. The
targets are `10x10` (`{"grid": [10, 10]}`, a 19 x 36 result) and `O32` (`{"grid": "O32"}`, 5248
values); `GRIB_COMBINATIONS` runs three of the six interpolation × target possibilities —
nearest-neighbour and linear onto `10x10`, nearest-neighbour onto `O32`. That is **120 triples**,
one `long_test` case each in `test_regrid_fieldlist_grib.py`, with the id
`<interpolation>-<file>-<target>`. The module's other 30 tests are the older fieldlist use cases,
the checks of the table and of the reference-free rules, one smoke test and three fast
regressions (the two spectral ones and the §E array-path one); the matrix's 120 cases, the 6 of
the negative set and the matrix's files check are the other 127 of its 157. Four of the 30
read nothing but the table and numpy — `test_regrid_fieldlist_grib_table`,
`..._nearest_ties_in_coordinates`, `..._nearest_values_are_input_values` and
`..._nearest_packing_allowance` — so they still run on a machine without earthkit-data, where the
120 matrix cases collapse onto one placeholder parameter, the 6 negative cases still collect
(`invalid_cases()` resolves nothing), and the rest of the module skips (4 passed, 34 skipped).

```
EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
    tests/regrid/regrid_mir/test_regrid_fieldlist_grib.py -E long --no-cov -q -k regular_ll.grib1
# 3 passed, 154 deselected
```

**What one case checks**, against its entry in `regrid_grib_cases.EXPECTED` — and, for a triple
the table records as `Expected`, against MIR's reference file; a `KnownBug` triple is instead
checked against the **correct** result its entry records, with no reference comparison (its
stored npz is the stale output of the §C bug): the number of result fields, and per
field the value shape, the `gridType`, the gridspec `geography.grid_spec()` returns and the NaN
count; then the values, with the comparison rule of the interpolation
([How a result is compared](#how-a-result-is-compared)); then that `param`, `level`, `date` and
`time` are still the input's. For a `KnownBug` triple the property checks are only those the
entry spells out — the two sub-area ones record n_fields and an unrotated result, because the
post-fix clip is unknown while the bug hides it — but the reference-free nearest rule below
still applies to its nearest-neighbour cases. For nearest-neighbour it also checks, with no
reference involved,
that every non-NaN value of a field **is** a value of that input field — except for a spectral
input, whose array holds spherical-harmonic coefficients, not grid values, so the rule is skipped
for the six spectral triples (the spectral-input bug is fixed: those six are `Expected` now) —
exact membership, the only allowance being what the two GRIB messages can account for — half the
value step of the coarser of them (`2**binaryScaleFactor * 10**-decimalScaleFactor`, the lattice simple
and complex packing, also with spatial differencing, round onto; the lossy compressions report no
such thing, so they give no allowance) — plus one float32 ulp of the value. Of the 178 result fields
of the 73 non-spectral nearest triples, 168 need neither. The ten that do: the two `mercator.grib`
fields are copies re-packed onto
the complex lattice of the result, 0.4 of its step away (`binaryScaleFactor` -7 onto `10x10` and -5
onto `O32`; the input's `4.3` comes back as `4.30312502`, its `8.2` as `8.1875`); the two
`reduced_rotated_gg_subarea_*` ones deviate by 7.6e-06 and 5.7e-06, far inside the 0.015625 their
simple step of 0.03125 allows; the two `reduced_gg_*.grib2` ones by the same 7.3e-12 — less
than half an ulp, which is all they get, their `grid_ccsds` packing reporting no step; and the four
`mercator_pacific.grib2` / `polar_stereographic_alaska.grib2` ones sit 0.05 off the value they
copied, because the re-packed differenced result message uses a coarser lattice (step 0.25 and
0.125 against the inputs' 0.1) than the input, 0.2 and 0.4 of the result's own step, inside the
0.125 and 0.0625 half step its packing allows
(`_assert_values_from_input`). There is no tolerance for "near
enough": the distance between two input values is not a packing effect, and half of it would
accept values MIR never copied. 62 of the 115 outcome-recording triples have more than one field
(56 of the 109 `Expected` ones, and all 6 `KnownBug` ones), and one npz holds them all: `arr_0` is
the field values stacked into an array of shape `(n_fields, *field shape)`.

**Nothing of the input is read before the regridding.** The matrix test builds the fieldlist,
calls `regrid` and reads values only once that is done, because unpacking a file's own values
can fail with an error of its own while the table records the regridding's. That is how the two
files with self-contradicting grid metadata behave: reading `regular_gg_wrong_increment.grib`
raises a `ValueError` about a reshape, while regridding it raises the `gribapi.errors.WrongGridError`
its `INVALID` entry records (`regular_ll_wrong_increment.grib` reads fine, only MIR refuses it).
They are the negative set below now, and its test keeps the same order for the same reason.

**The table records what MIR did, odd or not — except where MIR is provably wrong.** A result
onto a target that is the input's own frame came back as a no-op copy that kept the rotation
while a real interpolation drops it (§C): that is a **MIR bug**, not MIR semantics — the target
is unrotated, so the correct output is unrotated — and the six affected triples are `KnownBug`
entries asserting the correct result, not `Expected` ones (see below). The remaining oddities
are recorded as observed: a
sub-area source clips the target to its area — 2 x 3 of the 684 nodes of `10x10` for
`alternate-scanning.grib`, 6 x 19 with an all-NaN pole row for the rotated sub-area Gaussians,
1 x 2 spelled `{"grid": [10, 0]}` for `rotated_ll_subarea.grib2` and for `icon_ch1.grib2` (§D), the
ICON one also coming back from `O32` as a 13-node clip whose spec reads `grid: "N32"` together
with that area (§G); 18 triples come back with
NaN, and each of them is the input's own missing points — the `9999` fill of `mercator.grib`, 3.43 M
of its 4.51 M points, the 40% of `polar_stereographic_north.grib1` its bitmap marks missing, or the
140 and 1596 nodes of `polar_stereographic_alaska.grib2`'s nearest-neighbour results whose
spherical-nearest input point is one of the 73 its bitmap masks (legitimate NaN, §F) —, the target outside
the source area, or MIR finding no element among the scattered points of a projected input (682 of
684 for `linear` onto `lambert_conformal.grib`, 491 for `lambert_conformal_na.grib2`, 341 for
`mercator_pacific.grib2` and 522 for `polar_stereographic_alaska.grib2`, onto `10x10`);
nearest-neighbour onto those projected inputs does not leave the far nodes NaN but fills the whole
globe with regional values — 0 NaN for `lambert_conformal_na.grib2` and `mercator_pacific.grib2`,
§F — and `eORCA025_T.grib`, whose bitmap masks 44,442 of its points, comes back with only 36 of 684
and 25 of 5248 NaN because MIR silently resamples its unstructured points onto O400 first and each
target node snaps to a finite O400 neighbour unless only NaN-ring ones are near (§H); and the six projection files
(`lambert_conformal.grib`, `lambert_conformal_na.grib2`, `mercator.grib`, `mercator_pacific.grib2`,
`polar_stereographic_alaska.grib2`, `polar_stereographic_north.grib1`) answer onto `O32` with the
target's own gridspec, `{"grid": "O32"}`: the array path re-encodes a message ecCodes can read,
while a projected input
stays on the array path because MIR's GRIB path cannot regrid it (§E). The `_field_points`
fallback that says a result's grid cannot be read back is a guard: no case of the table needs it.
Each of them is an `Expected` entry, commented with the section of the regrid-grib findings (§A–§H)
that explains it.

**5 BROKEN triples**, of 3 files, in two groups:

- the two **one-point** files (4): no stencil can be built from a single node and MIR asserts
  instead of reporting a degenerate grid, §B;
- `linear-mercator.grib-10x10` (1): not a hang, just 233 s the first time, while MIR builds the
  weights of 4.5 M scattered points (2.4 s with its weight cache) — too slow to be a test, §E.

**The negative set** is what became of the first group the matrix used to have — the two
**wrong-increment** files, all three combinations each (6): each contradicts its own grid
metadata, so no regridding of it can be right and a reference would be meaningless.
`regular_gg_wrong_increment.grib` says `Ni*Nj` = 192\*64 = 12288 for 18432 data points and ecCodes
itself refuses the grid — `gribapi.errors.WrongGridError: Grid description is wrong or inconsistent`; `regular_ll_wrong_increment.grib` encodes 2° increments on a 5° grid and MIR's own
GRIB reader asserts on it — `RuntimeError: Assertion failed: ni == ni_`. The files are like that
on purpose, §A; responsible: the data file, in both cases. They live in
`regrid_grib_cases.INVALID_GRIB_FILES` instead of `GRIB_FILES`, with no reference and never
probed by the generator; `INVALID` records the exception, message and reason per combination, and
`test_regrid_fieldlist_grib_invalid_file` runs the same three combinations on each file and
asserts the regridding fails with `pytest.raises` — the refusal comes from the grid, not from the
interpolation or the target.

There used to be a fourth group, the two **spectral** files, all three combinations each (6):
`AttributeError: 'NoneType' object has no attribute 'shape'`, because earthkit-geo asked such a
field for its latitudes — ecCodes gives no gridspec for a spectral truncation, and earthkit-data
raised on the `None` its own `latlons()` had no answer for — before MIR ever
saw the message, which regrids spectral GRIB natively. **Fixed in earthkit-geo on 2026-11-20**:
`regrid_grib_field` sends a field with no grid points straight to the backend's GRIB path;
the six triples are `Expected` entries with references now, plus two fast regression tests.
(`latlons()` answers `(None, None)` for a grid with no coordinates now.)

Four of the five BROKEN triples are marked `xfail(raises=<exception>, strict=True)` and **really
run**, so a fix turns into an `XPASS` and a failure that says so — unlike the numpy matrix, where
a BROKEN pair is `xfail(run=False)` and stays invisible
([Combinations MIR cannot do (BROKEN)](#combinations-mir-cannot-do-broken)). Only the slow
mercator one is `run=False`; and unlike those 4 it does produce a result (19 x 36,
`regular_ll`, spec `{"grid": [10, 10]}`, 595 NaN — in the comment above its entry),
but no reference is kept for a triple that never runs: the files check counts exactly one npz per
triple that runs — 115 files: 109 for the `Expected` entries plus the 6 stale ones the `KnownBug`
triples keep (see next paragraph, and [The files](#the-grib-matrix); 847 KB). Every long case
carries `pytest.mark.timeout(300, method="thread")`, so a MIR hang ends the run with a message
instead of stalling it.

**6 KNOWN MIR BUG triples**, of 4 files, all one bug (§C): regridding a **rotated** input onto
the **bare target that spells its own geometry** — the rotated 10x10 `rotated_ll.grib1` and
`rotated_ll_subarea.grib1` onto `{"grid": [10, 10]}`, the rotated
`reduced_rotated_gg_O32.grib1` and `reduced_rotated_gg_subarea_O32.grib1` onto `{"grid": "O32"}`
— comes back **byte-for-byte the input, rotation kept**, reproduced with plain mir-python without
earthkit-geo. Its six are the two `10x10` files with nn and linear plus the two `O32` files with
nn, because `GRIB_COMBINATIONS` runs only nearest-neighbour onto `O32` (see above); the bug is not
method-specific — MIR copies the rotated input back with linear too. The target is unrotated, so
the correct output is unrotated: this is a MIR bug, not MIR semantics. (A target of *other*
geometry interpolates and drops the rotation correctly, and the no-op copies of the *unrotated*
`regular_ll_subarea.grib1` / `regular_ll_single_point.grib1` onto `10x10` are correct identities —
their nodes already lie on the target grid — so those stay `Expected`.) Each of the six gets a
`KnownBug` entry recording the CORRECT result — full field lists for the two whole-grid files,
n_fields plus "unrotated" only for the two sub-area ones, whose post-fix clip the bug hides —
and `cases()` marks it `xfail(strict=True)`: no `raises=`, no `run=False`, the case really runs
and fails the correct-expectation checks while MIR is buggy. The day MIR is fixed it turns into
an `XPASS` and a failure that says so: the alert to regenerate the six references (their stored
npz is the stale no-op output, which the matrix test deliberately never compares against) and
paste the entries back as `Expected`. Two things to keep in mind: the strict mark alerts in the
fix direction only, so **any unrelated breakage of these six cases hides behind the expected
failure — read the xfail reasons in the report, not just the counts**; and since they fail before
the reference is read, they xfail identically with and without the reference tree.

**The files** are the `regrid/grib/` subtree of the same tree as the numpy data
([The data files](#the-data-files)): `<root>/regrid/grib/<interpolation>/<file>/<target>.npz`,
cached in `tests/data/regrid/grib/`, while the GRIB **inputs** keep coming from `test-data/grib/`
into `tests/data/grib/`. Each matrix checks its own half: `test_regrid_fieldlist_grib_files`
looks at `regrid/grib/` alone, and `test_regrid_numpy_grids_files` ignores it. Those 115 references
are **not on the download server yet**, so a plain long run currently fails on the 109 that are
actually compared: `FileNotFoundError`s naming the generator, next to 628 passes
([Running the tests](#running-the-tests)) — the six `KnownBug` references are not read and do not
add failures, and they are stale by design until the §C bug is fixed. Until they are uploaded, a
long GRIB run has to point
at the generated tree instead — `EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data`, which is green.

### Adding a GRIB file

1. Add its name to `GRIB_FILES` in `regrid_grib_cases.py`. It must exist under `URL_ROOT/grib/` on
   the download server: that is where `data_path("legacy", <name>, subfolder="grib")` fetches it
   from. The matrix regrids the whole fieldlist it holds, so a file with several fields or several
   grids needs no extra work.

1. Generate its references, from the repository root (every triple in its own subprocess, because
   MIR can abort or hang; `--timeout` is per triple, 120 s by default):

   ```
   python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file <new-name.grib>
   ```

1. Paste what it prints into `EXPECTED`: one `Expected(...)` line per triple that produced a
   result — n_fields, and per field the shape, gridType, gridspec and NaN count — or a
   `Broken(kind="raises", exception=..., reason="BROKEN: ...")` entry for one that did not. The
   summary's `raised` becomes a `Broken(kind="raises")` with the exception type it names; a `hang`
   or `abort` becomes `Broken(kind="hang")`, whose case is then never run. The same lines are in
   the JSON report, `<root>/regrid-grib-expected.json` (`--report` chooses where, `-` prints it):
   the table is the record, the generator only observes it. One exception: the generator observes
   the MIR *result*, so for the six `KnownBug` triples it prints `Expected` lines for the buggy
   no-op copy — never paste those over a `KnownBug` entry; a `KnownBug` entry is only replaced,
   and its stale reference regenerated, once MIR is fixed and the case XPASSes.

1. Run it against the local tree — until the upload the server has no reference for it, and a case
   that cannot read one fails:

   ```
   EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
       tests/regrid/regrid_mir/test_regrid_fieldlist_grib.py -E long --no-cov -q -k <new-name>
   ```

1. Upload the new `regrid/grib/` files to `URL_ROOT/regrid/grib/` — the maintainer's step, and the
   only way a plain long run goes green again. Keep the report JSON out of the upload: the tree is
   uploaded as it is, no sidecar files.

A file whose **own grid metadata contradicts itself** — no regridding of it can be right, §A, like
the two wrong-increment files — is not matrix material. Add its name to `INVALID_GRIB_FILES`
instead of `GRIB_FILES` (it must still exist under `URL_ROOT/grib/`), with one
`Invalid(exception=..., match=..., reason="INVALID INPUT: ...")` entry in `INVALID` per
`GRIB_COMBINATIONS` pair: the exception module-qualified unless a builtin, `match` a regex for the
message, `reason` naming the finding and the responsible side. Observe the errors by regridding the
file by hand before writing the entries — the generator does not probe negative files, they get no
references, and `test_regrid_fieldlist_grib_invalid_file` asserts each combination's failure with
`pytest.raises`. The fast table check keeps a file out of both lists at once and demands one
`INVALID` entry per combination of every negative file.

A new **target** is `GribTarget(id, spec)` in `GRIB_TARGETS`, a new **combination** an
`(interpolation, target id)` pair in `GRIB_COMBINATIONS` (the interpolation must be registered in
`regrid_common.INTERPOLATIONS`); `cases()` and the generator both read those lists, new pairs have
no references and no `EXPECTED` entries until the generator has run for them, and each negative
file needs one `INVALID` entry for every new pair.

### Regenerating the GRIB references

```
python tests/regrid/regrid_mir/generate_regrid_grib_data.py --help
python tests/regrid/regrid_mir/generate_regrid_grib_data.py                 # all 120 triples
python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file regular_ll.grib1 --interpolation nn
python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file regular_ll.grib1 --output-dir /tmp/gribdata
```

The same output-root logic as the numpy generator, which it imports: `--output-dir PATH`, else
`EARTHKIT_GEO_REGRID_GENERATE_DIR`, else `_dev/test-data`; the interpolation filter takes aliases
(`nn`); `--force` overwrites references that exist, and a triple that is skipped for want of it has
no new observation in the report. `--target` narrows a run to one target. The script never uploads
and never downloads a reference — it sets `EARTHKIT_GEO_REGRID_TEST_DATA` to its own output root,
so its `grib_reference` reads are local — while the GRIB inputs keep coming from the server into
`tests/data/grib/`. Its exit code is 0 even when triples raise: a failure is a finding for the
table, and the run prints them.

Regenerate when a new MIR version changes results, exactly as for the numpy matrix
([Regenerating](#regenerating)): `--force` over what it touched, re-run the tests through the
override, and each failure is either a change to accept or a MIR bug to report — then paste the new
`Expected` lines over the old ones and upload.

## Known behaviour of MIR in the matrix

From running every pair of the table once (`.ecmwf-multi-agents/notes/regrid-matrix-probe.md`,
2026-11-19). These are recorded by the tests, not fixed here.

- **A regional source clips the output.** Onto a regular lat-lon or a reduced Gaussian target,
  MIR returns only the part of the target inside the source area: `swisslv95` onto `O32` gives
  1 point, `LAEA-EFAS-5km` onto `ll-10x10` gives 5 x 7. All 34 such pairs are listed in
  `regrid_cases.CLIPPED` with the gridspec they return, and the test asserts that gridspec
  instead of the target. Only the 12 covered regional sources ever clip.
- **A HEALPix target is not clipped: it comes back whole, NaN-filled outside the source
  area** (12 pairs, 698 to 767 of the 768 values NaN). The nearest comparison therefore first
  requires result and reference to have NaN at exactly the same positions, and only then
  compares the finite values.
- **Nearest-neighbour ties are machine dependent.** A target point that is exactly between
  several source points has more than one nearest neighbour, and MIR picks one: `ll-reference`
  onto `ll-1x1-area` has all 3,408 target points in that situation (a 1° lattice of integer
  coordinates over a source of cell centres), so every one of its values may legitimately come
  out differently elsewhere. The comparison is strict and the difference is then checked, not
  tolerated, see [How a result is compared](#how-a-result-is-compared).
- **A clip that is too small for the target resolution gets a degenerate gridspec**:
  `ICON-CH2` onto `ll-10x10` returns 1 x 2 as `{"grid": [10, 0]}` and
  `xy-laea-grid-bbox-reverse` returns 1 x 1 as `{"grid": [0, 0]}`. A spacing of 0 is not a
  legal regular lat-lon resolution; the point counts are right, the spelling is not. Such a
  result is compared **whole**: the pole-row exemption of `linear` and `grid-box-average` is not
  derived from the shape, and a clipped pair never gets it — one row of values is exactly the
  result a regional source produced, and it is the values that matter.
- **`nearest-neighbour` onto a coarse global target fails when no target point lies in the
  source area**: the three Swiss sources onto `ll-10x10` raise
  `Serious bug: KNearestNeighbours: failed to interpolate` (3 BROKEN pairs). It is the target
  lattice missing the source entirely, not "regional onto global" and not the projection: the
  same sources succeed onto `O32`, `H8-ring` and `ll-1x1-area`.
- **`grid-box-average` needs grid boxes, which the unstructured meshes do not have**: `ORCA2_T`
  and `pi_C` raise `Representation::gridBoxes() not implemented` (4 BROKEN pairs). Expect the
  same of any ORCA, FESOM or ICON source if that selection grows; nearest-neighbour and linear
  are fine on those grids.
- **The returned gridspec can keep the rotation of the source**: the probe found three pairs
  doing it; one is left — `O32-rotated` onto `O32` returns the target shape described as a
  rotated grid, which is a different grid than the one asked for. Its values pass every check,
  so it lives in `CLIPPED` with the rotated spec and stays visible as a regression test rather
  than as an xfail. The other two, `ll-rotation` and `ll-rotated_ll-alias` onto `ll-1x1-area`,
  stopped leaking the rotation when the backend began passing the target's area inside the grid
  spec (2026-11-26): they return the target spec itself, are ordinary matrix pairs with no
  `CLIPPED` entry, and their references were regenerated for the new result. The GRIB matrix's
  six same-geometry rotated cases are the opposite situation: MIR does not interpolate at all but
  copies the rotated input back
  unchanged, so the output grid *and* values are wrong — that is a verified MIR bug, recorded as
  strict-xfail `KnownBug` entries ([The GRIB matrix](#the-grib-matrix)).
- **Onto a reduced Gaussian target, a clip is not exactly the target points inside the source
  box**: `O32-area` onto `O32` returns 274 points where 250 lie inside the source box, because
  the area cut keeps whole Gaussian rows. `CLIPPED` records what comes back.
- Nothing aborted and nothing hung in the 412 pairs; the slowest took 47.6 s (a 717k-point
  ICON grid onto `O32`).

## Other data files

The older tests of this folder read their reference npz files through
`regrid_common.data_path("legacy", ...)`, which delegates to
`earthkit.geo.utils.testing.get_test_data`: files under `tests/data/` in the subfolders
`global_0_360/`, `orca/`, `local/` and `xr/`, which that helper fetches into the subfolder when
a file is missing. (The GRIB inputs of the fieldlist tests are streamed from the same server
with `get_test_data_path`, and `tests/data/o32.grib2` is tracked in git.) Those files are
**shared with `tests/regrid/regrid_backends/` and `tests/xarray_grid/`**: do not rename, move,
regenerate or edit them, and do not change that helper for existing callers. The new matrix
data is kept in its own tree (`regrid/` on the server, `regrid/` in the cache
`tests/data/regrid/`) precisely to stay clear of them.
