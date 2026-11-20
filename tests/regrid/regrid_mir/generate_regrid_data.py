# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Generate the input grids and the MIR reference results of the regrid matrix tests.

Not a test file; run by hand from the repository root (never during pytest):

    python tests/regrid/regrid_mir/generate_regrid_data.py [options]

It writes the new data hierarchy of ``regrid_common.data_path`` under
``<output root>/regrid/`` (output root: ``--output-dir``, default
``_dev/test-data``): an input per covered source to
``regrid/input/<source>.npz`` and a reference per (interpolation, source, target)
pair to ``regrid/<interpolation>/<source>/<target>.npz``, each a compressed npz
with the single array key ``arr_0``. The relative layout mirrors ``data_path``
(the :func:`_input_relpath` / :func:`_reference_relpath` below are the generator
side of that one-place path logic). The tree is then uploaded to the download
server (``https://sites.ecmwf.int/repository/earthkit-geo/test-data/regrid/``,
the maintainer's step); the tests fetch it from there into their cache
``tests/data/regrid/``. To run the tests on freshly generated files BEFORE the
upload, run pytest with ``EARTHKIT_GEO_REGRID_TEST_DATA=<output root>`` (in the
default setup that is ``_dev/test-data``) —
``data_path`` then reads this script's output folder instead of downloading
(this script sets that variable for its own process and its workers anyway, so
generating never goes to the network). No sidecar files go into the tree -- it
is uploaded as it is -- the per-pair metadata (returned gridspec, shapes,
seconds) is printed by the workers as JSON lines and saved in a run report at
``<output root>/regrid-run-report.json``, outside ``regrid/``.

Output root
-----------
The output root is the folder that contains the ``regrid/`` tree the script
writes. ``--output-dir PATH`` sets it; without the option the environment
variable ``EARTHKIT_GEO_REGRID_GENERATE_DIR`` is the default, and without that
either the built-in default ``_dev/test-data``. A relative path is taken
relative to the current working directory, and the folder is created when it
does not exist (a ``--dry-run`` writes nothing and creates nothing). The option
wins over the variable, the variable over the built-in default. The root is
resolved once in ``main()`` and passed down; every worker subprocess gets the
same root (``_run_pair`` passes ``--output-dir <root>`` on its command line,
and the same folder ends up in the ``EARTHKIT_GEO_REGRID_TEST_DATA`` the parent
sets), so a custom root never sends the generator's own reads to the download
server.

Usage
-----
Everything (94 inputs, 376 + 18 + 18 reference pairs; takes hours, run parts of
it with the filters below):

    python tests/regrid/regrid_mir/generate_regrid_data.py

One source: its input, and every reference of every interpolation that runs it:

    python tests/regrid/regrid_mir/generate_regrid_data.py --source O32

One pair:

    python tests/regrid/regrid_mir/generate_regrid_data.py --source pi_C --target O32 --interpolation linear

Inputs only, four parallel workers, overwriting existing files:

    python tests/regrid/regrid_mir/generate_regrid_data.py --inputs-only --jobs 4 --force

Print what a run would do (counts, the output root in use and one line per
file) without writing:

    python tests/regrid/regrid_mir/generate_regrid_data.py --dry-run

Write somewhere other than ``_dev/test-data`` (the folder is created, relative
paths count from the current working directory; the same default can come from
``EARTHKIT_GEO_REGRID_GENERATE_DIR``):

    python tests/regrid/regrid_mir/generate_regrid_data.py --output-dir /tmp/regrid-data --source O32

``--source``, ``--target`` and ``--interpolation`` are repeatable and filter the
run (``--interpolation`` takes canonical names and aliases alike; with no
filters it generates inputs for every covered source and references for every
pair of ``regrid_cases.matrix(interpolation)`` for every interpolation in
``SELECTION``). Existing files are never overwritten without ``--force``. The
exit code is 0 even when pairs fail -- a failure is a finding for the probe and
for ``BROKEN`` in ``regrid_cases.py`` -- and non-zero only for usage errors.

Options
-------
``--output-dir PATH``   the folder that will contain the ``regrid/`` tree
                        (default: ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, or
                        ``_dev/test-data``; relative to the current working
                        directory, created if missing; the option wins over
                        the variable).
``--source ID``         only this source (repeatable); must be covered.
``--target ID``         only this target (repeatable).
``--interpolation NAME`` only this interpolation (repeatable, alias ok).
``--inputs-only``       write the input files and nothing else.
``--force``             overwrite files that already exist.
``--timeout SECONDS``   per-pair subprocess timeout (default 120).
``--jobs N``            run N pair subprocesses at a time (default 1, serial).
``--include-broken``    also run pairs listed in ``BROKEN`` (for the probe).
``--dry-run``           print what would be done, write nothing.
``--worker``            (internal) regrid exactly one pair, see "Safety".

Adding data for a new grid, target or interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Add the source to ``regrid_cases.SOURCES`` (or the target to ``TARGETS``, or,
for a new interpolation, the registry entry in ``regrid_common.INTERPOLATIONS``
plus the ``SELECTION`` entry -- each module docstring explains how), then run
the generator for just the new data and check the report:

    python tests/regrid/regrid_mir/generate_regrid_data.py --source <new-id>
    python tests/regrid/regrid_mir/generate_regrid_data.py --interpolation <new-name>

Test the new files through the override BEFORE the upload (the tests would
otherwise download the old tree from the server), and once they pass, upload
the regenerated ``regrid/`` tree — from whichever output folder was used — to
the server; that step is the maintainer's, this script never uploads:

    EARTHKIT_GEO_REGRID_TEST_DATA=_dev/test-data python -m pytest \
        tests/regrid/regrid_mir/test_regrid_numpy_grids.py -E long --no-cov -k <new-id>

(With ``--output-dir`` or ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, point
``EARTHKIT_GEO_REGRID_TEST_DATA`` at that folder instead of ``_dev/test-data``.)

Safety
^^^^^^
MIR can abort the interpreter (an uncaught C++ ``eckit::NotImplemented`` =
SIGABRT) or hang on some grid combinations, so every pair is regridded in a
separate subprocess (this very script, restarted with ``--worker``) guarded by
``--timeout``; an abort, timeout or exception is recorded in the results and the
reference file is not written (the worker writes to a temporary file and renames
it only after a successful regridding). Inputs are computed in the parent from
``Grid(spec).to_latlons()`` for covered sources only, which by the table's own
rule means below ``MAX_POINTS`` points. ``Grid.order`` is never read (it aborts
for ``sh``/ORCA/ICON/FESOM/unstructured grids) and no coordinates are ever built
for a source that is not covered.

Regenerating an input with ``--force`` invalidates the references computed from
the old values: after ``--force`` on inputs, regenerate the references that
depend on them (the same command without ``--inputs-only`` skips nothing when
``--force`` is given, so use the filters to keep it cheap).
"""

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STATUSES = ("ok", "skipped", "failed", "aborted", "timeout")

# The default generation output root: the folder whose regrid/ tree this script writes
# (data_path's relative layout below it). --output-dir overrides it, then the environment
# variable GENERATE_DIR_ENV below, then this constant. main() resolves the root once, sets
# EARTHKIT_GEO_REGRID_TEST_DATA to it for the process and, through the inherited environment,
# for every worker subprocess, and _run_pair also passes it to each worker as --output-dir,
# so the generator's own reads never go to the download server, whatever the root is.
# The tree is uploaded to the server by the maintainer; see the module docstring.
GENERATED_DATA_ROOT = Path(__file__).resolve().parents[3] / "_dev" / "test-data"

# Environment variable holding the default output root for runs without --output-dir
# (the option wins over the variable, the variable over GENERATED_DATA_ROOT).
GENERATE_DIR_ENV = "EARTHKIT_GEO_REGRID_GENERATE_DIR"


def _resolve_output_root(option):
    """The output root as an absolute path: --output-dir, else the env var, else the default.

    ``option`` is the value of --output-dir or None. Relative paths are taken
    relative to the current working directory. Nothing is created here: the
    folder is made by :func:`_execute` unless the run is a --dry-run.
    """
    value = option or os.environ.get(GENERATE_DIR_ENV) or str(GENERATED_DATA_ROOT)
    return Path(value).expanduser().resolve()


# The file name suffixes and the tree layout, as regrid_common.data_path builds them
# for the new hierarchy. Keep the three _*relpath helpers in step with that one function.
REPORT_NAME = "regrid-run-report.json"


def _npz(name: str) -> str:
    """Return ``name`` with a guaranteed ``.npz`` suffix (as data_path does)."""
    return name if name.endswith(".npz") else f"{name}.npz"


def _input_relpath(source_id: str) -> Path:
    """The path of ``data_path("input", source_id)`` relative to the data root."""
    return Path("regrid", "input", _npz(source_id))


def _reference_relpath(interpolation: str, source_id: str, target_id: str) -> Path:
    """The path of ``data_path("reference", ...)`` relative to the data root.

    ``interpolation`` must already be the canonical name, which is what
    ``data_path`` would make of it through ``get_interpolation``.
    """
    return Path("regrid", interpolation, source_id, _npz(target_id))


def smooth_field(lat, lon):
    """The deterministic test field: a smooth function of (lat, lon) in degrees.

    With ``phi`` and ``lam`` the latitude and longitude converted to radians::

        value = 15 + 10*sin(2*phi) + 4*cos(3*lam)*cos(phi) + 1.5*sin(5*phi)*sin(2*lam)

    float64, finite and positive everywhere (about 1.2 to 28.8), no NaN. The
    value of a point depends on its coordinates only, never on the order the
    points come in, so two runs write the same file. A nearest-neighbour result
    is then a value of the input, and the field varies on all scales the targets
    have, without the discontinuities that would make linear and grid-box-average
    results ill-conditioned.
    """
    import numpy as np

    phi = np.radians(np.asarray(lat, dtype="float64"))
    lam = np.radians(np.asarray(lon, dtype="float64"))
    return 15.0 + 10.0 * np.sin(2 * phi) + 4.0 * np.cos(3 * lam) * np.cos(phi) + 1.5 * np.sin(5 * phi) * np.sin(2 * lam)


def _json_safe(obj):
    """Convert a gridspec (or any MIR output) to JSON-safe types.

    numpy scalars and arrays go through ``tolist()``, tuples and sets become
    lists, dicts keep their values and get string keys, anything else is
    stringified so the JSON line can always be printed.
    """
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if hasattr(obj, "tolist"):
        return _json_safe(obj.tolist())
    if isinstance(obj, dict):
        return {str(key): _json_safe(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple, set, frozenset)):
        return [_json_safe(item) for item in obj]
    return str(obj)


def _last_json_line(text: str):
    """Return the last parseable JSON object printed on ``text``, or None."""
    for line in reversed((text or "").splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                payload = json.loads(line)
            except ValueError:
                continue
            if isinstance(payload, dict):
                return payload
    return None


def _last_line(text: str) -> str:
    """The last non-empty line of worker stderr, the crash message of an abort."""
    for line in reversed((text or "").splitlines()):
        if line.strip():
            return line.strip()
    return ""


# --------------------------------------------------------------------------------------
# the parent: plan, inputs, workers, report
# --------------------------------------------------------------------------------------


def _canonical_interpolations(names, parser):
    """Map the requested interpolation names to canonical names, or use all of SELECTION."""
    from regrid_cases import SELECTION
    from regrid_common import get_interpolation

    if not names:
        return list(SELECTION)
    canonical = []
    for name in names:
        try:
            interpolation = get_interpolation(name)
        except ValueError as error:
            parser.error(str(error))
        if interpolation.name not in SELECTION:
            parser.error(
                f"interpolation {interpolation.name!r} is registered but has no SELECTION entry "
                "in regrid_cases.py, so the generator does not know which pairs to run for it"
            )
        if interpolation.name not in canonical:
            canonical.append(interpolation.name)
    return canonical


def _selected_ids(given, entries, what, parser):
    """Validate the --source/--target ids against the table and return the id list to use."""
    by_id = {entry.id: entry for entry in entries}
    ids = list(dict.fromkeys(given or []))
    for id in ids:
        entry = by_id.get(id)
        if entry is None:
            parser.error(f"unknown {what} {id!r} for --{what}: not in regrid_cases.py")
        covered = getattr(entry, "covered", True)
        if not covered:
            parser.error(f"{what} {id!r} is not a matrix {what}: {entry.reason}")
    if ids:
        return ids
    return [entry.id for entry in entries if getattr(entry, "covered", True)]


def _plan(interpolations, selected_sources, selected_targets, root, force, include_broken):
    """The pairs of the run as dicts, with status "run", "skipped" or reason broken/existing."""
    from regrid_cases import BROKEN, matrix

    source_ids = {source.id for source in selected_sources}
    target_ids = {target.id for target in selected_targets}
    pairs = []
    for interpolation in interpolations:
        for source, target in matrix(interpolation):
            if source.id not in source_ids or target.id not in target_ids:
                continue
            key = (interpolation, source.id, target.id)
            path = root / _reference_relpath(interpolation, source.id, target.id)
            pair = {
                "kind": "reference",
                "interpolation": interpolation,
                "source": source.id,
                "target": target.id,
                "path": str(path.relative_to(root)),
            }
            if key in BROKEN and not include_broken:
                pair.update(status="skipped", reason="broken", message=f"skipped: BROKEN ({BROKEN[key]})")
            elif path.exists() and not force:
                pair.update(
                    status="skipped", reason="existing", message="skipped, the file exists (--force to overwrite)"
                )
            else:
                pair["status"] = "run"
            pairs.append(pair)
    return pairs


def _input_ids(args, pairs, source_ids):
    """The sources whose input files this run writes.

    With --inputs-only the (filtered) selection itself; otherwise the sources of the
    planned pairs, skipped ones included -- without filters that is every covered
    source, and a pair that is skipped today still tells us which input belongs here
    once BROKEN changes.
    """
    if args.inputs_only:
        return list(source_ids)
    ids = []
    for pair in pairs:
        if pair["source"] not in ids:
            ids.append(pair["source"])
    return ids


def _make_input(source, root, force):
    """Build and save the smooth input field of one covered source; return the result dict."""
    from regrid_common import save_values

    result = {"kind": "input", "source": source.id, "path": str(_input_relpath(source.id))}
    path = root / _input_relpath(source.id)
    if path.exists() and not force:
        result.update(status="skipped", reason="existing", message="skipped, the file exists (--force to overwrite)")
        return result

    from earthkit.geo.grids import Grid  # the heavy import, lazy on purpose

    start = time.perf_counter()
    try:
        grid = Grid(source.spec)
        lat, lon = grid.to_latlons()
        values = smooth_field(lat, lon).reshape(grid.shape)
        save_values(path, values)
    except Exception as error:  # a broken source is a finding, not a crash of the run
        result.update(status="failed", error=f"{type(error).__name__}: {error}", seconds=time.perf_counter() - start)
        return result
    result.update(status="ok", shape=list(values.shape), seconds=time.perf_counter() - start)
    return result


def _run_pair(pair, root, timeout):
    """Regrid one pair in a worker subprocess and turn its JSON line into a result dict."""
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--interpolation",
        pair["interpolation"],
        "--source",
        pair["source"],
        "--target",
        pair["target"],
        # the worker writes its reference under the same root the parent planned with
        "--output-dir",
        str(root),
    ]
    result = dict(pair)
    start = time.perf_counter()
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        result.update(status="timeout", error=f"timeout: killed after {timeout} s")
        return result

    payload = _last_json_line(completed.stdout)
    seconds = time.perf_counter() - start
    if payload and payload.get("status") == "ok":
        result.update(status="ok")
        for key in ("shape", "gridspec", "gridspec_shape", "seconds"):
            if key in payload:
                result[key] = payload[key]
        result.setdefault("seconds", round(seconds, 3))
    elif payload and payload.get("status") == "error":
        result.update(status="failed", error=payload.get("error", "the worker reported an error"))
    elif completed.returncode < 0:
        result.update(
            status="aborted",
            error=f"aborted: signal {-completed.returncode}: {_last_line(completed.stderr)} (no reference written)",
        )
    else:
        result.update(status="failed", error=_last_line(completed.stderr) or f"worker exited {completed.returncode}")
    return result


def _sweep_temp_files(root):
    """Remove worker leftovers (``*.tmp.npz``) from the data tree; returns their names."""
    leftovers = []
    regrid_root = root / "regrid"
    if regrid_root.exists():
        for path in regrid_root.rglob("*.tmp.npz"):
            try:
                path.unlink()
                leftovers.append(str(path.relative_to(root)))
            except OSError:
                pass
    return leftovers


def _write_report(root, results, counts, timeout):
    """Save the run report as JSON outside the regrid/ tree; return its path."""
    report = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command": " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
        "python": sys.version,
        "timeout": timeout,
        "counts": counts,
        "results": results,
    }
    path = root / REPORT_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as stream:
        json.dump(report, stream, indent=1)
        stream.write("\n")
    return path


def _print_results(results):
    """One line per result, the way the run discovers its findings."""
    for result in results:
        head = (
            f"input {result['source']}"
            if result["kind"] == "input"
            else f"reference {result['interpolation']} {result['source']} -> {result['target']}"
        )
        status = result["status"]
        if status == "ok":
            seconds = result.get("seconds")
            detail = f"shape {result.get('shape')}"
            if seconds is not None:
                detail += f", {seconds:.1f} s"
            print(f"{head}: ok ({detail})")
        else:
            message = result.get("message") or result.get("error") or status
            print(f"{head}: {status}: {message}")


def _summarise(results):
    """Print the counts and the failing pairs; return the counts dict."""
    counts = {status: 0 for status in STATUSES}
    skipped_broken = 0
    problems = []
    for result in results:
        counts[result["status"]] += 1
        if result["status"] == "skipped" and result.get("reason") == "broken":
            skipped_broken += 1
        if result["status"] not in ("ok", "skipped"):
            head = (
                f"input {result['source']}"
                if result["kind"] == "input"
                else f"{result['interpolation']} {result['source']} -> {result['target']}"
            )
            problems.append(f"  {head}: {result['status']}: {result.get('error', '')}")

    print(
        "summary: "
        + ", ".join(f"{counts[status]} {status}" for status in STATUSES)
        + f"  (skipped: {counts['skipped'] - skipped_broken} existing, {skipped_broken} broken)"
    )
    if problems:
        print("failing cases (findings for the probe / regrid_cases.BROKEN):")
        print("\n".join(problems))
    return counts


def _execute(args, parser, root):
    """The parent process: validate, plan, write the inputs, run the pair workers, report.

    ``root`` is the resolved output root (:func:`_resolve_output_root`); the ``regrid/``
    tree and the run report are written below it.
    """
    from regrid_cases import SOURCES, TARGETS

    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")

    interpolations = _canonical_interpolations(args.interpolation, parser)
    source_ids = _selected_ids(args.source, SOURCES, "source", parser)
    target_ids = _selected_ids(args.target, TARGETS, "target", parser)
    sources = {source.id: source for source in SOURCES}
    targets = {target.id: target for target in TARGETS}

    pairs = _plan(
        interpolations,
        [sources[id] for id in source_ids],
        [targets[id] for id in target_ids],
        root,
        args.force,
        args.include_broken,
    )
    if args.inputs_only:
        pairs = []  # write the input files and nothing else
    input_ids = _input_ids(args, pairs, source_ids)
    to_run = sum(pair["status"] == "run" for pair in pairs)

    print(f"output root: {root}")
    print(
        f"generator: {len(input_ids)} inputs, {to_run} reference pairs to run, "
        f"{len(pairs) - to_run} skipped, for {', '.join(interpolations)}"
    )

    if args.dry_run:
        for id in input_ids:
            exists = (root / _input_relpath(id)).exists()
            would = "would skip, the file exists" if (exists and not args.force) else "would write"
            print(f"input {id}: {would}")
        for pair in pairs:
            state = pair.get("message") or "would run"
            print(f"reference {pair['interpolation']} {pair['source']} -> {pair['target']}: {state}")
        per_interpolation = {name: sum(pair["interpolation"] == name for pair in pairs) for name in interpolations}
        print(f"dry run, nothing written to {root}:")
        print(
            f"  {len(input_ids)} inputs, "
            + " + ".join(f"{count} {name}" for name, count in per_interpolation.items())
            + f" = {len(pairs)} reference pairs"
        )
        return 0

    root.mkdir(parents=True, exist_ok=True)
    results = [_make_input(sources[id], root, args.force) for id in input_ids]
    running = [pair for pair in pairs if pair["status"] == "run"]
    executed = {}
    if running:
        if args.jobs > 1:
            with ThreadPoolExecutor(max_workers=min(args.jobs, len(running))) as pool:
                executed = {
                    (r["interpolation"], r["source"], r["target"]): r
                    for r in pool.map(lambda pair: _run_pair(pair, root, args.timeout), running)
                }
        else:
            executed = {
                (r["interpolation"], r["source"], r["target"]): r
                for r in (_run_pair(pair, root, args.timeout) for pair in running)
            }
    # every planned pair appears in the results, the skipped ones straight from the plan
    results += [executed.get((p["interpolation"], p["source"], p["target"]), p) for p in pairs]

    for leftover in _sweep_temp_files(root):
        print(f"cleaned up worker temp file {leftover}")

    _print_results(results)
    counts = _summarise(results)
    report = _write_report(root, results, counts, args.timeout)
    print(f"run report: {report}")
    return 0


# --------------------------------------------------------------------------------------
# the worker: exactly one pair, in its own process
# --------------------------------------------------------------------------------------


def _worker(args, parser, root):
    """Regrid one pair with MIR and print the one JSON line the parent collects.

    Restarted by :func:`_run_pair` as
    ``python generate_regrid_data.py --worker --interpolation I --source S --target T
    --output-dir ROOT``: the root travels on the command line, so the worker writes
    its reference and reads its input (through the ``EARTHKIT_GEO_REGRID_TEST_DATA``
    ``main()`` sets from the same root) under exactly the folder the parent planned.
    Loads the input through ``data_path`` (which names the generator command when it
    is missing), calls ``earthkit.geo.grids.array.regrid`` and saves the result to the
    reference path -- to a temporary file first, renamed only on success, so a crash
    mid-write cannot leave a half file in the data tree. Prints
    ``{"status": "ok", "shape": ..., "gridspec": ..., "gridspec_shape": ...,
    "seconds": ...}`` or ``{"status": "error", "error": ...}`` as the last stdout
    line; exits 0 for ok and 1 otherwise. An abort or a hang never reaches this code:
    the parent watches the return code and the subprocess timeout instead.
    """
    if not (args.interpolation and args.source and args.target):
        parser.error("--worker needs --interpolation, --source and --target")

    from regrid_cases import SOURCES, TARGETS
    from regrid_common import data_path, load_values, save_values

    interpolation = args.interpolation[-1]
    source = next((entry for entry in SOURCES if entry.id == args.source[-1]), None)
    target = next((entry for entry in TARGETS if entry.id == args.target[-1]), None)
    payload = {"status": "error", "interpolation": interpolation, "source": args.source[-1], "target": args.target[-1]}
    if source is None or target is None:
        payload["error"] = f"unknown source {args.source[-1]!r} or target {args.target[-1]!r}"
        print(json.dumps(payload), flush=True)
        return 1

    path = root / _reference_relpath(interpolation, source.id, target.id)
    start = time.perf_counter()
    try:
        import numpy as np

        from earthkit.geo.grids import Grid
        from earthkit.geo.grids.array import regrid as regrid_array

        values = load_values(data_path("input", source.id))
        result, returned = regrid_array(values, source.spec, target.spec, interpolation=interpolation)
        result = np.asarray(result)
        temp = path.with_name(f"{path.stem}.tmp.npz")
        save_values(temp, result)
        os.replace(temp, path)

        payload["shape"] = list(result.shape)
        payload["gridspec"] = _json_safe(returned)
        try:
            payload["gridspec_shape"] = list(Grid(returned).shape)
        except Exception as error:  # a gridspec that will not build is reported, not fatal
            payload["gridspec_shape"] = None
            payload["gridspec_error"] = f"{type(error).__name__}: {error}"
        payload["status"] = "ok"
    except Exception as error:
        payload["error"] = f"{type(error).__name__}: {error}"
    payload["seconds"] = round(time.perf_counter() - start, 3)
    print(json.dumps(payload, default=str), flush=True)
    return 0 if payload["status"] == "ok" else 1


# --------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------


def _build_parser():
    """The command line of the module docstring's Usage section."""
    parser = argparse.ArgumentParser(
        prog="python tests/regrid/regrid_mir/generate_regrid_data.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        metavar="PATH",
        help=(
            "root folder for the regrid/ tree and the run report (relative to the current working "
            f"directory, created if missing); default: ${GENERATE_DIR_ENV} or {GENERATED_DATA_ROOT}"
        ),
    )
    parser.add_argument("--source", action="append", metavar="ID", help="only this source (repeatable)")
    parser.add_argument("--target", action="append", metavar="ID", help="only this target (repeatable)")
    parser.add_argument(
        "--interpolation",
        action="append",
        metavar="NAME",
        help="only this interpolation, canonical name or alias (repeatable)",
    )
    parser.add_argument("--inputs-only", action="store_true", help="write the input files and no references")
    parser.add_argument("--force", action="store_true", help="overwrite files that already exist (default: skip them)")
    parser.add_argument(
        "--timeout", type=float, default=120.0, metavar="SECONDS", help="per-pair timeout (default 120)"
    )
    parser.add_argument(
        "--jobs", type=int, default=1, metavar="N", help="pair subprocesses to run at a time (default 1)"
    )
    parser.add_argument("--include-broken", action="store_true", help="also run the pairs listed in BROKEN")
    parser.add_argument("--dry-run", action="store_true", help="print what would be done, write nothing")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    return parser


def main(argv=None):
    """Parse the command line and run the parent or the worker.

    Resolves the output root (:func:`_resolve_output_root`: ``--output-dir``,
    else ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, else :data:`GENERATED_DATA_ROOT`)
    once, then sets ``EARTHKIT_GEO_REGRID_TEST_DATA`` to it before anything
    reads a data file: this process and (through the inherited environment)
    every worker subprocess then reads inputs and checks files in the
    generation output folder, never through the download server. The root is
    also passed to each worker as ``--output-dir`` (see :func:`_run_pair`).
    """
    from regrid_common import REGRID_TEST_DATA_ENV

    parser = _build_parser()
    args = parser.parse_args(argv)
    root = _resolve_output_root(args.output_dir)
    os.environ[REGRID_TEST_DATA_ENV] = str(root)
    if args.worker:
        return _worker(args, parser, root)
    return _execute(args, parser, root)


if __name__ == "__main__":
    sys.exit(main())
