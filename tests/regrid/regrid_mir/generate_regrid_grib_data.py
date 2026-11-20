# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Generate the MIR reference results of the GRIB regrid matrix tests.

Not a test file; run by hand from the repository root (never during pytest):

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py [options]

For every (interpolation, file, target) triple of ``regrid_grib_cases`` it reads the GRIB
file with ``earthkit.data.from_source("file", ...).to_fieldlist()`` -- the file comes from
the download server through ``regrid_common.data_path("legacy", filename, subfolder="grib")``,
which caches it in ``tests/data/grib/`` -- regrids it with ``earthkit.geo.regrid`` onto the
target gridspec and writes the result to

    <output root>/regrid/grib/<interpolation>/<filename>/<target id>.npz

the layout ``regrid_common.data_path("grib_reference", interpolation, filename, target_id)``
reads. ``arr_0`` holds the values of all result fields stacked into one array of shape
``(n_fields, *field shape)``; per field the run also observes ``shape``, ``gridType``,
``geography.grid_spec()`` (None where MIR returns no spec) and the NaN count, and those
observations are what ``regrid_grib_cases.EXPECTED`` records: the script prints them as
pasteable ``Expected(...)`` / ``Broken(...)`` lines and writes them as JSON, by default to
``<output root>/regrid-grib-expected.json`` (``--report`` to choose, ``-`` for stdout).
Keep the JSON out of the ``regrid/`` tree -- that tree is uploaded to the server as it is,
no sidecar files. Uploading the regenerated ``regrid/grib/`` subtree is the maintainer's
step; this script never uploads. Only the matrix is here: the two wrong-increment files of
``regrid_grib_cases.INVALID_GRIB_FILES`` are the negative set, tested for the error their
regridding raises, so they have no reference and ``--file`` does not accept them.

The output root is the folder that will contain the ``regrid/`` tree: ``--output-dir PATH``,
else ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, else ``_dev/test-data`` -- the same resolution
``generate_regrid_data.py`` uses (the script imports it), relative paths counting from the
current working directory, the folder created when missing. The variable
``EARTHKIT_GEO_REGRID_TEST_DATA`` is set to that root for this process and every worker, so
the generator's own ``grib_reference`` reads stay local whatever the root is; the legacy
GRIB inputs are unaffected by it and keep coming from the server.

Usage
-----
Everything (40 files x the three combinations; parts of it with the filters below):

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py

One file, or one interpolation of it (``--interpolation`` takes aliases; all three options
are repeatable, unknown names are refused):

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file regular_ll.grib2
    python tests/regrid/regrid_mir/generate_regrid_grib_data.py --file regular_ll.grib2 --interpolation nn

Write somewhere other than ``_dev/test-data``, overwriting existing references, with a
300 s per-triple timeout and the report to a fixed path:

    python tests/regrid/regrid_mir/generate_regrid_grib_data.py --output-dir /tmp/gribdata \
        --force --timeout 300 --report /tmp/gribdata/expected.json

Existing reference files are never overwritten without ``--force``; a skipped triple has no
new observation, so it is absent from the report. The exit code is 0 even when triples fail
-- a failure is a finding for ``regrid_grib_cases.EXPECTED`` -- and non-zero only for usage
errors.

Options
-------
``--output-dir PATH``   root folder for the regrid/ tree (default:
                        ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, or ``_dev/test-data``).
``--file NAME``         only this GRIB file (repeatable; must be in ``GRIB_FILES``).
``--interpolation NAME`` only this interpolation (repeatable, alias ok; must be part of
                        ``GRIB_COMBINATIONS``).
``--target ID``         only this target (repeatable; must be in ``GRIB_TARGETS``).
``--timeout SECONDS``   per-triple subprocess timeout (default 120).
``--force``             overwrite references that already exist.
``--report PATH``       where to write the observed expectations JSON (default
                        ``<output root>/regrid-grib-expected.json``; ``-`` prints it to
                        stdout instead of a file).
``--worker``            (internal) regrid exactly one triple, see "Safety".

Safety
^^^^^^
MIR can abort the interpreter (an uncaught C++ exception is a SIGABRT) or hang on some
grids, so every triple is regridded in a separate subprocess (this very script, restarted
with ``--worker``) guarded by ``--timeout``. A timeout is reported as a hang, a signal kill
or a worker that died without reporting as an abort, a caught exception as a raise with its
module-qualified type and message; in all three cases no reference is written (the worker
saves to a temporary file and renames it only after a successful regridding). The stacked
values require all result fields to have the same shape; a triple whose result fields differ
is reported as a raise of that ``ValueError`` -- a finding for the table, not a crash.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from generate_regrid_data import (
    GENERATE_DIR_ENV,
    GENERATED_DATA_ROOT,
    _json_safe,
    _last_json_line,
    _last_line,
    _npz,
    _resolve_output_root,
)

# The four observation statuses of a run, plus "skipped" for an existing reference. STATUSES
# is the order of the summary; the report maps them onto regrid_grib_cases entries: ok ->
# Expected, raised -> Broken(kind="raises"), hang/abort -> Broken(kind="hang").
STATUSES = ("ok", "raised", "hang", "abort", "skipped")

# Default name of the observed-expectations report, written outside the regrid/ tree.
REPORT_NAME = "regrid-grib-expected.json"


def _grib_reference_relpath(interpolation: str, filename: str, target_id: str) -> Path:
    """The path of ``data_path("grib_reference", ...)`` relative to the data root.

    ``interpolation`` must already be the canonical name, which is what ``data_path``
    would make of it through ``get_interpolation``. Keep in step with that one function.
    """
    return Path("regrid", "grib", interpolation, filename, _npz(target_id))


def _exception_name(error_type: type) -> str:
    """The exception type name as a table entry records it: module-qualified unless builtin."""
    if error_type.__module__ in ("builtins", "__main__"):
        return error_type.__qualname__
    return f"{error_type.__module__}.{error_type.__qualname__}"


# --------------------------------------------------------------------------------------
# the parent: plan, workers, report
# --------------------------------------------------------------------------------------


def _canonical_interpolations(names, parser):
    """The requested interpolation names, canonicalised, or all of GRIB_COMBINATIONS."""
    from regrid_common import get_interpolation
    from regrid_grib_cases import GRIB_COMBINATIONS

    wanted = list(dict.fromkeys(interpolation for interpolation, _ in GRIB_COMBINATIONS))
    if not names:
        return wanted
    canonical = []
    for name in names:
        try:
            interpolation = get_interpolation(name)
        except ValueError as error:
            parser.error(str(error))
        if interpolation.name not in wanted:
            parser.error(
                f"interpolation {interpolation.name!r} is not part of the GRIB matrix: "
                f"GRIB_COMBINATIONS runs {', '.join(wanted)}"
            )
        if interpolation.name not in canonical:
            canonical.append(interpolation.name)
    return canonical


def _selected_targets(names, parser):
    """The requested target ids, validated against GRIB_TARGETS, or all of them."""
    from regrid_grib_cases import GRIB_TARGETS

    known = [target.id for target in GRIB_TARGETS]
    for name in dict.fromkeys(names or []):
        if name not in known:
            parser.error(f"unknown target {name!r} for --target: GRIB_TARGETS has {', '.join(known)}")
    return [id for id in known if not names or id in set(names)]


def _selected_files(names, parser):
    """The requested file names, validated against GRIB_FILES, or all of them, in table order."""
    from regrid_grib_cases import GRIB_FILES

    for name in dict.fromkeys(names or []):
        if name not in GRIB_FILES:
            parser.error(f"unknown GRIB file {name!r} for --file: not in regrid_grib_cases.GRIB_FILES")
    return [name for name in GRIB_FILES if not names or name in set(names)]


def _plan(interpolations, target_ids, files, root, force):
    """The triples of the run as dicts, with status "run" or "skipped" (existing reference)."""
    from regrid_grib_cases import GRIB_COMBINATIONS

    pairs = []
    for interpolation, target_id in GRIB_COMBINATIONS:
        if interpolation not in interpolations or target_id not in target_ids:
            continue
        for filename in files:
            path = root / _grib_reference_relpath(interpolation, filename, target_id)
            pair = {
                "kind": "grib_reference",
                "interpolation": interpolation,
                "file": filename,
                "target": target_id,
                "path": str(path.relative_to(root)),
            }
            if path.exists() and not force:
                pair.update(
                    status="skipped", reason="existing", message="skipped, the file exists (--force to overwrite)"
                )
            else:
                pair["status"] = "run"
            pairs.append(pair)
    return pairs


def _run_triple(pair, root, timeout):
    """Regrid one triple in a worker subprocess and turn its JSON line into a result dict."""
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--interpolation",
        pair["interpolation"],
        "--file",
        pair["file"],
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
        result.update(status="hang", error=f"timeout: killed after {timeout} s, no reference written")
        return result

    payload = _last_json_line(completed.stdout)
    seconds = round(time.perf_counter() - start, 3)
    if payload and payload.get("status") == "ok":
        result.update(status="ok")
        for key in ("n_fields", "shapes", "grid_types", "grid_specs", "nan_counts"):
            result[key] = payload[key]
        result["seconds"] = payload.get("seconds", seconds)
    elif payload and payload.get("status") == "error":
        result.update(
            status="raised",
            error_type=payload.get("error_type", "Error"),
            error=payload.get("error", "the worker reported an error"),
            seconds=seconds,
        )
    elif completed.returncode < 0:
        result.update(
            status="abort",
            error=f"aborted: signal {-completed.returncode}: {_last_line(completed.stderr)} (no reference written)",
        )
    else:
        # a worker that died before reporting anything (a failed import, a bad argument):
        # nothing observed, so it hangs/aborts like a crash -- xfail(run=False), never raises=
        result.update(
            status="abort",
            error=f"worker exited {completed.returncode} without reporting: {_last_line(completed.stderr)}",
        )
    return result


def _expected_entries(results):
    """The ok/raised/hang/abort observations as EXPECTED-pasteable entries, nested by key."""
    entries = {}
    for result in results:
        if result["status"] == "ok":
            value = {
                "n_fields": result["n_fields"],
                "shapes": result["shapes"],
                "grid_types": result["grid_types"],
                "grid_specs": result["grid_specs"],
                "nan_counts": result["nan_counts"],
            }
        elif result["status"] == "raised":
            value = {
                "kind": "raises",
                "exception": result["error_type"],
                "reason": f"BROKEN: {result['error_type']}: {result['error']}",
            }
        else:  # hang / abort: never executed, no exception to name
            value = {"kind": "hang", "reason": f"BROKEN: {result['error']}"}
        entries.setdefault(result["interpolation"], {}).setdefault(result["file"], {})[result["target"]] = value
    return entries


def _python_entry(result):
    """One result as the pasteable ``regrid_grib_cases.EXPECTED`` line (a 4-space item)."""
    key = f'    ("{result["interpolation"]}", "{result["file"]}", "{result["target"]}"): '
    if result["status"] == "ok":
        return (
            key
            + "Expected("
            + f"n_fields={result['n_fields']}, "
            + f"shapes={tuple(tuple(shape) for shape in result['shapes'])}, "
            + f"grid_types={tuple(result['grid_types'])}, "
            + f"grid_specs={tuple(result['grid_specs'])}, "
            + f"nan_counts={tuple(result['nan_counts'])}),"
        )
    entry = _expected_entries([result])[result["interpolation"]][result["file"]][result["target"]]
    arguments = ", ".join(f"{name}={value!r}" for name, value in entry.items())
    return f"{key}Broken({arguments}),"


def _write_report(root, args, results, counts):
    """Save the observed expectations as JSON outside the regrid/ tree; return its path."""
    report = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "command": " ".join([Path(sys.argv[0]).name, *sys.argv[1:]]),
        "python": sys.version,
        "timeout": args.timeout,
        "counts": counts,
        "expected": _expected_entries([r for r in results if r["status"] != "skipped"]),
    }
    if args.report == "-":
        print(json.dumps(report, indent=1))
        return None
    path = Path(args.report).expanduser().resolve() if args.report else root / REPORT_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as stream:
        json.dump(report, stream, indent=1)
        stream.write("\n")
    return path


def _print_results(results):
    """One line per triple, the way the run discovers its findings."""
    for result in results:
        head = f"{result['interpolation']} {result['file']} -> {result['target']}"
        status = result["status"]
        if status == "ok":
            print(f"{head}: ok ({result['n_fields']} fields, shape {result['shapes'][0]}, {result['seconds']:.1f} s)")
        elif status == "skipped":
            print(f"{head}: {status}: {result['message']}")
        else:
            detail = result.get("error") or status
            if result.get("error_type"):
                detail = f"{result['error_type']}: {detail}"
            print(f"{head}: {status}: {detail}")


def _summarise(results):
    """Print the counts and the failing triples; return the counts dict."""
    counts = {status: 0 for status in STATUSES}
    for result in results:
        counts[result["status"]] += 1
    print("summary: " + ", ".join(f"{counts[status]} {status}" for status in STATUSES))
    problems = [
        f"  {result['interpolation']} {result['file']} -> {result['target']}: "
        f"{result['status']}: {result.get('error', '')}"
        for result in results
        if result["status"] in ("raised", "hang", "abort")
    ]
    if problems:
        print("failing triples (findings for the probe / regrid_grib_cases.EXPECTED):")
        print("\n".join(problems))
    return counts


def _execute(args, parser, root):
    """The parent process: validate, plan, run the triple workers, report, print the entries."""
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    interpolations = _canonical_interpolations(args.interpolation, parser)
    target_ids = _selected_targets(args.target, parser)
    files = _selected_files(args.file, parser)
    pairs = _plan(interpolations, target_ids, files, root, args.force)
    to_run = sum(pair["status"] == "run" for pair in pairs)

    print(f"output root: {root}")
    print(
        f"generator: {to_run} of {len(pairs)} (interpolation, file, target) triples to run, "
        f"{len(pairs) - to_run} skipped"
    )

    root.mkdir(parents=True, exist_ok=True)
    results = [pair if pair["status"] != "run" else _run_triple(pair, root, args.timeout) for pair in pairs]

    _print_results(results)
    counts = _summarise(results)
    report = _write_report(root, args, results, counts)
    print(f"expected report: {report if report else '(stdout)'}")

    observed = [result for result in results if result["status"] != "skipped"]
    if observed:
        print("paste into regrid_grib_cases.EXPECTED (keys are (interpolation, file, target) triples):")
        for result in observed:
            print(_python_entry(result))
    return 0


# --------------------------------------------------------------------------------------
# the worker: exactly one triple, in its own process
# --------------------------------------------------------------------------------------


def _worker(args, parser, root):
    """Regrid one triple with MIR and print the one JSON line the parent collects.

    Restarted by :func:`_run_triple` as
    ``python generate_regrid_grib_data.py --worker --interpolation I --file F --target T
    --output-dir ROOT``, the pattern of ``generate_regrid_data.py``. The input is fetched
    through ``data_path("legacy", filename, subfolder="grib")``, the fieldlist regridded
    with ``earthkit.geo.regrid``, and the stacked values of all result fields saved to the
    ``grib_reference`` path -- to a temporary file first, renamed only on success, so a
    crash mid-write cannot leave a half file in the data tree. Prints
    ``{"status": "ok", "n_fields": ..., "shapes": [...], "grid_types": [...],
    "grid_specs": [...], "nan_counts": [...], "seconds": ...}`` or
    ``{"status": "error", "error_type": ..., "error": ...}`` as the last stdout line; exits
    0 for ok and 1 otherwise. An abort or a hang never reaches this code: the parent
    watches the return code and the subprocess timeout instead.
    """
    if not (args.interpolation and args.file and args.target):
        parser.error("--worker needs --interpolation, --file and --target")

    interpolation = args.interpolation[-1]
    filename = args.file[-1]
    target_id = args.target[-1]
    payload = {"status": "error", "interpolation": interpolation, "file": filename, "target": target_id}

    start = time.perf_counter()
    try:
        import numpy as np
        from earthkit.data import from_source
        from regrid_common import data_path, get_interpolation, save_values
        from regrid_grib_cases import GRIB_TARGETS

        from earthkit.geo import regrid

        interpolation = get_interpolation(interpolation).name
        target = next((entry for entry in GRIB_TARGETS if entry.id == target_id), None)
        if target is None:
            raise ValueError(f"unknown target {target_id!r} for --target")

        input_path = data_path("legacy", filename, subfolder="grib")
        ds = from_source("file", str(input_path)).to_fieldlist()
        out = regrid(ds, out_grid=target.spec, interpolation=interpolation)

        fields = []
        for f in out:
            values = np.asarray(f.to_numpy())
            fields.append({
                "values": values,
                "shape": list(values.shape),
                "grid_type": f.metadata("gridType"),
                "grid_spec": _json_safe(f.geography.grid_spec()),
                "nan": int(np.isnan(values).sum()),
            })
        if not fields:
            raise ValueError("the regridded fieldlist holds no field")
        stacked = np.stack([entry["values"] for entry in fields])

        path = root / _grib_reference_relpath(interpolation, filename, target_id)
        temp = path.with_name(f"{path.stem}.tmp.npz")
        save_values(temp, stacked)
        os.replace(temp, path)

        payload.update(
            status="ok",
            interpolation=interpolation,
            n_fields=len(fields),
            shapes=[entry["shape"] for entry in fields],
            grid_types=[entry["grid_type"] for entry in fields],
            grid_specs=[entry["grid_spec"] for entry in fields],
            nan_counts=[entry["nan"] for entry in fields],
        )
    except Exception as error:
        payload["error_type"] = _exception_name(type(error))
        payload["error"] = str(error) or repr(error)
    payload["seconds"] = round(time.perf_counter() - start, 3)
    print(json.dumps(payload, default=str), flush=True)
    return 0 if payload["status"] == "ok" else 1


# --------------------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------------------


def _build_parser():
    """The command line of the module docstring's Usage section."""
    parser = argparse.ArgumentParser(
        prog="python tests/regrid/regrid_mir/generate_regrid_grib_data.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-dir",
        metavar="PATH",
        help=(
            "root folder for the regrid/ tree and the expected report (relative to the current working "
            f"directory, created if missing); default: ${GENERATE_DIR_ENV} or {GENERATED_DATA_ROOT}"
        ),
    )
    parser.add_argument("--file", action="append", metavar="NAME", help="only this GRIB file (repeatable)")
    parser.add_argument(
        "--interpolation",
        action="append",
        metavar="NAME",
        help="only this interpolation, canonical name or alias (repeatable)",
    )
    parser.add_argument("--target", action="append", metavar="ID", help="only this target (repeatable)")
    parser.add_argument(
        "--timeout", type=float, default=120.0, metavar="SECONDS", help="per-triple timeout (default 120)"
    )
    parser.add_argument("--force", action="store_true", help="overwrite references that already exist")
    parser.add_argument(
        "--report",
        metavar="PATH",
        help=(
            f"where to write the observed expectations JSON (default <output root>/{REPORT_NAME}; "
            '"-" prints it to stdout instead)'
        ),
    )
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    return parser


def main(argv=None):
    """Parse the command line and run the parent or the worker.

    Resolves the output root the way ``generate_regrid_data.py`` does (that function is
    imported here: ``--output-dir``, else ``EARTHKIT_GEO_REGRID_GENERATE_DIR``, else
    ``_dev/test-data``), then sets ``EARTHKIT_GEO_REGRID_TEST_DATA`` to it before anything
    reads a data file, so this process and every worker subprocess resolve new-hierarchy
    files in the generation output folder, never through the download server. The GRIB
    inputs are legacy files and keep being fetched from the server either way.
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
