# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import pytest

pytest.importorskip("earthkit.cli.standard_args")
pytest.importorskip("earthkit.data")

from click.testing import CliRunner  # noqa: E402
from earthkit.cli.main import earthkit  # noqa: E402
from earthkit.cli.standard_args import SOURCE_HELP, TARGET_HELP  # noqa: E402

from earthkit.cli.geo import grids as grids_cli  # noqa: E402
from earthkit.geo.utils.testing import NO_MIR, earthkit_test_data_path  # noqa: E402

O32_GRIB = earthkit_test_data_path("o32.grib2")


def _invoke(*args, exit_code=0):
    result = CliRunner().invoke(earthkit, [str(a) for a in args])
    assert result.exit_code == exit_code, result.output + repr(result.exception)
    return result


def test_cli_registers_regrid():
    assert earthkit.get_command(None, "regrid") is grids_cli.regrid


def test_cli_info_lists_geo_commands():
    output = _invoke("info").output
    assert "earthkit-geo" in output
    assert "(commands: regrid)" in output


def test_cli_regrid_help():
    output = _invoke("regrid", "--help").output
    assert "[OPTIONS] SOURCE TARGET\n" in output
    for option in ("-g, --target-grid-spec", "--interpolation"):
        assert option in output
    assert "-i," not in output
    for text in ("--source", "--target ", "--profile"):
        assert text not in output
    # The shared descriptions from earthkit-utils, rewrapped by click
    usage = " ".join(output.split())
    for text in (SOURCE_HELP, TARGET_HELP):
        assert " ".join(text.split()) in usage


@pytest.mark.parametrize(
    "value, expected",
    (
        ("O96", "O96"),
        (" 5/5 ", "5/5"),
        ('{"grid": [1, 1]}', {"grid": [1, 1]}),
        ('"O96"', "O96"),
        ({"grid": [1, 1]}, {"grid": [1, 1]}),
    ),
)
def test_grid_spec(value, expected):
    assert grids_cli.GRID_SPEC.convert(value, None, None) == expected


@pytest.mark.parametrize(
    "value, message",
    (
        ("", "must not be empty"),
        ("{", "is not valid JSON"),
        ("[1, 1]", "must be a JSON object or a grid name"),
    ),
)
def test_grid_spec_invalid(value, message):
    result = _invoke("regrid", O32_GRIB, "out.grib", "-g", value, exit_code=2)
    assert message in result.output


@pytest.mark.parametrize(
    "options, interpolation",
    (
        ([], "linear"),
        (["--interpolation", "nearest-neighbour"], "nearest-neighbour"),
    ),
)
def test_cli_regrid_calls_regrid(tmp_path, monkeypatch, options, interpolation):
    import earthkit.data as ekd

    import earthkit.geo as ekg

    calls = []

    def _regrid(data, **kwargs):
        calls.append((type(data).__name__, kwargs))
        return "regridded"

    monkeypatch.setattr(ekg, "regrid", _regrid)
    monkeypatch.setattr(ekd, "to_target", lambda *args, **kwargs: calls.append((args, kwargs)))
    out_path = tmp_path / "out.grib"
    _invoke("regrid", O32_GRIB, out_path, "-g", '{"grid": [5, 5]}', *options)
    assert calls[0][1] == {"out_grid": {"grid": [5, 5]}, "interpolation": interpolation}
    assert calls[1] == (("file", str(out_path)), {"data": "regridded"})


def test_cli_regrid_missing_input(tmp_path):
    result = _invoke("regrid", tmp_path / "missing.grib", tmp_path / "out.grib", "-g", "O96", exit_code=2)
    assert "Invalid value for 'SOURCE'" in result.output


@pytest.mark.skipif(NO_MIR, reason="No mir available")
def test_cli_regrid(tmp_path):
    import earthkit.data as ekd

    out_path = tmp_path / "out.grib"
    _invoke("regrid", O32_GRIB, out_path, "-g", '{"grid": [5, 5]}')
    result = ekd.from_source("file", str(out_path)).to_fieldlist()
    assert len(result) == len(ekd.from_source("file", O32_GRIB).to_fieldlist())
    assert result[0].shape == (37, 72)
