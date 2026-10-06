# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Commands contributed by earthkit-geo to the shared ``earthkit`` command line interface.
The ``earthkit`` console script itself lives in :mod:`earthkit.utils.cli`. The commands
defined here are registered with it through the ``earthkit.cli`` entry point group in
``pyproject.toml``, so ``earthkit regrid <...>`` becomes available once earthkit-geo is installed.
"""

import json

import click


class GridSpecParamType(click.ParamType):
    """Click parameter type for a grid spec given as JSON or as a plain grid name."""

    name = "grid_spec"

    def convert(self, value, param, ctx):
        if isinstance(value, dict):
            return value

        value = value.strip()
        if not value:
            self.fail("grid spec must not be empty", param, ctx)

        # anything that is not JSON is taken as a grid name, e.g. "O96" or "5/5"
        if value[0] not in "{[\"":
            return value

        try:
            spec = json.loads(value)
        except json.JSONDecodeError as e:
            self.fail(f"{value!r} is not valid JSON: {e}", param, ctx)

        if not isinstance(spec, (dict, str)):
            self.fail(
                f"{value!r} must be a JSON object or a grid name, got {type(spec).__name__}",
                param,
                ctx,
            )

        return spec


GRID_SPEC = GridSpecParamType()


@click.command()
@click.argument("source-file", type=click.Path(exists=True, dir_okay=False))
@click.argument("target-file", type=click.Path(exists=False, dir_okay=False))
@click.option(
    "-g",
    "--target-grid-spec",
    type=GRID_SPEC,
    required=True,
    help='Target grid specification, either as JSON, e.g. \'{"grid": [1, 1]}\', '
    "or as a grid name, e.g. O96.",
)
@click.option(
    "-i",
    "--interpolation",
    type=str,
    required=False,
    default="linear",
    help='Interpolation method, e.g. "nearest-neighbour".',
)
def regrid(source_file, target_file, target_grid_spec, interpolation):
    """Regrids the input data to a new grid."""

    import earthkit.data as ekd
    import earthkit.geo as ekg

    in_data = ekd.from_source("file", source_file).to_fieldlist()

    out_data = ekg.regrid(in_data, out_grid=target_grid_spec, interpolation=interpolation)

    out_data.to_target("file", target_file)


COMMANDS = {
    "regrid": regrid,
}