# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Grid commands of earthkit-geo for the shared ``earthkit`` command line interface.

Provides ``earthkit regrid``, which wraps :func:`earthkit.geo.regrid`.
"""

import json

import click
from earthkit.cli.main import earthkit
from earthkit.cli.standard_args import SOURCE_HELP, TARGET_HELP, add_options, source_options, target_options


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
        if value[0] not in '{["':
            return value

        try:
            spec = json.loads(value)
        except json.JSONDecodeError as e:
            self.fail(f"{value!r} is not valid JSON: {e}", param, ctx)

        if not isinstance(spec, (dict, str)):
            self.fail(f"{value!r} must be a JSON object or a grid name, got {type(spec).__name__}", param, ctx)

        return spec


GRID_SPEC = GridSpecParamType()


@earthkit.command(
    help=f"""Regrid SOURCE to a new grid and write the result to TARGET.

SOURCE: {SOURCE_HELP}

TARGET: {TARGET_HELP}

\b
Example:
    earthkit regrid input.grib output.grib --target-grid-spec O96
    earthkit regrid input.grib output.grib --target-grid-spec '{{"grid": [1, 1]}}' --interpolation nearest-neighbour
"""
)
@add_options([source_options(positional=True), target_options(positional=True)])
@click.option(
    "--target-grid-spec",
    type=GRID_SPEC,
    required=True,
    help="Target grid specification, either as JSON, e.g. '{\"grid\": [1, 1]}', or as a grid name, e.g. O96.",
)
@click.option(
    "--interpolation",
    default="linear",
    show_default=True,
    help="Interpolation method, e.g. 'linear', 'nearest-neighbour' or 'grid-box-average'.",
)
def regrid(source, target, target_grid_spec, interpolation):
    import earthkit.geo as ekg

    target.to_target(ekg.regrid(source.to_fieldlist(), out_grid=target_grid_spec, interpolation=interpolation))
