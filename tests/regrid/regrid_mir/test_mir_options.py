# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""The temporary ``area`` / ``rotation`` gridspec workaround of the MIR backend.

``MirBackend.adjust_options`` moves these keys from the output grid spec into the MIR
job options. These tests need no ``mir`` module: they only check the dict handling.
"""

import pytest

from earthkit.geo.grids._regrid.backends.mir import MirBackend


@pytest.mark.parametrize(
    "grid, expected_grid, expected_kwargs",
    [
        ({"grid": [1, 1], "area": [10, -10, -10, 10]}, {"grid": [1, 1]}, {"area": "10/-10/-10/10"}),
        ({"grid": [1, 1], "rotation": [-43, 10]}, {"grid": [1, 1]}, {"rotation": "-43/10"}),
        (
            {"grid": [1, 1], "area": "10/-10/-10/10", "rotation": "-43/10"},
            {"grid": [1, 1]},
            {"area": "10/-10/-10/10", "rotation": "-43/10"},
        ),
    ],
)
def test_mir_adjust_options_moves_keys(grid, expected_grid, expected_kwargs):
    grid_in = dict(grid)
    kwargs_in = {"interpolation": "linear"}

    with pytest.warns(DeprecationWarning):
        new_grid, new_kwargs = MirBackend.adjust_options(grid_in, kwargs_in)

    assert new_grid == expected_grid
    assert new_kwargs == {"interpolation": "linear", **expected_kwargs}
    # the inputs are not modified
    assert grid_in == grid
    assert kwargs_in == {"interpolation": "linear"}


def test_mir_adjust_options_no_keys():
    grid = {"grid": [1, 1]}
    kwargs = {"interpolation": "linear"}
    new_grid, new_kwargs = MirBackend.adjust_options(grid, kwargs)
    assert new_grid is grid
    assert new_kwargs is kwargs


@pytest.mark.parametrize("rotation", [[-43], [-43, 10, 0], -43, None])
def test_mir_adjust_options_bad_rotation(rotation):
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError, match="Invalid rotation"):
        MirBackend.adjust_options({"grid": [1, 1], "rotation": rotation}, {})
