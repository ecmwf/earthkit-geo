# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


import pytest

from earthkit.geo import Grid


@pytest.mark.parametrize(
    "spec, shape",
    [
        (dict(grid=[1, 1]), (181, 360)),
    ],
)
def test_grid(spec, shape):
    grid = Grid(spec)
    assert grid.shape == shape
