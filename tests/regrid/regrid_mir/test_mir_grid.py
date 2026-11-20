# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Smoke test that any gridspec accepted by ``mir.Grid`` can regrid to another one.

The interpolation names and the ``mir`` skip mark come from ``regrid_common``. Each
case loops over all interpolations internally and only checks the shape of the
result (no data files). Exact duplicate parameter tuples of the old list
({"grid": [5, 5]} -> {"grid": [10, 10]} appeared four times, {"grid": "O32"} and
{"grid": "N32"} -> {"grid": [10, 10]} twice each) are kept once.
"""

import numpy as np
import pytest
from regrid_common import INTERPOLATION_NAMES, requires_mir

from earthkit.geo.grids.array import regrid as regrid_array


@requires_mir
@pytest.mark.parametrize(
    "gs_in, gs_out",
    [
        ({"grid": [5, 5]}, {"grid": [10, 10]}),
        ({"grid": [5, 5], "area": [90, 0, -90, 360]}, {"grid": [10, 10]}),
        ({"grid": [5, 5], "area": [90, 0, -90, 355]}, {"grid": [10, 10]}),
        ({"grid": [5, 5]}, {"grid": [10, 10], "area": [90, 0, -90, 360]}),
        ({"grid": [5, 5]}, {"grid": [10, 10], "area": [90, 0, -90, 350]}),
        (
            {"grid": [5, 5], "area": [90, 0, -90, 360]},
            {"grid": [10, 10], "area": [90, 0, -90, 360]},
        ),
        ({"grid": "O32"}, {"grid": [10, 10]}),
        ({"grid": "o32"}, {"grid": [10, 10]}),
        ({"grid": "O32", "area": [90, 0, -90, 360]}, {"grid": [10, 10]}),
        ({"grid": "O32", "area": [87.8638, 0, -87.8638, 357.5]}, {"grid": [10, 10]}),
        ({"grid": "O32", "global": 1}, {"grid": [10, 10]}),
        (
            {"grid": "O32", "global": 1, "area": [87.8638, 0, -87.8638, 357.5]},
            {"grid": [10, 10]},
        ),
        ({"grid": "N32"}, {"grid": [10, 10]}),
        ({"grid": "N32", "area": [90, 0, -90, 360]}, {"grid": [10, 10]}),
        ({"grid": "N32", "area": [87.8638, 0, -87.8638, 357.188]}, {"grid": [10, 10]}),
        ({"grid": "N32", "area": [87.8638, 0, -87.8638, 357.189]}, {"grid": [10, 10]}),
        (
            {"grid": "O32", "area": [87.8638, 0, -87.8638, 357.6]},
            {"grid": [10, 10]},
        ),
        ({"grid": "N32", "global": 1}, {"grid": [10, 10]}),
        ({"grid": "n32"}, {"grid": [10, 10]}),
        (
            {"grid": "N32", "global": 1, "area": [87.8638, 0, -87.8638, 357.188]},
            {"grid": [10, 10]},
        ),
        ({"grid": "H4"}, {"grid": [10, 10]}),
        ({"grid": "H4", "order": "ring"}, {"grid": [10, 10]}),
        ({"grid": "H4", "order": "nested"}, {"grid": [10, 10]}),
        ({"grid": "h4", "order": "nested"}, {"grid": [10, 10]}),
        # Disabled, as before the refactor: the O32/N32 areas starting at 0.01
        # ({"grid": "O32", "area": [87.8638, 0.01, -87.8638, 357.5]} and the N32
        # equivalent), and the slow ({"grid": "eORCA025_T"}, {"grid": "O96"}).
    ],
)
def test_regrid_with_mir_gridspec(gs_in, gs_out):
    from mir import Grid

    in_grid = Grid(**gs_in)
    values = np.random.rand(in_grid.size())

    for interpolation in INTERPOLATION_NAMES:
        if gs_in["grid"] == "eORCA025_T" and interpolation == "grid-box-average":
            continue

        # TODO: make this code work
        values_res, _ = regrid_array(values, gs_in, gs_out, interpolation=interpolation)

        result_grid = Grid(gs_out)  # NOTE: not necessarily true
        assert values_res.shape == result_grid.shape
