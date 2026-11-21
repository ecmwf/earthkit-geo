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


# The 3 x 3 grid of the tests below, in the plane of a LAEA projection centred on (52 N, 10 E).
# PROJ_SOURCE_SPEC spells its projection the way the xarray loader detects it from a CF
# proj4_params attribute.
PROJ = "+proj=laea +lat_0=52 +lon_0=10 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
PROJ_SOURCE_SPEC = {
    "type": "regular_xy",
    "x": [-50000, 50000, 50000],
    "y": [50000, -50000, -50000],
    "projection": {"type": "proj", "source": PROJ},
}


def test_mir_backend_grid_spec_passthrough():
    """The MIR backend returns a grid spec as it is, projected gridspecs included.

    Converting a PROJ-source projected input grid into its unstructured lat/lon points
    is the job of the data handlers (see
    :func:`earthkit.geo.grids._regrid.data.utils.scatter_projected_grid`), not of the
    backend: :meth:`~earthkit.geo.grids._regrid.backends.mir.MirBackend.get_grid_spec`
    only extracts the spec from a ``Grid`` and passes a plain spec through unchanged.
    """
    from eckit.geo import Grid

    from earthkit.geo.grids._regrid.backends.mir import MirBackend

    projected = Grid(PROJ_SOURCE_SPEC)
    assert MirBackend.get_grid_spec(projected) == projected.spec
    assert MirBackend.get_grid_spec(PROJ_SOURCE_SPEC) == PROJ_SOURCE_SPEC


@requires_mir
def test_regrid_array_projected_input_to_ll():
    """A PROJ-source projected array regrids onto a lat/lon grid, its points where they are."""
    # the 3 x 3 input, its values numbered in the row-major (y, x) order of the projected grid
    values = np.arange(9, dtype=np.float64)
    out_grid = {"grid": [1, 1], "area": [53, 9, 51, 11]}

    result, result_spec = regrid_array(values, PROJ_SOURCE_SPEC, out_grid, interpolation="nearest-neighbour")

    assert result.shape == (3, 3)
    assert result_spec.items() >= out_grid.items()
    # the output points 53/52/51 N by 9/10/11 E, each nearest to the input point of the same
    # corner of the grid: this only holds when MIR was given the input points in the order
    # of the values, and where the projection says they are
    assert np.array_equal(result.reshape(-1), values)

    # linear interpolation over the three middle longitudes keeps the centre value
    result, _ = regrid_array(values, PROJ_SOURCE_SPEC, out_grid, interpolation="linear")
    assert np.isclose(result[1, 1], 4.0)
