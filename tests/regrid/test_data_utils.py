# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Tests for the regrid data utils.

Ported from ``tests/regrid/regrid_mir/test_mir_grid.py``: the conversion of a
PROJ-source projected input grid into the unstructured lat/lon grid of its points
now lives in :func:`earthkit.geo.grids._regrid.data.utils.scatter_projected_grid`,
which returns the ``Grid`` object itself instead of the uid dict the MIR backend
wraps it in.
"""

import numpy as np
from eckit.geo import Grid

from earthkit.geo.grids._regrid.data.utils import scatter_projected_grid

# The 3 x 3 grid of the tests below, in the plane of a LAEA projection centred on (52 N, 10 E).
# PROJ_SOURCE_SPEC spells its projection the way the xarray loader detects it from a CF
# proj4_params attribute, and GRID_SPEC the way eckit reads a projection from a structured
# spec: both build the same points, in the same order.
PROJ = "+proj=laea +lat_0=52 +lon_0=10 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
PROJ_SOURCE_SPEC = {
    "type": "regular_xy",
    "x": [-50000, 50000, 50000],
    "y": [50000, -50000, -50000],
    "projection": {"type": "proj", "source": PROJ},
}
GRID_SPEC = {
    "type": "regular_xy",
    "grid": [50000, 50000],
    "bounding_box_xy": [-50000, -50000, 50000, 50000],
    "projection": {"type": "laea", "lat_0": 52, "lon_0": 10, "reverse": True, "figure": "grs80"},
}


def test_scatter_projected_input_grid():
    """A PROJ-source projected grid is scattered to the unstructured lat/lon grid of its points."""
    projected = Grid(PROJ_SOURCE_SPEC)
    scattered = projected.to_unstructured_ll()

    # the grid, as a Grid object or as a plain spec, with or without its "type"
    for grid in (projected, PROJ_SOURCE_SPEC, projected.spec):
        result = scatter_projected_grid(grid)
        assert isinstance(result, Grid)
        assert result is not grid
        assert result.type == "unstructured_ll"
        assert result.uid == scattered.uid

        # the same points, in the same order: the projected grid's x axis scanned fastest,
        # which is the order data laid out on the grid is in
        lat, lon = (np.asarray(v) for v in result.to_latlons())
        expected_lat, expected_lon = (np.asarray(v) for v in projected.to_latlons())
        assert lat.shape == expected_lat.shape == (projected.size(),)
        assert np.allclose(lat, expected_lat)
        assert np.allclose(lon, expected_lon)


def test_scatter_projected_grid_keeps_other_grids():
    """Grids that are not PROJ-source projected are returned unchanged, as the same object."""
    # a projected grid whose projection is a structured spec, and an unprojected one
    named = Grid(GRID_SPEC)
    assert scatter_projected_grid(named) is named
    plain = Grid({"x": [0, 10, 10], "y": [10, 0, -10]})
    assert scatter_projected_grid(plain) is plain

    # the grids that have always been taken as they are
    ll = Grid({"grid": [10, 10]})
    assert scatter_projected_grid(ll) is ll
    unstructured = Grid({"latitudes": [0.0, 10.5, -20.25], "longitudes": [0.0, 90.0, 180.0]})
    assert scatter_projected_grid(unstructured) is unstructured


def test_scatter_projected_grid_keeps_other_specs():
    """Specs that are not of the PROJ kind, JSON strings and None are passed through."""
    assert scatter_projected_grid(GRID_SPEC) is GRID_SPEC
    assert scatter_projected_grid('{"grid": [10, 10]}') == '{"grid": [10, 10]}'
    assert scatter_projected_grid(None) is None


def test_scatter_projected_grid_keeps_other_proj_specs():
    """A dict spec with a PROJ projection that is not a ``regular_xy`` grid is returned as the same dict."""
    spec = {"grid": [10, 10], "projection": {"type": "proj", "source": "+proj=longlat +datum=WGS84 +no_defs"}}
    assert scatter_projected_grid(spec) is spec
