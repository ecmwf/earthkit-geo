# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


def get_input(backend, in_grid, out_grid):
    # for precomputed backend we need to build a special gridspec object
    # to match the matrix inventory items
    # TODO: remove this limitation
    if backend.name == "precomputed":
        from earthkit.geo.grids._regrid.backends.precomputed.gridspec import _GridWrapper

        in_grid = _GridWrapper.from_any(in_grid)
        out_grid = _GridWrapper.from_any(out_grid)
    else:
        from eckit.geo import Grid

        if in_grid is not None and not isinstance(in_grid, Grid):
            in_grid = Grid(in_grid)
        if out_grid is not None and not isinstance(out_grid, Grid):
            out_grid = Grid(out_grid)

    return in_grid, out_grid
