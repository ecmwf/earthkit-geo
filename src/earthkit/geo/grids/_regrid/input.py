# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


def get_input(backend, in_grid, out_grid):
    """Normalise ``in_grid``/``out_grid`` into the objects expected by ``backend``.

    Parameters
    ----------
    backend : Backend
        The regrid backend the grids will be used with.
    in_grid : dict, str, eckit.geo.Grid or None
        The input grid spec.
    out_grid : dict, str, eckit.geo.Grid or None
        The output grid spec.

    Returns
    -------
    Tuple[Any, Any]
        The normalised ``(in_grid, out_grid)``: a pair of
        :class:`~earthkit.geo.grids._regrid.backends.precomputed.gridspec._GridWrapper`
        when ``backend.name == "precomputed"`` (to match the matrix
        inventory items), otherwise a pair of ``eckit.geo.Grid`` objects
        (left as None where the input was None).
    """
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
