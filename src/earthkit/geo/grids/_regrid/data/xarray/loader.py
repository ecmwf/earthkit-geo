# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Discovery of geographical variables and their grid in an xarray dataset.

The main entry point is :func:`variables`, which inspects an
``xarray.Dataset``, uses :class:`~.guesser.DefaultCoordinateGuesser` to
classify each variable's coordinates, and returns a :class:`Variable` for
every 2D+ geographical variable, carrying its :class:`~.grid.XarrayGrid`,
grid dimension names, and the corresponding ``eckit.geo.Grid`` used for
regridding.
"""

import itertools
import logging
from typing import Any

LOG = logging.getLogger(__name__)


class Variable:
    """A geographical xarray variable paired with its detected grid.

    Attributes
    ----------
    name : Hashable
        The variable's name.
    variable : xr.DataArray
        The underlying xarray variable.
    geo_dims : Optional[List[str]]
        The names of the variable's geographical grid dimensions.
    xr_grid : Optional[XarrayGrid]
        The grid derived from the variable's xarray coordinates.
    ek_grid : Optional[eckit.geo.Grid]
        The corresponding eckit-geo grid used for regridding.
    """

    def __init__(self, variable):
        """Initialise the wrapper.

        Parameters
        ----------
        variable : xr.DataArray
            The xarray variable to wrap.
        """
        self.name = variable.name
        self.variable = variable
        self.geo_dims = None
        self.xr_grid = None
        self.ek_grid = None


def grid_from_earthkit(ds):
    """Get the eckit-geo grid from an earthkit-data accessor, if available.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to inspect.

    Returns
    -------
    Optional[eckit.geo.Grid]
        The grid, or None if the earthkit accessor or its grid_spec is
        not available.
    """
    if hasattr(ds, "earthkit") and hasattr(ds.earthkit, "grid_spec"):
        gs = ds.earthkit.grid_spec
        if gs is not None:
            from eckit.geo import Grid

            return Grid(gs)
    return None


def grid_from_xr_grid(xr_grid):
    """Build an eckit-geo unstructured lat/lon grid from an :class:`XarrayGrid`.

    Parameters
    ----------
    xr_grid : XarrayGrid
        The grid to convert, providing flat ``latlons``.

    Returns
    -------
    eckit.geo.Grid
        An unstructured lat/lon grid built from the flattened points.
    """
    from eckit.geo import Grid

    lat, lon = xr_grid.latlons
    grid_spec = {"latitudes": lat.flatten().tolist(), "longitudes": lon.flatten().tolist()}
    return Grid(grid_spec)


def adjust_variable_dim_order(variable_dims, coordinates):
    """Reorder a variable's two grid dimensions to ``(lat, lon)`` if given as ``(lon, lat)``.

    Parameters
    ----------
    variable_dims : List[str]
        The two grid dimension names of the variable, in their current order.
    coordinates : List[Coordinate]
        The variable's guessed coordinates, used to identify which
        dimension is latitude and which is longitude.

    Returns
    -------
    List[str]
        The dimension names, reordered to ``(lat, lon)`` if needed.

    Notes
    -----
    This is needed to ensure the order is (lat, lon) when data is passed to
    MIR. When scanning mode support will be added, this function may not be needed.
    """
    if variable_dims and len(variable_dims) == 2:
        order = [None, None]
        for i in range(2):
            for c in coordinates:
                if c.name == variable_dims[i]:
                    if c.is_lat:
                        order[i] = "lat"
                    elif c.is_lon:
                        order[i] = "lon"
                    break

        if order == ["lon", "lat"]:
            variable_dims = [variable_dims[1], variable_dims[0]]

    return variable_dims


def variables(ds, user_ek_grid=None):
    """Find the geographical data variables in a dataset and their grid.

    Skips coordinate/auxiliary variables (referenced via ``coordinates``,
    ``bounds``, ``climatology`` or ``grid_mapping`` attributes) and data
    variables with fewer than two grid coordinates. For each remaining
    variable, guesses its coordinates and grid via
    :class:`~.guesser.DefaultCoordinateGuesser` and resolves the
    corresponding eckit-geo grid.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to inspect.
    user_ek_grid : Any, optional
        A user-provided grid spec or ``eckit.geo.Grid`` used as the input
        grid for every variable instead of guessing it from the dataset.

    Returns
    -------
    List[Variable]
        The geographical variables found, each with its ``ek_grid``,
        ``xr_grid`` and ``geo_dims`` populated.
    """
    ek_grid = None

    if user_ek_grid is not None:
        from eckit.geo import Grid

        ek_grid = Grid(user_ek_grid)

    from .guesser import DefaultCoordinateGuesser

    guess = DefaultCoordinateGuesser(ds)

    skip = set()

    def _skip_attr(v: Any, attr_name: str) -> None:
        attr_val: str = getattr(v, attr_name, "")
        if isinstance(attr_val, str):
            v = attr_val.split()
            if v:
                skip.update(v)

    for name in itertools.chain(ds.coords, ds.data_vars):
        variable = ds[name]
        _skip_attr(variable, "coordinates")
        _skip_attr(variable, "bounds")
        _skip_attr(variable, "climatology")
        _skip_attr(variable, "grid_mapping")

    result = []

    # Select only geographical variables
    for name in ds.data_vars:
        if name in skip:
            continue

        variable = ds[name]
        coordinates = []

        for coord in variable.coords:
            c = guess.guess(ds[coord], coord)
            assert c, f"Could not guess coordinate for {coord}"
            if coord not in variable.dims:
                LOG.debug("%s: coord=%s (not a dimension): dims=%s", variable, coord, variable.dims)
                c.is_dim = False
            coordinates.append(c)

        grid_coords: int = sum(1 for c in coordinates if c.is_grid)

        if grid_coords < 2:
            LOG.debug("Skipping %s (not 2D): %s", variable, [(c, c.is_grid, c.is_dim) for c in coordinates])
            continue

        xr_grid = guess.grid(coordinates, variable)

        variable_dims = xr_grid.variable_dims
        variable_dims = adjust_variable_dim_order(variable_dims, coordinates)

        def _check_values_geo() -> None:
            """Handle the case where lat, lon is a coordinate but not a dimension and their
            dimension is not recognised as a grid coordinate.

            E.g.:

            Dimensions:  (level: 2, values: 9)
            Coordinates:
                * level    (level) int64 16B 700 500
                    lat      (values) int64 72B 50 50 50 40 40 40 30 30 30
                    lon      (values) int64 72B 0 10 20 0 10 20 0 10 20
                    values   (values) int64 72B 0 1 2 3 4 5 6 7 8
                Dimensions without coordinates: values
                Data variables:
                    a        (level, values) int64 144B 11 12 13 21 22 23 ... 24 25 26 34 35 36

            """
            from .coordinates import UnsupportedCoordinate

            g = [c for c in coordinates if c.is_grid]

            for c in coordinates:
                if c.is_dim and isinstance(c, UnsupportedCoordinate):
                    for cx in g:
                        if c.name in cx.variable.sizes:
                            c.is_grid = True
                            break

        _check_values_geo()

        ek_grid = None
        if user_ek_grid is not None:
            if not isinstance(user_ek_grid, Grid):
                ek_grid = Grid(user_ek_grid)
            else:
                ek_grid = user_ek_grid
        else:
            ek_grid = grid_from_earthkit(ds)
            if ek_grid is None:
                ek_grid = grid_from_xr_grid(xr_grid)

        v = Variable(variable)
        v.ek_grid = ek_grid
        v.xr_grid = xr_grid
        v.geo_dims = variable_dims

        result.append(v)

    return result
