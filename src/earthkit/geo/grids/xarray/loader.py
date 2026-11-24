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

from earthkit.geo.grids._regrid.data.utils import create_grid_object, create_unstructured_grid_object

from .patch import patch_dataset

LOG = logging.getLogger(__name__)


class Variable:
    """A geographical xarray variable paired with its detected grid."""

    def __init__(self, variable):
        """Initialise the wrapper.

        Parameters
        ----------
        variable : xr.DataArray
            The xarray variable to wrap.
        """
        self.name = variable.name
        """The variable's name (Hashable)."""

        self.variable = variable
        """The underlying xarray variable (``xr.DataArray``)."""

        self.geo_dims = None
        """The names of the variable's geographical grid dimensions (Optional[List[str]])."""

        self.xr_grid = None
        """The grid derived from the variable's xarray coordinates (Optional[XarrayGrid])."""

        self.ek_grid = None
        """The corresponding eckit-geo grid used for regridding (Optional[eckit.geo.Grid])."""


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
        return create_grid_object(gs)
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
    gs = xr_grid.grid_spec
    if gs is not None:
        return create_grid_object(gs)
    lat, lon = xr_grid.latlons
    return create_unstructured_grid_object(lat, lon)


# the eckit-geo grid types that are scanned axis by axis; whether they iterate
# their points with the x axis fastest is decided by their scan order, see
# :func:`_scans_x_fastest`
_STRUCTURED_GRID_TYPES = ("regular_xy", "regular_ll")

# the eckit default scan order, used by a gridspec that does not name one
_DEFAULT_ORDER = "i+j-"


def _scans_x_fastest(grid) -> bool:
    """Whether the points of an eckit-geo grid iterate with its x axis fastest.

    Parameters
    ----------
    grid : eckit.geo.Grid
        The grid to inspect.

    Returns
    -------
    bool
        True for a ``regular_xy`` grid and for a ``regular_ll`` grid whose
        scan order leaves the x axis (the longitude of a ``regular_ll``)
        fastest, False for any other grid type, for an unrecognised grid and
        for the orders that iterate the grid with its x axis outermost.

    Notes
    -----
    eckit spells a scan order ``<outer><dir><inner><dir>``, naming first the
    axis whose points are furthest apart (the rows, ``i``: y for a
    ``regular_xy``, the latitude for a ``regular_ll``) and then the axis whose
    points are consecutive (the columns, ``j``); an order without a name is
    :data:`_DEFAULT_ORDER`, which is i-outermost. The points therefore sit x
    fastest exactly when the order starts with ``i``: only then does data
    laid out (x, y) have to be transposed to (y, x) before it is flattened
    onto the grid (see :func:`adjust_variable_dim_order`).
    """
    if getattr(grid, "type", None) not in _STRUCTURED_GRID_TYPES:
        return False
    order = getattr(grid, "order", None) or grid.spec.get("order") or _DEFAULT_ORDER
    return order.startswith("i")


def adjust_variable_dim_order(variable_dims, coordinates, swap_xy=False):
    """Reorder a variable's two grid dimensions to ``(lat, lon)``, or ``(y, x)``, if needed.

    Parameters
    ----------
    variable_dims : List[str]
        The two grid dimension names of the variable, in their current order.
    coordinates : List[Coordinate]
        The variable's guessed coordinates, used to identify which
        dimension is latitude and which is longitude, or x and y.
    swap_xy : bool
        Whether the variable's input grid is an explicitly given eckit grid
        that iterates its points with the x axis fastest (see
        :func:`_scans_x_fastest`), given by the user or by an earthkit-data
        accessor. Only then are ``(x, y)`` dimensions swapped to ``(y, x)``:
        an unstructured explicit grid and a structured one scanned with the
        x axis outermost hold their points in the order the data was laid
        out, and a detected grid is built to match that order in the first
        place.

    Returns
    -------
    List[str]
        The dimension names, reordered to ``(lat, lon)``, and to ``(y, x)``
        when ``swap_xy``, if needed.

    Notes
    -----
    This is needed to ensure the order is (lat, lon), and (y, x) for a
    projected grid, when data is passed to MIR. When scanning mode support
    will be added, this function may not be needed.
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
                    elif c.is_y:
                        order[i] = "y"
                    elif c.is_x:
                        order[i] = "x"
                    break

        if order == ["lon", "lat"]:
            variable_dims = [variable_dims[1], variable_dims[0]]
        elif swap_xy and order == ["x", "y"]:
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
    patch = {"latlon": None}
    ds = patch_dataset(ds, patch)

    from .guesser import DefaultCoordinateGuesser

    guess = DefaultCoordinateGuesser(ds)

    # the input grid, when it is given explicitly rather than detected from the
    # dataset: by the user, or by an earthkit-data accessor. When it iterates
    # its points with the x axis fastest, the dimensions of data laid out
    # (x, y) are reordered to (y, x) to match it (below); a grid scanned the
    # other way round, like an unstructured one, takes the data in its own
    # dimension order, and so does a detected grid, which is built to match it
    ek_grid = create_grid_object(user_ek_grid)
    if ek_grid is None:
        ek_grid = grid_from_earthkit(ds)
    swap_xy = _scans_x_fastest(ek_grid)

    skip = set()

    def _skip_attr(v: Any, attr_name: str) -> None:
        """Add the whitespace-separated names in ``v``'s ``attr_name`` attribute to ``skip``.

        Parameters
        ----------
        v : Any
            The xarray variable to inspect.
        attr_name : str
            The name of the attribute to read (e.g. ``"coordinates"``,
            ``"bounds"``).
        """
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
        variable_dims = adjust_variable_dim_order(variable_dims, coordinates, swap_xy=swap_xy)

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

        variable_ek_grid = ek_grid
        if variable_ek_grid is None:
            variable_ek_grid = grid_from_xr_grid(xr_grid)

        v = Variable(variable)
        v.ek_grid = variable_ek_grid
        v.xr_grid = xr_grid
        v.geo_dims = variable_dims

        result.append(v)

    return result
