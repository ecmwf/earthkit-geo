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
from typing import Hashable, List

from earthkit.geo.grids._regrid.data.utils import create_grid_object, create_unstructured_grid_object
from earthkit.geo.utils import ensure_list

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

        self.point_order = None
        """The permutation taking the values to the point order of :attr:`ek_grid` (Optional[np.ndarray]).

        Set only for a variable on a CF reduced Gaussian sub-area whose grid was
        detected from the dataset itself, and whose points eckit iterates in
        another order than the file holds them: ``values[point_order]`` are the
        values in the order of the eckit grid (see
        :attr:`~.grid.ReducedGaussianGrid.point_order`). ``None``, the default,
        whenever the values already sit in the grid's order — which is always
        the case when the grid was not detected but given explicitly, by the
        user's ``in_grid`` or by an earthkit-data accessor: an explicitly given
        grid takes the data exactly as it stands (decision A8).
        """


def auxiliary_names(ds) -> set:
    """The names of a dataset's variables which another variable declares as auxiliary.

    The ``coordinates``, ``bounds``, ``climatology`` and ``grid_mapping`` attributes
    name, whitespace-separated, the variables which describe the referring one:
    its auxiliary coordinates, its cell bounds, the variable its climatology
    lives in, the mapping which locates its points. Those describe data, they
    are not data of their own: the loader never regrids a variable it finds
    under one of these names (see :func:`variables`), whichever dimensions
    happen to index it.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to inspect.

    Returns
    -------
    Set[Hashable]
        Every name referenced by one of the four attributes, in any order.
    """
    names = set()

    for variable in itertools.chain(ds.coords.values(), ds.data_vars.values()):
        for attr_name in ("coordinates", "bounds", "climatology", "grid_mapping"):
            attr_val = getattr(variable, attr_name, "")
            if isinstance(attr_val, str):
                names.update(attr_val.split())

    return names


def reduced_gaussian_indexed(ds) -> List[Hashable]:
    """The names of a dataset's data variables indexed by a CF reduced Gaussian point index.

    A variable is named when one of its 1D coordinates carries the CF
    ``reduced_gaussian_index`` standard name and indexes one of its dimensions,
    whichever ``grid_mapping`` the variable declares and whether or not the
    dataset holds the mapping and points-per-latitude tables that would place
    the numbered points. Variables another variable declares auxiliary (see
    :func:`auxiliary_names`) are not named: the loader skips them as data on
    purpose, so the clear error which reports a field the grid tables cannot
    place is not their business, tables or none.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to inspect.

    Returns
    -------
    List[Hashable]
        The names of the data variables so indexed, in dataset order.
    """
    skip = auxiliary_names(ds)
    names = []

    for name in ds.data_vars:
        if name in skip:
            continue
        variable = ds[name]
        for coord in variable.coords:
            c = variable.coords[coord]
            if c.ndim != 1 or c.dims[0] not in variable.dims:
                continue
            if str(c.attrs.get("standard_name", "")).strip().lower() != "reduced_gaussian_index":
                continue
            names.append(name)
            break

    return names


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


def variables(ds, user_ek_grid=None, user_in_dims=None):
    """Find the geographical data variables in a dataset and their grid.

    Skips coordinate/auxiliary variables (referenced via ``coordinates``,
    ``bounds``, ``climatology`` or ``grid_mapping`` attributes) and data
    variables with fewer than two grid coordinates. A variable is only kept
    without two grid coordinates when the dataset describes its points with a
    CF reduced Gaussian grid mapping instead of lat/lon coordinates, or when
    an input grid is given explicitly together with the dimensions the values
    are laid out along (see ``user_in_dims``); a variable which carries both
    the mapping and auxiliary coordinates is located by the mapping, which
    says more about the grid than auxiliary coordinates of its own (see
    :meth:`~.guesser.CoordinateGuesser.reduced_gaussian_grid`).
    For each remaining variable, guesses its coordinates and grid via
    :class:`~.guesser.DefaultCoordinateGuesser` and resolves the
    corresponding eckit-geo grid.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to inspect.
    user_ek_grid : Any, optional
        A user-provided grid spec or ``eckit.geo.Grid`` used as the input
        grid for every variable instead of guessing it from the dataset.
    user_in_dims : List[str], optional
        The user-provided grid dimension names (the ``in_dims`` of a regrid
        call). A variable which is neither located by the mapping of a CF
        reduced Gaussian grid nor by two grid coordinates is kept as
        index-only data when a grid is given and all of these dimensions
        index it: its values then number the points of the given grid in the
        order the named dimensions hold them, and are regridded without any
        reordering (decision A8). The sizes of the dimensions are not
        checked here.

    Returns
    -------
    List[Variable]
        The geographical variables found, each with its ``ek_grid``,
        ``xr_grid``, ``geo_dims`` and ``point_order`` populated.
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

    # the dimensions a user has named as the grid dimensions of the data
    # (the "in_dims" of a regrid call), which index-only data of an explicitly
    # given grid can be laid out along (see the fallback in the loop below)
    explicit_in_dims = ensure_list(user_in_dims) if user_in_dims is not None else []

    # the variables the dataset declares as auxiliary: they describe other data,
    # and are never regridded as data of their own (see :func:`auxiliary_names`)
    skip = auxiliary_names(ds)

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

        # the points of a CF reduced Gaussian grid are numbered by a grid mapping
        # rather than described by lat/lon coordinates, and a variable which
        # carries both is located by the mapping: it is the stronger statement
        # about the grid the data sits on. Only when the mapping describes no
        # grid does the variable fall back to the coordinates, and one with
        # fewer than two of them is then kept only as index-only data of an
        # explicitly given grid, and skipped as it was before otherwise
        xr_grid = guess.reduced_gaussian_grid(variable)
        if xr_grid is None:
            if grid_coords < 2:
                if ek_grid is not None and explicit_in_dims and all(d in variable.dims for d in explicit_in_dims):
                    # the data numbers the points of a grid it does not describe:
                    # regriddable when the caller gives that grid with "in_grid"
                    # and names its dimensions with "in_dims" (decision A8). The
                    # grid dimensions are the user's, so the regridder reads them
                    # off "in_dims" and checks the sizes; the values are taken as
                    # they stand, with no reordering.
                    v = Variable(variable)
                    v.ek_grid = ek_grid
                    result.append(v)
                    continue
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
        # only a grid detected from the dataset's own mapping knows the order
        # the file holds its points in; an explicitly given grid takes the
        # data as it stands (decision A8)
        if ek_grid is None:
            v.point_order = xr_grid.point_order

        result.append(v)

    return result
