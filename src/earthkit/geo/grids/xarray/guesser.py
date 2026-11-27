# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Coordinate and grid guessing for xarray datasets.

Defines :class:`CoordinateGuesser`, an abstract base that inspects an
``xarray.Dataset``'s coordinate variables (using CF-convention attributes
such as ``standard_name``, ``axis``, ``long_name`` and ``units``, plus
common naming fallbacks) and classifies each one into a
:class:`~.coordinates.Coordinate` subclass (latitude, longitude, x, y,
point, scalar, or unsupported). It also assembles the resulting
coordinates into a :class:`~.grid.Grid`/:class:`~.grid.XarrayGrid` for a
given variable; a variable whose points are not described by coordinates of
their own, but numbered by a CF reduced Gaussian grid mapping, gets its grid
from :meth:`CoordinateGuesser.reduced_gaussian_grid`.
:class:`DefaultCoordinateGuesser` is the concrete
implementation used by :mod:`~.loader`; :class:`RulesCoordinateGuesser`
classifies the geography coordinates from user-supplied rules (the format
of an earthkit-data flavour ``rules`` section). Subclasses extend the
classification through the :meth:`CoordinateGuesser._guess_non_geo`,
:meth:`CoordinateGuesser._scalar`, :meth:`CoordinateGuesser._unsupported`
and :meth:`CoordinateGuesser._grid_cache_key` hooks.
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Hashable, Optional, Tuple

import xarray as xr

from earthkit.geo.utils.dotdict import DotDict

from .coordinates import (
    Coordinate,
    LatitudeCoordinate,
    LongitudeCoordinate,
    PointCoordinate,
    ScalarCoordinate,
    UnsupportedCoordinate,
    XCoordinate,
    YCoordinate,
    is_scalar,
)
from .grid import (
    Grid,
    MeshedGrid,
    MeshedXYGrid,
    MeshProjectionGrid,
    ReducedGaussianGrid,
    SelectiveReducedGaussianGrid,
    UnstructuredGrid,
    UnstructuredProjectionGrid,
    UnstructuredXYGrid,
    XarrayGrid,
)
from .reduced_gaussian import SELECTIVE, cf_reduced_gaussian_grid_name, cf_reduced_gaussian_subset

# CoordinateAttributes = namedtuple("CoordinateAttributes", ["axis", "name", "long_name", "standard_name", "units"])


LOG = logging.getLogger(__name__)

# the two CF attributes by which a reduced Gaussian grid mapping names the variable
# counting its points per latitude, and whether that variable holds the running
# total of them rather than the counts themselves (Appendix F, "Reduced Gaussian")
_POINTS_PER_LATITUDE = (("points_per_latitude", False), ("accumulated_points_per_latitude", True))


class CoordinateAttributes(DotDict):
    """Dict-with-attribute-access holding a coordinate's CF-relevant attributes.

    Populated with ``axis``, ``name``, ``long_name``, ``standard_name`` and
    ``units`` by :meth:`CoordinateGuesser._guess`.
    """

    pass


class CoordinateGuesser(ABC):
    """Class to guess the type of coordinates in a dataset."""

    def __init__(self, ds: xr.Dataset) -> None:
        """Initializes the CoordinateGuesser.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset to guess coordinates from.
        """
        self.ds = ds
        self._coordinate_cache: Dict[Hashable, Coordinate] = {}
        # the grids are keyed by their coordinate names, their dimensions and the
        # variable's grid mapping (see _mapping_key): two variables sharing one
        # coordinate pair only share a grid when they declare the same mapping
        self._grid_cache: Dict[Hashable, Grid] = {}
        # the grids of the variables whose points a reduced Gaussian grid mapping
        # numbers, keyed the same way; a mapping is read and its point index
        # checked once, whatever number of variables sit on it
        self._reduced_gaussian_cache: Dict[Hashable, Optional[XarrayGrid]] = {}

    def guess(self, c: xr.DataArray, coord: Hashable) -> Coordinate:
        """Guesses the type of a coordinate.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to guess.
        coord : Hashable
            The name of the coordinate.

        Returns
        -------
        Coordinate
            The guessed coordinate type.
        """
        if coord not in self._coordinate_cache:
            self._coordinate_cache[coord] = self._guess(c, coord)
        return self._coordinate_cache[coord]

    def _guess(self, coordinate: xr.DataArray, coord: Hashable) -> Coordinate:
        """Internal method to guess the type of a coordinate.

        Parameters
        ----------
        coordinate : xr.DataArray
            The coordinate to guess.
        coord : Hashable
            The name of the coordinate.

        Returns
        -------
        Coordinate
            The guessed coordinate type.
        """
        name = coordinate.name
        standard_name = getattr(coordinate, "standard_name", "").lower()
        axis = getattr(coordinate, "axis", "")
        long_name = getattr(coordinate, "long_name", "").lower()
        units = getattr(coordinate, "units", "")

        attributes = CoordinateAttributes(
            axis=axis,
            name=name,
            long_name=long_name,
            standard_name=standard_name,
            units=units,
        )

        d: Optional[Coordinate] = None

        d = self._is_point(coordinate, attributes)
        if d is not None:
            return d

        d = self._is_longitude(coordinate, attributes)
        if d is not None:
            return d

        d = self._is_latitude(coordinate, attributes)
        if d is not None:
            return d

        d = self._is_x(coordinate, attributes)
        if d is not None:
            return d

        d = self._is_y(coordinate, attributes)
        if d is not None:
            return d

        d = self._guess_non_geo(coordinate, attributes)
        if d is not None:
            return d

        if coordinate.shape in ((1,), tuple()):
            return self._scalar(coordinate)

        # LOG.warning(
        #     f"Coordinate {coord} not supported\n{axis=}, {name=},"
        #     f" {long_name=}, {standard_name=}, units\n\n{coordinate}\n\n{type(coordinate.values)} {coordinate.shape}"
        # )

        return self._unsupported(coordinate)

    def _guess_non_geo(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[Coordinate]:
        """Classify a coordinate that is not geography.

        Extension hook called by :meth:`_guess` after all geography
        matchers have declined the coordinate and before the scalar /
        unsupported fallback. The base implementation classifies nothing;
        subclasses may return a :class:`~.coordinates.Coordinate` (e.g. a
        time, step, level or ensemble coordinate of their own) or None to
        let the fallback apply.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to classify.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[Coordinate]
            The classified coordinate, or None.
        """
        return None

    def _scalar(self, c: xr.DataArray) -> ScalarCoordinate:
        """Build the coordinate for an unmatched scalar coordinate.

        Factory used by :meth:`_guess`; subclasses may override it to
        return their own :class:`~.coordinates.ScalarCoordinate` subclass.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to wrap.

        Returns
        -------
        ScalarCoordinate
            The scalar coordinate.
        """
        return ScalarCoordinate(c)

    def _unsupported(self, c: xr.DataArray) -> UnsupportedCoordinate:
        """Build the coordinate for an unmatched non-scalar coordinate.

        Factory used by :meth:`_guess`; subclasses may override it to
        return their own :class:`~.coordinates.UnsupportedCoordinate` subclass.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to wrap.

        Returns
        -------
        UnsupportedCoordinate
            The unsupported coordinate.
        """
        return UnsupportedCoordinate(c)

    def _grid_cache_key(
        self, first: Coordinate, second: Coordinate, dim_vars: Tuple[Hashable, ...], variable: Any
    ) -> tuple:
        """Build the key under which a coordinate pair's grid is cached.

        Factory used by :meth:`_lat_lon_provided` and :meth:`_x_y_provided`.
        The default combines the two coordinate names, the variable's grid
        dimensions and the variable's grid-mapping identity (see
        :meth:`_mapping_key`); a subclass that lets the grid depend on
        anything else from the variable must override this so the cache
        cannot mix up two variables whose grids differ.

        Parameters
        ----------
        first : Coordinate
            The first coordinate of the pair, in cache key order: ``lat``
            for a latitude/longitude pair, ``x`` for an x/y pair.
        second : Coordinate
            The second coordinate of the pair (``lon``, respectively ``y``).
        dim_vars : Tuple[Hashable, ...]
            The variable's grid dimensions.
        variable : Any
            The variable the grid is being derived for.

        Returns
        -------
        tuple
            A hashable cache key.
        """
        return (first.name, second.name, dim_vars, self._mapping_key(variable))

    def grid(self, coordinates: Any, variable: Any) -> Any:
        """Determines the grid type for the given coordinates and variable.

        Parameters
        ----------
        coordinates : Any
            The coordinates to determine the grid from.
        variable : Any
            The variable to determine the grid from.

        Returns
        -------
        Any
            The determined grid type.
        """
        lat = [c for c in coordinates if c.is_lat]
        lon = [c for c in coordinates if c.is_lon]

        latlon_grid = None
        if len(lat) == 1 and len(lon) == 1:
            latlon_grid = self._lat_lon_provided(lat, lon, variable)

        x = [c for c in coordinates if c.is_x]
        y = [c for c in coordinates if c.is_y]

        xy_grid = None
        if len(x) == 1 and len(y) == 1:
            xy_grid = self._x_y_provided(x, y, variable, strict=latlon_grid is None)

        if latlon_grid is not None or xy_grid is not None:
            return XarrayGrid(latlon_grid=latlon_grid, xy_grid=xy_grid)

        raise NotImplementedError(f"Cannot establish grid {coordinates}")

    def _check_dims(self, variable: Any, x_or_lon: Any, y_or_lat: Any) -> Tuple[Any, bool]:
        """Checks the dimensions of the variable against the coordinates.

        Parameters
        ----------
        variable : Any
            The variable to check.
        x_or_lon : Any
            The x or longitude coordinate.
        y_or_lat : Any
            The y or latitude coordinate.

        Returns
        -------
        Tuple[Any, bool]
            The checked dimensions and a flag indicating if the grid is unstructured.
        """
        x_dims = set(x_or_lon.variable.dims)
        y_dims = set(y_or_lat.variable.dims)
        variable_dims = set(variable.dims)

        if not (x_dims <= variable_dims) or not (y_dims <= variable_dims):
            raise ValueError(
                f"Dimensions do not match {variable.name}{variable.dims} !="
                f" {x_or_lon.name}{x_or_lon.variable.dims} and {y_or_lat.name}{y_or_lat.variable.dims}"
            )

        variable_dims = tuple(v for v in variable.dims if v in (x_dims | y_dims))
        if x_dims == y_dims:
            # It's unstructured
            return variable_dims, True

        if len(x_dims) == 1 and len(y_dims) == 1:
            # It's a mesh
            return variable_dims, False

        raise ValueError(
            f"Cannot establish grid for {variable.name}{variable.dims},"
            f" {x_or_lon.name}{x_or_lon.variable.dims},"
            f" {y_or_lat.name}{y_or_lat.variable.dims}"
        )

    def _grid_mapping(self, variable: Any) -> Optional[Dict[str, Any]]:
        """Resolve the CF grid mapping attributes of a variable.

        The variable's ``grid_mapping`` attribute names the (scalar)
        variable in the dataset carrying the mapping's attributes, e.g. a
        ``rotated_latitude_longitude`` or a ``lambert_azimuthal_equal_area``
        grid description.

        Parameters
        ----------
        variable : Any
            The variable to inspect.

        Returns
        -------
        Dict[str, Any] or None
            The grid mapping attributes, or None if the variable carries no
            usable ``grid_mapping`` reference.
        """
        name = variable.attrs.get("grid_mapping")
        if not isinstance(name, str) or not name:
            return None

        if name not in self.ds.variables:
            LOG.debug(f"Grid mapping variable '{name}' not found in dataset")
            return None

        return dict(self.ds[name].attrs)

    def _mapping_key(self, variable: Any) -> tuple:
        """Hashable identity of the grid mapping a variable declares, for the grid cache key.

        The grid derived for a coordinate pair can depend on the variable that uses
        it: an unstructured lat/lon grid carries the variable's resolved grid
        mapping (see :meth:`_grid_mapping`), and the projected or plain x/y grid
        follows the variable's own ``grid_mapping`` attribute, falling back to its
        ``crs`` attribute or to a candidate probed from the dataset. Everything the
        resolution takes from the variable is in those two attributes, and
        everything it takes from the dataset is the same for every variable of one
        call, so a key carrying them both cannot mix up two variables whose
        mappings differ, whatever order they are looked up in.

        Parameters
        ----------
        variable : Any
            The variable the grid is being derived for.

        Returns
        -------
        tuple
            A hashable standing for the variable's grid mapping declaration.
        """
        attrs = variable.attrs
        return (repr(attrs.get("grid_mapping")), repr(attrs.get("crs")))

    def _lat_lon_provided(self, lat: Any, lon: Any, variable: Any) -> Any:
        """Determines the grid type when latitude and longitude are provided.

        Parameters
        ----------
        lat : Any
            The latitude coordinate.
        lon : Any
            The longitude coordinate.
        variable : Any
            The variable to determine the grid from.

        Returns
        -------
        Any
            The determined grid type.
        """
        lat = lat[0]
        lon = lon[0]

        dim_vars, unstructured = self._check_dims(variable, lon, lat)

        # the grid of a coordinate pair keeps the variable's grid mapping, so the
        # cache cannot be keyed by the coordinates alone
        key = self._grid_cache_key(lat, lon, dim_vars, variable)

        if key in self._grid_cache:
            return self._grid_cache[key]

        grid: Grid = (
            UnstructuredGrid(lat, lon, dim_vars, grid_mapping=self._grid_mapping(variable))
            if unstructured
            else MeshedGrid(lat, lon, dim_vars)
        )

        self._grid_cache[key] = grid

        return grid

    def _x_y_provided(self, x: Any, y: Any, variable: Any, strict: bool = False) -> Any:
        """Determines the grid type when x and y coordinates are provided.

        Parameters
        ----------
        x : Any
            The x coordinate.
        y : Any
            The y coordinate.
        variable : Any
            The variable to determine the grid from.

        Returns
        -------
        Any
            The determined grid type.
        """
        x = x[0]
        y = y[0]

        dim_vars, unstructured = self._check_dims(variable, x, y)

        # the grid of a coordinate pair keeps the variable's grid mapping, so the
        # cache cannot be keyed by the coordinates alone
        key = self._grid_cache_key(x, y, dim_vars, variable)

        if key in self._grid_cache:
            return self._grid_cache[key]

        grid_mapping = variable.attrs.get("grid_mapping", None)

        if grid_mapping is None:
            LOG.debug(f"No 'grid_mapping' attribute provided for '{variable.name}'")
            LOG.debug("Trying to guess...")

            PROBE = {
                "prime_meridian_name",
                "reference_ellipsoid_name",
                "crs_wkt",
                "horizontal_datum_name",
                "semi_major_axis",
                "spatial_ref",
                "inverse_flattening",
                "semi_minor_axis",
                "geographic_crs_name",
                "GeoTransform",
                "grid_mapping_name",
                "longitude_of_prime_meridian",
            }
            candidate = None
            for v in self.ds.variables:
                var = self.ds[v]
                if not is_scalar(var):
                    continue

                if PROBE.intersection(var.attrs.keys()):
                    if candidate:
                        raise ValueError(f"Multiple candidates for 'grid_mapping': {candidate} and {v}")
                    candidate = v

            if candidate:
                LOG.debug(f"Using '{candidate}' as 'grid_mapping'")
                grid_mapping = candidate
            else:
                LOG.debug("Could not find a candidate for 'grid_mapping'")

        if grid_mapping is None:
            if "crs" in self.ds[variable.name].attrs:
                grid_mapping = self.ds[variable.name].attrs["crs"]
                LOG.debug(f"Using CRS {grid_mapping} from variable '{variable.name}' attributes")

        if grid_mapping is None:
            if "crs" in self.ds.attrs:
                grid_mapping = self.ds.attrs["crs"]
                LOG.debug(f"Using CRS {grid_mapping} from global attributes")

        grid: Optional[Grid] = None
        if grid_mapping is not None:
            if grid_mapping in self.ds.variables:
                grid_mapping = dict(self.ds[grid_mapping].attrs)
                if unstructured:
                    grid = UnstructuredProjectionGrid(x, y, grid_mapping, variable_dims=dim_vars)
                else:
                    grid = MeshProjectionGrid(x, y, grid_mapping, variable_dims=dim_vars)
            else:
                LOG.debug(f"Grid mapping variable '{grid_mapping}' not found in dataset")
        else:
            if unstructured:
                grid = UnstructuredXYGrid(x, y, dim_vars)
            else:
                grid = MeshedXYGrid(x, y, dim_vars)

        if grid is not None:
            self._grid_cache[key] = grid
            return grid

        LOG.debug("Could not find a candidate for 'grid_mapping'")

        if strict:
            raise NotImplementedError(f"Unstructured grid {x.name} {y.name}")
        else:
            return None

    def reduced_gaussian_grid(self, variable: Any) -> Optional[XarrayGrid]:
        """Determines the grid of a variable whose points a reduced Gaussian grid mapping numbers.

        A CF reduced Gaussian grid (Appendix F, "Reduced Gaussian") describes its
        points with a ``reduced_gaussian`` grid mapping and a
        ``reduced_gaussian_index`` coordinate instead of the latitude/longitude
        coordinates :meth:`grid` derives a grid from. The mapping and the table of
        points per latitude it names are read off the dataset, the point index
        among the coordinates of the variable itself; the two together locate the
        variable's points on the *global* grid the mapping describes. The whole
        index yields the :class:`~.grid.ReducedGaussianGrid` of that grid, the
        indices of a rectangular sub-area the same grid with an ``area`` in its
        spec — and with a :attr:`~.grid.ReducedGaussianGrid.point_order`
        permutation when eckit iterates the area in another order than the file —
        and any other valid subset the unstructured points of a
        :class:`~.grid.SelectiveReducedGaussianGrid`. Metadata an eckit grid
        cannot be built from — an unusable mapping or table, a point index which
        is no subset of the named grid — yields ``None``, like a dataset which
        does not describe such a grid at all.

        Parameters
        ----------
        variable : Any
            The variable to determine the grid from.

        Returns
        -------
        Optional[XarrayGrid]
            The grid the variable's points sit on, or ``None`` when the variable
            is not on a recognised reduced Gaussian grid.
        """
        mapping = self._grid_mapping(variable)
        if not isinstance(mapping, dict) or mapping.get("grid_mapping_name") != "reduced_gaussian":
            return None

        found = self._reduced_gaussian_index(variable)
        if found is None:
            return None
        name, dim = found

        # the mapping and the dimension the points are numbered by say all there is
        # to the grid, so two variables on one grid only check their index once:
        # a coordinate of one name holds one array of index values in one dataset,
        # and two variables indexing different point sets do it through
        # differently named coordinates, which the key tells apart
        key = (self._mapping_key(variable), name, dim)
        if key in self._reduced_gaussian_cache:
            return self._reduced_gaussian_cache[key]

        grid: Optional[XarrayGrid] = None
        pl = self._points_per_latitude(mapping)
        if pl is not None:
            counts, accumulated = pl

            # a few files spell the meridian their first point sits on as
            # longitude_of_first_point rather than the CF name
            # longitude_of_first_meridian: whichever a file uses, eckit builds an
            # O<n> grid from the Greenwich meridian, so it must be that one
            if "longitude_of_first_meridian" not in mapping and "longitude_of_first_point" in mapping:
                mapping = dict(mapping, longitude_of_first_meridian=mapping["longitude_of_first_point"])

            named = cf_reduced_gaussian_grid_name(mapping, counts, accumulated=accumulated)
            if named is None:
                LOG.debug("The mapping of '%s' describes no global reduced Gaussian grid", variable.name)
            else:
                subset = cf_reduced_gaussian_subset(named[0], self.ds[name].values)
                if subset is None:
                    LOG.debug("The points numbered by '%s' are no subset of the grid '%s'", name, named[0])
                elif subset.kind == SELECTIVE:
                    grid = SelectiveReducedGaussianGrid(subset.lats, subset.lons, (dim,))
                else:
                    # GLOBAL or RECTANGULAR: the subset's spec places its points,
                    # and its order says whether eckit iterates them the way the
                    # file holds them
                    grid = ReducedGaussianGrid(subset.spec, (dim,), point_order=subset.order)

        self._reduced_gaussian_cache[key] = grid
        return grid

    def _reduced_gaussian_index(self, variable: Any) -> Optional[Tuple[Hashable, Hashable]]:
        """The coordinate numbering the points of a reduced Gaussian variable, and the dimension it sits on.

        A CF reduced Gaussian grid numbers its points west to east and north to
        south with a ``reduced_gaussian_index`` (Appendix F, "Reduced Gaussian"):
        the 1D coordinate of the variable carrying that standard name, on a
        dimension the variable is indexed by. It is looked for among the
        coordinates the variable is associated with rather than among the
        dataset's variables: a variable of the dataset which happens to hold the
        same standard name says nothing about how the points of this one are
        numbered, and neither do two of them, with no way to tell which one the
        data follows.

        Parameters
        ----------
        variable : Any
            The variable to look the index coordinate up for.

        Returns
        -------
        Optional[Tuple[Hashable, Hashable]]
            The name of the index coordinate and the name of its dimension, or
            ``None`` when the variable is associated with no such coordinate, or
            with more than one.
        """
        candidates = []

        for name in variable.coords:
            if name not in self.ds.variables:
                continue

            c = self.ds[name]
            if c.ndim != 1:
                continue

            if str(c.attrs.get("standard_name", "")).strip().lower() != "reduced_gaussian_index":
                continue

            if c.dims[0] not in variable.dims:
                continue

            candidates.append((name, c.dims[0]))

        if len(candidates) != 1:
            if len(candidates) > 1:
                names = [str(name) for name, _ in candidates]
                LOG.debug(f"Variable '{variable.name}' is indexed by more than one reduced Gaussian index: {names}")
            return None

        return candidates[0]

    def _points_per_latitude(self, mapping: Dict[str, Any]) -> Optional[Tuple[Any, bool]]:
        """The longitude points per latitude a reduced Gaussian mapping names, and whether they are accumulated.

        Either of the two CF attributes names the variable the mapping counts the
        points of its lines with, along the dimension its ``latitude_dimension``
        names, when the mapping carries one.

        Parameters
        ----------
        mapping : Dict[str, Any]
            The attributes of the grid mapping, as resolved by :meth:`_grid_mapping`.

        Returns
        -------
        Optional[Tuple[Any, bool]]
            The values of the table as a 1D array and whether they are the
            accumulated points per latitude, or ``None`` when the mapping names no
            usable table.
        """
        for attr, accumulated in _POINTS_PER_LATITUDE:
            name = mapping.get(attr)
            if not isinstance(name, str) or not name:
                continue

            if name not in self.ds.variables:
                LOG.debug(f"Points per latitude variable '{name}' not found in dataset")
                return None

            pl = self.ds[name]
            if pl.ndim != 1:
                LOG.debug(f"Points per latitude variable '{name}' is not 1D")
                return None

            latitude = mapping.get("latitude_dimension")
            if isinstance(latitude, str) and latitude and latitude not in pl.dims:
                LOG.debug(f"Points per latitude variable '{name}' is not along '{latitude}'")
                return None

            return pl.values, accumulated

        return None

    @abstractmethod
    def _is_point(self, c: xr.DataArray, attributes: CoordinateAttributes) -> PointCoordinate | None:
        """Check if the coordinate identifies point/station data.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[PointCoordinate]
            The PointCoordinate if matched, else None.
        """
        pass

    @abstractmethod
    def _is_longitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LongitudeCoordinate]:
        """Checks if the coordinate is a longitude.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LongitudeCoordinate]
            The LongitudeCoordinate if matched, else None.
        """
        pass

    @abstractmethod
    def _is_latitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LatitudeCoordinate]:
        """Checks if the coordinate is a latitude.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LatitudeCoordinate]
            The LatitudeCoordinate if matched, else None.
        """
        pass

    @abstractmethod
    def _is_x(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[XCoordinate]:
        """Checks if the coordinate is an x coordinate.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[XCoordinate]
            The XCoordinate if matched, else None.
        """
        pass

    @abstractmethod
    def _is_y(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[YCoordinate]:
        """Checks if the coordinate is a y coordinate.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[YCoordinate]
            The YCoordinate if matched, else None.
        """
        pass


class DefaultCoordinateGuesser(CoordinateGuesser):
    """Default implementation of CoordinateGuesser."""

    def __init__(self, ds: xr.Dataset) -> None:
        """Initializes the DefaultCoordinateGuesser.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset to guess coordinates from.
        """
        super().__init__(ds)

    def _is_point(self, c: xr.DataArray, attributes: CoordinateAttributes) -> PointCoordinate | None:
        """Check if the coordinate identifies point/station data.

        Matches on ``standard_name`` or ``name`` being one of
        ``"location"``, ``"cell"``, ``"id"``, ``"station"``, ``"poi"`` or
        ``"point"``.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[PointCoordinate]
            The PointCoordinate if matched, else None.
        """
        if attributes.standard_name in ["location", "cell", "id", "station", "poi", "point"]:
            return PointCoordinate(c)

        if attributes.name in ["location", "cell", "id", "station", "poi", "point"]:  # WeatherBench
            return PointCoordinate(c)

        return None

    def _is_longitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LongitudeCoordinate]:
        """Checks if the coordinate is a longitude.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LongitudeCoordinate]
            The LongitudeCoordinate if matched, else None.
        """
        # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#longitude-coordinate

        if attributes.standard_name == "longitude":
            return LongitudeCoordinate(c)

        if attributes.long_name == "longitude" and attributes.units == "degrees_east":
            return LongitudeCoordinate(c)

        if attributes.name == "longitude":  # WeatherBench
            return LongitudeCoordinate(c)

        if attributes.name in ("lon", "grid_longitude"):
            return LongitudeCoordinate(c)

        return None

    def _is_latitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LatitudeCoordinate]:
        """Checks if the coordinate is a latitude.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LatitudeCoordinate]
            The LatitudeCoordinate if matched, else None.
        """
        # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html#latitude-coordinate

        if attributes.standard_name == "latitude":
            return LatitudeCoordinate(c)

        if attributes.long_name == "latitude" and attributes.units == "degrees_north":
            return LatitudeCoordinate(c)

        if attributes.name == "latitude":  # WeatherBench
            return LatitudeCoordinate(c)

        if attributes.name in ("lat", "grid_latitude"):
            return LatitudeCoordinate(c)

        return None

    def _is_x(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[XCoordinate]:
        """Checks if the coordinate is an x coordinate.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[XCoordinate]
            The XCoordinate if matched, else None.
        """
        if attributes.standard_name in ["projection_x_coordinate", "grid_longitude"]:
            return XCoordinate(c)

        if attributes.name == "x":
            return XCoordinate(c)

        return None

    def _is_y(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[YCoordinate]:
        """Checks if the coordinate is a y coordinate.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[YCoordinate]
            The YCoordinate if matched, else None.
        """
        if attributes.standard_name in ["projection_y_coordinate", "grid_latitude"]:
            return YCoordinate(c)

        if attributes.name == "y":
            return YCoordinate(c)

        return None


class RulesCoordinateGuesser(CoordinateGuesser):
    """Implementation of CoordinateGuesser that matches geography coordinates against rules.

    The ``rules`` mapping has the same format as the ``rules`` section of an
    earthkit-data flavour (a mapping from a coordinate key to a single rule
    or a list of rules, each rule a mapping from coordinate attributes --
    ``axis``, ``name``, ``long_name``, ``standard_name``, ``units`` -- to the
    string values they must equal for the rule to match). Only the geography
    keys ``point``, ``longitude``, ``latitude``, ``x`` and ``y`` are consulted
    here; other keys (e.g. ``time`` or ``level``) are ignored and should be
    handled by a subclass through :meth:`~.CoordinateGuesser._guess_non_geo`
    or its own matchers.

    Unlike :class:`DefaultCoordinateGuesser`, the geography matchers use
    *only* the rules: a coordinate that matches no rule is classified by the
    :meth:`~.CoordinateGuesser._guess_non_geo` hook or falls back to
    scalar/unsupported, without the CF-convention heuristics.
    """

    def __init__(self, ds: xr.Dataset, rules: Optional[Dict[str, Any]] = None) -> None:
        """Initializes the RulesCoordinateGuesser.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset to guess coordinates from.
        rules : Dict[str, Any], optional
            The rules to use for guessing, keyed by coordinate role
            (``"point"``, ``"longitude"``, ``"latitude"``, ``"x"``, ``"y"``
            and, for subclasses, further non-geography keys).
        """
        super().__init__(ds)
        self.rules: Dict[str, Any] = dict(rules) if rules else {}

    def _match(self, c: xr.DataArray, key: str, attributes: CoordinateAttributes) -> Optional[Dict[str, Any]]:
        """Matches the coordinate against the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to match.
        key : str
            The key to match in the rules.
        attributes : CoordinateAttributes
            The values to match against.

        Returns
        -------
        Optional[Dict[str, Any]]
            The matched rule if any, else None.
        """
        if key not in self.rules:
            return None

        rules = self.rules[key]

        if not isinstance(rules, list):
            rules = [rules]

        for rule in rules:
            ok = True
            for k, v in rule.items():
                if isinstance(v, str) and attributes.get(k) != v:
                    ok = False
            if ok:
                return rule

        return None

    def _is_point(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[PointCoordinate]:
        """Checks if the coordinate is a point coordinate using the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[PointCoordinate]
            The PointCoordinate if matched, else None.
        """
        if self._match(c, "point", attributes):
            return PointCoordinate(c)

        return None

    def _is_longitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LongitudeCoordinate]:
        """Checks if the coordinate is a longitude using the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LongitudeCoordinate]
            The LongitudeCoordinate if matched, else None.
        """
        if self._match(c, "longitude", attributes):
            return LongitudeCoordinate(c)

        return None

    def _is_latitude(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[LatitudeCoordinate]:
        """Checks if the coordinate is a latitude using the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[LatitudeCoordinate]
            The LatitudeCoordinate if matched, else None.
        """
        if self._match(c, "latitude", attributes):
            return LatitudeCoordinate(c)

        return None

    def _is_x(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[XCoordinate]:
        """Checks if the coordinate is an x coordinate using the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[XCoordinate]
            The XCoordinate if matched, else None.
        """
        if self._match(c, "x", attributes):
            return XCoordinate(c)

        return None

    def _is_y(self, c: xr.DataArray, attributes: CoordinateAttributes) -> Optional[YCoordinate]:
        """Checks if the coordinate is a y coordinate using the rules.

        Parameters
        ----------
        c : xr.DataArray
            The coordinate to check.
        attributes : CoordinateAttributes
            The attributes of the coordinate.

        Returns
        -------
        Optional[YCoordinate]
            The YCoordinate if matched, else None.
        """
        if self._match(c, "y", attributes):
            return YCoordinate(c)

        return None
