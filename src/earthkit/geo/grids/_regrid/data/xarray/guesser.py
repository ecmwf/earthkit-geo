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
given variable. :class:`DefaultCoordinateGuesser` is the concrete
implementation used by :mod:`~.loader`.
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
    UnstructuredGrid,
    UnstructuredProjectionGrid,
    UnstructuredXYGrid,
    XarrayGrid,
)

# CoordinateAttributes = namedtuple("CoordinateAttributes", ["axis", "name", "long_name", "standard_name", "units"])


LOG = logging.getLogger(__name__)


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
        self._grid_cache: Dict[Hashable, Grid] = {}

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

        if coordinate.shape in ((1,), tuple()):
            return ScalarCoordinate(coordinate)

        # LOG.warning(
        #     f"Coordinate {coord} not supported\n{axis=}, {name=},"
        #     f" {long_name=}, {standard_name=}, units\n\n{coordinate}\n\n{type(coordinate.values)} {coordinate.shape}"
        # )

        return UnsupportedCoordinate(coordinate)

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

        if (lat.name, lon.name, dim_vars) in self._grid_cache:
            return self._grid_cache[(lat.name, lon.name, dim_vars)]

        grid: Grid = UnstructuredGrid(lat, lon, dim_vars) if unstructured else MeshedGrid(lat, lon, dim_vars)

        self._grid_cache[(lat.name, lon.name, dim_vars)] = grid

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

        if (x.name, y.name, dim_vars) in self._grid_cache:
            return self._grid_cache[(x.name, y.name, dim_vars)]

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
                    grid = UnstructuredProjectionGrid(x, y, grid_mapping)
                else:
                    grid = MeshProjectionGrid(x, y, grid_mapping)
            else:
                LOG.debug(f"Grid mapping variable '{grid_mapping}' not found in dataset")
        else:
            if unstructured:
                grid = UnstructuredXYGrid(x, y, dim_vars)
            else:
                grid = MeshedXYGrid(x, y, dim_vars)

        if grid is not None:
            self._grid_cache[(x.name, y.name, dim_vars)] = grid
            return grid

        LOG.error("Could not find a candidate for 'grid_mapping'")

        if strict:
            raise NotImplementedError(f"Unstructured grid {x.name} {y.name}")
        else:
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
