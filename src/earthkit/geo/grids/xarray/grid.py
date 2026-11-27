# (C) Copyright 2024 Anemoi contributors.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


"""Grid classes deriving lat/lon and x/y point arrays from xarray coordinates.

Given the :class:`~.coordinates.Coordinate` objects identified by
:mod:`~.guesser`, these classes compute the flat ``(lat, lon)`` and/or
``(x, y)`` point arrays describing a variable's geographical grid, covering
meshed (structured) and unstructured layouts, plain x/y grids, and
projected (e.g. CRS-based) x/y grids. :class:`XarrayGrid` is the top-level
wrapper combining an optional lat/lon grid with an optional x/y grid, and
:class:`ReducedGaussianGrid` the grid of a CF reduced Gaussian variable whose
points are described by a grid mapping rather than by coordinates.
"""

import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Tuple

import numpy as np
from earthkit.utils.decorators import thread_safe_cached_property

from .grid_spec import (
    healpix_grid_spec_from_points,
    proj_projection,
    reduced_gg_grid_spec_from_points,
    regular_ll_grid_spec,
    regular_ll_grid_spec_from_points,
    regular_xy_grid_spec,
    rotation_projection,
)

LOG = logging.getLogger(__name__)

# how far the points of the eckit grid built from a candidate gridspec may sit
# from the dataset's own coordinates before the spec is rejected, in degrees
_GRID_MATCH_ATOL = 1e-3

# the eckit regular_ll orders which iterate the grid with its longitude axis
# outermost, keyed by the order which iterates it with the latitude axis
# outermost in the same directions ("i+j+" for an ascending latitude axis,
# the absent default for a descending one)
_ROTATED_ORDER_SWAP = {"i+j+": "j+i+", None: "j-i+"}


def _rotated_role(coord: Any) -> Optional[str]:
    """The rotated frame role of a coordinate, ``"lat"`` / ``"lon"``, or None if it has none.

    The axes of a CF rotated grid carry the ``grid_latitude`` /
    ``grid_longitude`` standard names; failing that, an X/Y axis in angular
    units is one.
    """
    standard_name = str(coord.attrs.get("standard_name", "")).lower()
    axis = str(coord.attrs.get("axis", "")).upper()
    units = str(coord.attrs.get("units", "")).lower()

    if standard_name == "grid_latitude" or (axis == "Y" and "degree" in units):
        return "lat"

    if standard_name == "grid_longitude" or (axis == "X" and "degree" in units):
        return "lon"

    return None


def _projected_role(coord: Any) -> Optional[str]:
    """The projection plane role of a coordinate, ``"x"`` / ``"y"``, or None if it has none.

    The axes of a CF projected grid carry the ``projection_x_coordinate`` /
    ``projection_y_coordinate`` standard names; failing that, an X/Y axis in
    units that are not angular is one. A coordinate with no units at all is
    not: the roles are guessed from CF hints only.
    """
    standard_name = str(coord.attrs.get("standard_name", "")).lower()
    axis = str(coord.attrs.get("axis", "")).upper()
    units = str(coord.attrs.get("units", "")).lower()

    if standard_name == "projection_y_coordinate" or (axis == "Y" and units and "degree" not in units):
        return "y"

    if standard_name == "projection_x_coordinate" or (axis == "X" and units and "degree" not in units):
        return "x"

    return None


def _frame_axes(lat: Any, role_of: Any) -> Optional[Dict[str, Any]]:
    """Get the two 1D frame axes a 2D lat/lon coordinate is defined on, by role.

    A CF grid whose coordinates are given as a 2D field of true latitudes and
    longitudes is defined on a frame of its own, be it a rotated one or the
    plane of a projection: the frame is carried by the two 1D dimension
    coordinates of that field.

    Parameters
    ----------
    lat : Any
        The 2D latitude of the grid, as its xarray DataArray.
    role_of : Any
        Callable mapping a dimension coordinate to its role in the frame, or
        to ``None`` when it does not describe one.

    Returns
    -------
    Dict[str, Any] or None
        The frame axes as ``{role: coord}``, or ``None`` if the field is not
        2D over exactly two recognised, 1D dimension coordinates.
    """
    if lat.ndim != 2:
        return None

    axes: Dict[str, Any] = {}
    for dim in lat.dims:
        coord = lat[dim] if dim in lat.coords else None
        if coord is None or coord.ndim != 1:
            return None

        name = role_of(coord)
        if name is None or name in axes:
            return None
        axes[name] = coord

    if len(axes) != 2:
        return None

    return axes


def _grid_spec_matches(spec: Any, shape: Any, lats: Any, lons: Any) -> bool:
    """Whether the eckit-geo grid of ``spec`` reproduces the points of a grid.

    Parameters
    ----------
    spec : Any
        The candidate gridspec.
    shape : Any
        The expected ``(y, x)``, i.e. ``(latitude, longitude)``, shape of the
        grid built from it.
    lats : Any
        The flat latitudes to match, in the grid's own point order.
    lons : Any
        The flat longitudes to match, in the grid's own point order.

    Returns
    -------
    bool
        True if the spec builds a grid of ``shape`` whose points sit within
        :data:`_GRID_MATCH_ATOL` degrees of ``lats``/``lons`` (longitudes
        compared modulo 360), False if it does not or cannot be built. A
        non-finite coordinate on either side makes the check fail.

    Notes
    -----
    Every coordinate takes part in the comparison: a spec is only returned
    when each point of the eckit grid is within the tolerance of the data's,
    which a NaN cannot satisfy. A dataset whose coordinates carry a NaN, even
    at a single point, therefore keeps falling back to the unstructured
    treatment, as it did before any gridspec was inferred from it; the check
    cannot tell a coordinate that is missing from one that disagrees with the
    spec, and a masked-out point must not vouch for the whole grid.
    """
    try:
        from eckit.geo import Grid

        grid = Grid(spec)
        ek_lats, ek_lons = grid.to_latlons()
    except Exception:
        return False

    if tuple(grid.shape) != tuple(shape):
        return False

    ek_lats = np.asarray(ek_lats, dtype=np.float64).reshape(-1)
    ek_lons = np.asarray(ek_lons, dtype=np.float64).reshape(-1)
    lats = np.asarray(lats, dtype=np.float64).reshape(-1)
    lons = np.asarray(lons, dtype=np.float64).reshape(-1)

    if ek_lats.size != lats.size or ek_lons.size != lons.size:
        return False

    if not (np.isfinite(ek_lats).all() and np.isfinite(ek_lons).all()):
        return False

    if not (np.isfinite(lats).all() and np.isfinite(lons).all()):
        return False

    if not bool(np.all(np.abs(ek_lats - lats) <= _GRID_MATCH_ATOL)):
        return False

    dlon = np.abs((ek_lons - lons + 180.0) % 360.0 - 180.0)
    return bool(np.all(dlon <= _GRID_MATCH_ATOL))


class XarrayGrid:
    """Top-level grid combining an optional lat/lon grid and/or x/y grid.

    Wraps the :class:`LatLonGrid` and/or :class:`XYGrid` derived for a
    variable, preferring the lat/lon grid where both are available.
    """

    def __init__(self, latlon_grid=None, xy_grid=None) -> None:
        """Initialise the wrapper.

        Parameters
        ----------
        latlon_grid : LatLonGrid, optional
            The variable's lat/lon grid, if available.
        xy_grid : XYGrid, optional
            The variable's x/y grid, if available.
        """
        self.latlon_grid = latlon_grid
        self.xy_grid = xy_grid

    @property
    def latlons(self):
        """Tuple[Any, Any]: The flat ``(lat, lon)`` point arrays.

        From :attr:`latlon_grid` if available, otherwise from
        :attr:`xy_grid`.
        """
        if self.latlon_grid is not None:
            return self.latlon_grid.latlons
        else:
            return self.xy_grid.latlons

    @property
    def xys(self):
        """Tuple[Any, Any]: The flat ``(x, y)`` point arrays, from :attr:`xy_grid`.

        Raises
        ------
        NotImplementedError
            If no ``xy_grid`` is available.
        """
        if self.xy_grid is not None:
            return self.xy_grid.xys
        else:
            raise NotImplementedError("XarrayGrid: xy grid not available")

    # Properly implement it for all grid types
    @thread_safe_cached_property
    def bbox(self):
        """Tuple[float, float, float, float]: The ``(north, south, east, west)`` bounding box."""
        lat, lon = self.latlons
        lat = lat.flatten()
        lon = lon.flatten()

        north = np.amax(lat)
        west = np.amin(lon)
        south = np.amin(lat)
        east = np.amax(lon)

        return (north, south, east, west)

    @property
    def grid_type(self) -> str:
        """Str or None: The eckit-geo grid type, from :attr:`latlon_grid` when available."""
        if self.latlon_grid is not None:
            return self.latlon_grid.grid_type
        elif self.xy_grid is not None:
            return self.xy_grid.grid_type
        else:
            return None

    @property
    def variable_dims(self):
        """Tuple[str, ...] or None: The variable's grid dimension names.

        From :attr:`latlon_grid` if available, otherwise from
        :attr:`xy_grid`, or None if neither is available.
        """
        if self.latlon_grid is not None:
            return self.latlon_grid.variable_dims
        elif self.xy_grid is not None:
            return self.xy_grid.variable_dims
        else:
            return None

    @property
    def grid_spec(self):
        """Get the grid specification from :attr:`latlon_grid` or :attr:`xy_grid`."""
        if self.latlon_grid is not None:
            return self.latlon_grid.grid_spec
        elif self.xy_grid is not None:
            return self.xy_grid.grid_spec
        else:
            return None


class Grid(ABC):
    """Abstract base class for grid structures."""

    def __init__(self) -> None:
        pass

    @property
    @abstractmethod
    def latlons(self) -> Tuple[Any, Any]:
        """Get the grid points."""
        pass

    @property
    @abstractmethod
    def xys(self) -> Tuple[Any, Any]:
        """Get the grid points."""
        pass

    @property
    @abstractmethod
    def grid_type(self) -> str:
        """str: Get the eckit-geo grid type."""
        pass

    @property
    def grid_spec(self):
        """Get the grid specification."""
        return None


class LatLonGrid(Grid):
    """Grid class for latitude and longitude coordinates."""

    def __init__(self, lat: Any, lon: Any, variable_dims: Any) -> None:
        """Initialize the LatLonGrid class.

        Parameters
        ----------
        lat : Any
            The latitudes.
        lon : Any
            The longitudes.
        variable_dims : Any
            The variable dimensions.
        """
        super().__init__()
        self.lat = lat
        self.lon = lon
        self.variable_dims = variable_dims

    @property
    def xys(self) -> Tuple[Any, Any]:
        """Get the x and y points for the lat-lon grid."""
        raise NotImplementedError("LatLonGrid: xy grid not available")


class XYGrid(Grid):
    """Grid class for x and y coordinates."""

    def __init__(self, x: Any, y: Any) -> None:
        """Initialize the XYGrid class.

        Parameters
        ----------
        x : Any
            The x-coordinates.
        y : Any
            The y-coordinates.
        """
        self.x = x
        self.y = y

    @property
    def grid_type(self):
        """None: X/Y grids have no known eckit-geo grid type."""
        return None


class MeshedGrid(LatLonGrid):
    """Grid class for meshed latitude and longitude coordinates."""

    @thread_safe_cached_property
    def latlons(self) -> Tuple[Any, Any]:
        """Get the grid points for the meshed grid.

        Raises
        ------
        NotImplementedError
            If ``variable_dims`` is not ``(lon_name, lat_name)`` or
            ``(lat_name, lon_name)``.
        """
        if self.variable_dims == (self.lon.variable.name, self.lat.variable.name):
            lat, lon = np.meshgrid(
                self.lat.variable.values,
                self.lon.variable.values,
            )
        elif self.variable_dims == (self.lat.variable.name, self.lon.variable.name):
            lon, lat = np.meshgrid(
                self.lon.variable.values,
                self.lat.variable.values,
            )

        else:
            raise NotImplementedError(f"MeshedGrid.grid_points: unrecognised variable_dims {self.variable_dims}")

        return lat.flatten(), lon.flatten()

    @property
    def grid_type(self):
        """str: Always ``"regular-ll"``."""
        return "regular-ll"

    @thread_safe_cached_property
    def grid_spec(self):
        """Dict or None: The gridspec for a regular lat/lon grid, or ``None`` if not regular."""
        return regular_ll_grid_spec(self.lat.variable.values, self.lon.variable.values)


class UnstructuredGrid(LatLonGrid):
    """Grid class for unstructured latitude and longitude coordinates."""

    def __init__(self, lat: Any, lon: Any, variable_dims: Any, grid_mapping: Any = None) -> None:
        """Initialize the UnstructuredGrid class.

        Parameters
        ----------
        lat : Any
            The latitudes.
        lon : Any
            The longitudes.
        variable_dims : Any
            The variable dimensions.
        grid_mapping : Any, optional
            The CF grid mapping attributes of the variable, if it carries a
            ``grid_mapping`` reference (see :class:`~.guesser.CoordinateGuesser`).
        """
        super().__init__(lat, lon, variable_dims)
        assert len(lat) == len(lon), (len(lat), len(lon))
        self.variable_dims = variable_dims
        self.grid_dims = lat.variable.dims
        assert lon.variable.dims == self.grid_dims, (lon.variable.dims, self.grid_dims)
        assert set(self.variable_dims) == set(self.grid_dims), (self.variable_dims, self.grid_dims)
        self.grid_mapping = grid_mapping

    @thread_safe_cached_property
    def latlons(self) -> Tuple[Any, Any]:
        """Get the grid points for the unstructured grid."""
        assert 1 <= len(self.variable_dims) <= 2

        if len(self.variable_dims) == 1:
            return self.lat.variable.values.flatten(), self.lon.variable.values.flatten()

        if len(self.variable_dims) == 2 and self.variable_dims == self.grid_dims:
            return self.lat.variable.values.flatten(), self.lon.variable.values.flatten()

        LOG.warning(
            "UnstructuredGrid: variable indexing %s does not match grid indexing %s",
            self.variable_dims,
            self.grid_dims,
        )

        lat = self.lat.variable.values.transpose().flatten()
        lon = self.lon.variable.values.transpose().flatten()

        return lat, lon

    @property
    def grid_type(self):
        """str: Always ``"unstructured"``."""
        return "unstructured"

    @thread_safe_cached_property
    def grid_spec(self):
        """Dict or None: The gridspec if the points form a known grid, else ``None``.

        Tries, in order, a regular lat/lon grid, a HEALPix grid, a reduced
        Gaussian grid, a CF rotated latitude/longitude grid and a CF projected
        (regular x/y) grid.
        """
        lat, lon = self.latlons
        return (
            regular_ll_grid_spec_from_points(lat, lon)
            or healpix_grid_spec_from_points(lat, lon)
            or reduced_gg_grid_spec_from_points(lat, lon)
            or self._rotated_ll_grid_spec()
            or self._projected_xy_grid_spec()
        )

    def _rotated_frame_axes(self) -> Optional[Dict[str, Any]]:
        """Get the 1D rotated frame axes of the grid, by role.

        The 2D ``lat``/``lon`` coordinates of a CF rotated grid are defined
        on two 1D dimension coordinates carrying the rotated frame (see
        :func:`_rotated_role`).

        Returns
        -------
        Dict[str, Any] or None
            The frame axes as ``{"lat": rlat, "lon": rlon}`` xarray
            DataArrays, or ``None`` if the two axes are not found.
        """
        return _frame_axes(self.lat.variable, _rotated_role)

    def _projected_frame_axes(self) -> Optional[Dict[str, Any]]:
        """Get the 1D axes of the projection plane the grid sits in, by role.

        The 2D ``lat``/``lon`` coordinates of a CF projected grid are defined
        on two 1D dimension coordinates carrying the plane's x and y axes, in
        the units of that plane, typically metres (see :func:`_projected_role`).

        Returns
        -------
        Dict[str, Any] or None
            The frame axes as ``{"x": x, "y": y}`` xarray DataArrays, or
            ``None`` if the two axes are not found.
        """
        return _frame_axes(self.lat.variable, _projected_role)

    def _rotated_ll_grid_spec(self) -> Optional[dict]:
        """Build the gridspec of a CF rotated latitude/longitude grid, or ``None``.

        The variable's grid mapping must describe a rotated latitude/longitude
        grid (see :func:`~.grid_spec.rotation_projection`), its two 1D frame
        axes must be regular (see :func:`~.grid_spec.regular_ll_grid_spec`)
        and the eckit-geo grid built from the resulting spec must reproduce
        the grid's own :attr:`latlons` within :data:`_GRID_MATCH_ATOL`
        degrees. A spec is only returned once that check passes, so a grid
        whose metadata disagrees with its coordinates keeps falling back to
        the unstructured treatment.
        """
        projection = rotation_projection(self.grid_mapping)
        if projection is None:
            return None

        axes = self._rotated_frame_axes()
        if axes is None:
            return None

        rlat, rlon = axes["lat"], axes["lon"]
        # eckit's rotation projects the spec's latitude slot as the rotated
        # latitude, so the frame axes always sit in their own slots
        spec = regular_ll_grid_spec(rlat.values, rlon.values)
        if spec is None:
            return None

        # eckit iterates the points with the latitude axis outermost unless
        # the order says otherwise; swap the order so that the spec iterates
        # in the order the data indexes its grid (see _ROTATED_ORDER_SWAP)
        if self.variable_dims[0] == rlon.dims[0]:
            order = _ROTATED_ORDER_SWAP.get(spec.get("order"))
            if order is None:
                return None
            spec = dict(spec, order=order)
        elif self.variable_dims[0] != rlat.dims[0]:
            return None

        spec = dict(spec, projection=projection)

        lat, lon = self.latlons
        # the eckit grid is shaped (rotated latitudes, rotated longitudes),
        # whichever of them the data varies fastest
        if not _grid_spec_matches(spec, (len(rlat), len(rlon)), lat, lon):
            LOG.debug("UnstructuredGrid: rotated gridspec %s does not match the grid points", spec)
            return None

        return spec

    def _projected_xy_grid_spec(self) -> Optional[dict]:
        """Build the gridspec of a CF projected (regular x/y) grid, or ``None``.

        The variable's grid mapping must carry a PROJ description of the
        projection (see :func:`~.grid_spec.proj_projection`), its two 1D frame
        axes must be regular (see :func:`~.grid_spec.regular_xy_grid_spec`) and
        the eckit-geo grid built from the resulting spec must reproduce the
        grid's own :attr:`latlons` within :data:`_GRID_MATCH_ATOL` degrees. A
        spec is only returned once that check passes, so a grid whose projection
        metadata disagrees with its coordinates, or whose axes the data varies
        fastest along in the other order (eckit scans a ``regular_xy`` grid with
        its x axis fastest), keeps falling back to the unstructured treatment.
        """
        projection = proj_projection(self.grid_mapping)
        if projection is None:
            return None

        axes = self._projected_frame_axes()
        if axes is None:
            return None

        x, y = axes["x"], axes["y"]
        spec = regular_xy_grid_spec(x.values, y.values, projection)
        if spec is None:
            return None

        lat, lon = self.latlons
        # the eckit grid is shaped (y values, x values), with x fastest
        if not _grid_spec_matches(spec, (len(y), len(x)), lat, lon):
            LOG.debug("UnstructuredGrid: projected gridspec %s does not match the grid points", spec)
            return None

        return spec


class RawXYGrid(XYGrid):
    """Grid class for raw x and y coordinates."""

    def __init__(self, x: Any, y: Any, variable_dims: Any) -> None:
        """Initialize the RawXYGrid class.

        Parameters
        ----------
        x : Any
            The x-coordinates.
        y : Any
            The y-coordinates.
        variable_dims : Any
            The variable dimensions.
        """
        super().__init__(x, y)
        self.variable_dims = variable_dims

    @property
    def latlons(self) -> Tuple[Any, Any]:
        """Get the lat and lon points for the raw xy grid."""
        raise NotImplementedError("RawXYGrid: latlon not available")


class MeshedXYGrid(RawXYGrid):
    """Grid class for meshed x and y coordinates."""

    @thread_safe_cached_property
    def xys(self) -> Tuple[Any, Any]:
        """Get the grid points for the meshed grid."""
        if self.variable_dims == (self.x.variable.name, self.y.variable.name):
            y, x = np.meshgrid(
                self.y.variable.values,
                self.x.variable.values,
            )
        elif self.variable_dims == (self.y.variable.name, self.x.variable.name):
            x, y = np.meshgrid(
                self.x.variable.values,
                self.y.variable.values,
            )

        else:
            raise NotImplementedError(f"MeshedGrid.grid_points: unrecognised variable_dims {self.variable_dims}")

        return x.flatten(), y.flatten()


class UnstructuredXYGrid(RawXYGrid):
    """Grid class for unstructured x and y coordinates."""

    def __init__(self, x: Any, y: Any, variable_dims: Any) -> None:
        """Initialize the UnstructuredGrid class.

        Parameters
        ----------
        x : Any
            The x-coordinates.
        y : Any
            The y-coordinates.
        variable_dims : Any
            The variable dimensions.
        """
        super().__init__(x, y, variable_dims)
        assert len(x) == len(y), (len(x), len(y))
        self.variable_dims = variable_dims
        self.grid_dims = x.variable.dims
        assert y.variable.dims == self.grid_dims, (y.variable.dims, self.grid_dims)
        assert set(self.variable_dims) == set(self.grid_dims), (self.variable_dims, self.grid_dims)

    @thread_safe_cached_property
    def xys(self) -> Tuple[Any, Any]:
        """Get the grid points for the unstructured grid."""
        assert 1 <= len(self.variable_dims) <= 2

        if len(self.variable_dims) == 1:
            return self.x.variable.values.flatten(), self.y.variable.values.flatten()

        if len(self.variable_dims) == 2 and self.variable_dims == self.grid_dims:
            return self.x.variable.values.flatten(), self.y.variable.values.flatten()

        LOG.warning(
            "UnstructuredGrid: variable indexing %s does not match grid indexing %s",
            self.variable_dims,
            self.grid_dims,
        )

        x = self.x.variable.values.transpose().flatten()
        y = self.y.variable.values.transpose().flatten()

        return x, y


class ProjectionGrid(XYGrid):
    """Grid class for projected x and y coordinates."""

    def __init__(self, x: Any, y: Any, projection: Any, variable_dims: Any = None) -> None:
        """Initialize the ProjectionGrid class.

        Parameters
        ----------
        x : Any
            The x-coordinates.
        y : Any
            The y-coordinates.
        projection : Any
            The projection information.
        variable_dims : Any, optional
            The variable dimensions, in the order the variable indexes them;
            needed to lay the grid points out in the order the data holds them.
        """
        super().__init__(x, y)
        self.projection = projection
        self.variable_dims = variable_dims

    def transformer(self) -> Any:
        """Get the transformer for the projection.

        Returns
        -------
        Any
            The transformer.
        """
        from pyproj import CRS, Transformer

        if isinstance(self.projection, dict):
            data_crs = CRS.from_cf(self.projection)
        else:
            data_crs = self.projection

        wgs84_crs = CRS.from_epsg(4326)  # WGS84

        return Transformer.from_crs(data_crs, wgs84_crs, always_xy=True)


class MeshProjectionGrid(ProjectionGrid):
    """Grid class for meshed projected coordinates."""

    @thread_safe_cached_property
    def _scans_x_fastest(self) -> Optional[bool]:
        """Bool or None: whether the data holds the mesh with its x axis fastest.

        The mesh is the product of the two 1D axes: with the variable indexed
        ``(y, x)`` its x axis varies fastest, and with it indexed ``(x, y)`` its
        y axis does. The order is read off the dimensions the two axes are
        defined on, not their names, so a variable indexed by dimensions named
        differently from the coordinates carrying them (as after
        :meth:`xarray.Dataset.rename_dims`) keeps its own layout, as in
        :class:`UnstructuredGrid`. ``None``: neither order matches the
        variable's indexing, which cannot be told at all.
        """
        x_dims = tuple(self.x.variable.dims)
        y_dims = tuple(self.y.variable.dims)

        if self.variable_dims == y_dims + x_dims:
            return True
        if self.variable_dims == x_dims + y_dims:
            return False
        return None

    @thread_safe_cached_property
    def _points(self) -> Tuple[Any, Any]:
        """The flat ``(x, y)`` values of the mesh, in the order the data holds them.

        The order is the one :attr:`_scans_x_fastest` reads off the axes'
        dimensions. An indexing that cannot be told keeps the ``(y, x)`` order,
        the one a projected grid is usually stored in, and no
        :attr:`grid_spec` is claimed for it.
        """
        xv = self.x.variable.values
        yv = self.y.variable.values

        if self._scans_x_fastest is False:
            y, x = np.meshgrid(yv, xv)
        else:
            if self._scans_x_fastest is None:
                LOG.warning(
                    "MeshProjectionGrid: variable indexing %s matches neither the axes' dimensions %s and %s,"
                    " taking the points with the x axis fastest",
                    self.variable_dims,
                    self.y.variable.dims,
                    self.x.variable.dims,
                )
            x, y = np.meshgrid(xv, yv)

        return x.flatten(), y.flatten()

    @thread_safe_cached_property
    def latlons(self) -> Tuple[Any, Any]:
        """Get the grid points for the mesh projection grid."""
        transformer = self.transformer()
        xv, yv = self._points
        lon, lat = transformer.transform(xv, yv)
        return lat.flatten(), lon.flatten()

    @thread_safe_cached_property
    def xys(self) -> Tuple[Any, Any]:
        """Get the x and y points for the mesh projection grid."""
        return self._points

    @property
    def grid_type(self):
        """str: ``"regular_xy"`` when a gridspec is detected, else ``"unstructured"``."""
        return "regular_xy" if self.grid_spec is not None else "unstructured"

    @thread_safe_cached_property
    def grid_spec(self) -> Optional[dict]:
        """Dict or None: The gridspec if the axes form a regular projected grid, else ``None``.

        The grid mapping must carry a PROJ description of the projection (see
        :func:`~.grid_spec.proj_projection`) and the axes must be regular (see
        :func:`~.grid_spec.regular_xy_grid_spec`), each range in the direction
        its axis runs. The eckit-geo grid built from the resulting spec must also
        reproduce the grid's own :attr:`latlons`, in their point order, within
        :data:`_GRID_MATCH_ATOL` degrees. A spec is only returned once that check
        passes, so a grid whose PROJ description disagrees with the CF attributes
        it is made of, or whose data indexing does not hold the mesh with its x
        axis fastest, known or not (eckit scans a ``regular_xy`` grid with its x
        axis fastest and cannot be told otherwise), falls back to the
        unstructured points of :attr:`latlons`, as a grid with no projection at
        all did before.
        """
        projection = proj_projection(self.projection)
        if projection is None:
            return None

        if self._scans_x_fastest is not True:
            LOG.debug("MeshProjectionGrid: indexing %s does not scan the mesh x fastest", self.variable_dims)
            return None

        xv = self.x.variable.values
        yv = self.y.variable.values

        spec = regular_xy_grid_spec(xv, yv, projection)
        if spec is None:
            return None

        lat, lon = self.latlons
        # the eckit grid is shaped (y values, x values), with x fastest
        if not _grid_spec_matches(spec, (len(yv), len(xv)), lat, lon):
            LOG.debug("MeshProjectionGrid: projected gridspec %s does not match the grid points", spec)
            return None

        return spec


class UnstructuredProjectionGrid(ProjectionGrid):
    """Grid class for unstructured projected coordinates."""

    def __init__(self, x: Any, y: Any, projection: Any, variable_dims: Any = None) -> None:
        """Initialize the UnstructuredProjectionGrid class.

        Parameters
        ----------
        x : Any
            The x-coordinates.
        y : Any
            The y-coordinates.
        projection : Any
            The projection information.
        variable_dims : Any, optional
            The variable dimensions, in the order the variable indexes them;
            needed to lay the grid points out in the order the data holds them.
            Defaults to the dimensions the two coordinates are defined on.
        """
        super().__init__(x, y, projection, variable_dims=variable_dims)
        assert len(x) == len(y), (len(x), len(y))
        self.grid_dims = tuple(x.variable.dims)
        assert y.variable.dims == self.grid_dims, (y.variable.dims, self.grid_dims)
        if self.variable_dims is None:
            self.variable_dims = self.grid_dims
        # the point order is read off the dimensions, which are compared as the
        # tuple xarray spells them whatever way the grid is built
        self.variable_dims = tuple(self.variable_dims)
        assert set(self.variable_dims) == set(self.grid_dims), (self.variable_dims, self.grid_dims)

    @thread_safe_cached_property
    def _points(self) -> Tuple[Any, Any]:
        """Tuple[Any, Any]: The flat ``(x, y)`` values, one point per index of the grid dimension.

        The two coordinates are paired, index by index, and not crossed into a
        mesh of their values: a grid held as the two 1D axes of a plane is a
        :class:`MeshProjectionGrid`. The pairs come out in the order the data
        indexes them, which for a variable indexed in the other order than the
        coordinates are defined on means transposing them, as in
        :class:`UnstructuredXYGrid`.
        """
        assert 1 <= len(self.variable_dims) <= 2, self.variable_dims

        if len(self.variable_dims) == 1 or self.variable_dims == self.grid_dims:
            return self.x.variable.values.flatten(), self.y.variable.values.flatten()

        LOG.warning(
            "UnstructuredProjectionGrid: variable indexing %s does not match grid indexing %s",
            self.variable_dims,
            self.grid_dims,
        )

        x = self.x.variable.values.transpose().flatten()
        y = self.y.variable.values.transpose().flatten()

        return x, y

    @thread_safe_cached_property
    def latlons(self) -> Tuple[Any, Any]:
        """Tuple[Any, Any]: The flat ``(lat, lon)`` points, the :attr:`_points` in degrees.

        The plane coordinates are transformed with the CRS the projection
        describes (see :meth:`ProjectionGrid.transformer`).
        """
        transformer = self.transformer()
        xv, yv = self._points
        lon, lat = transformer.transform(xv, yv)
        return lat.flatten(), lon.flatten()

    @thread_safe_cached_property
    def xys(self) -> Tuple[Any, Any]:
        """Tuple[Any, Any]: The flat ``(x, y)`` points, one pair per index of the grid dimension."""
        if self.projection == "epsg:4326":
            # WGS84, no transformation needed; the pair has always been returned
            # the other way round for that projection, and latlons does not
            # follow it, taking its points from _points
            x, y = self._points
            return y, x

        return self._points

    @property
    def grid_type(self) -> str:
        """str: Always ``"unstructured"``: the points are a scatter over the plane."""
        return "unstructured"

    @property
    def grid_spec(self) -> None:
        """None: A scatter of projected points has no regular gridspec.

        The loader therefore builds the eckit-geo grid out of the
        :attr:`latlons`, as it does for every grid without a spec.
        """
        return None


class ReducedGaussianGrid(XarrayGrid):
    """Grid class for a CF reduced Gaussian grid described by a grid mapping.

    A reduced Gaussian grid written with the CF conventions (Appendix F,
    "Reduced Gaussian") holds no per-point latitudes and longitudes: its points
    are numbered by a ``reduced_gaussian_index`` dimension, and the grid they sit
    on is described by the mapping's table of points per latitude. The gridspec
    of a global grid so described, octahedral ``O<n>`` or classical ``N<n>``, is
    all eckit-geo needs to place the points, so :attr:`latlons` come from the
    grid built out of it rather than from the dataset.
    """

    def __init__(self, spec: Dict[str, Any], grid_dims: Any) -> None:
        """Initialize the ReducedGaussianGrid class.

        Parameters
        ----------
        spec : Dict[str, Any]
            The gridspec of the grid, e.g. ``{"grid": "O1280"}`` or
            ``{"grid": "N80"}`` (see
            :func:`~.grid_spec.cf_reduced_gaussian_grid_spec`).
        grid_dims : Any
            The names of the dimension the points are indexed by, in the order
            the variable indexes it.
        """
        super().__init__(latlon_grid=None, xy_grid=None)
        self.grid_dims = tuple(grid_dims)
        self._spec = dict(spec)

    @property
    def variable_dims(self) -> Tuple[str, ...]:
        """Tuple[str, ...]: The dimension the variable's points are indexed by."""
        return self.grid_dims

    @property
    def grid_spec(self) -> Dict[str, Any]:
        """Dict[str, Any]: The gridspec the grid mapping describes."""
        return self._spec

    @property
    def grid_type(self) -> str:
        """str: Always ``"reduced_gg"``: the eckit-geo grid type of a reduced Gaussian grid."""
        return "reduced_gg"

    @thread_safe_cached_property
    def latlons(self) -> Tuple[Any, Any]:
        """Tuple[Any, Any]: The flat ``(lat, lon)`` points, from the eckit-geo grid of :attr:`grid_spec`.

        The dataset does not hold the points of a reduced Gaussian grid it only
        numbers, so they are read off the grid built from the gridspec, in the
        order that grid iterates them: west to east and north to south, which is
        the order of the ``reduced_gaussian_index`` the data is indexed by.
        """
        from eckit.geo import Grid as EckitGrid

        lat, lon = EckitGrid(self._spec).to_latlons()
        return np.asarray(lat, dtype=np.float64).flatten(), np.asarray(lon, dtype=np.float64).flatten()
