# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging

from ..utils import create_grid_object

LOG = logging.getLogger(__name__)


# This is a temporary wrapper for the Grid interface
class GridWrapper:
    """Thin wrapper around an ``eckit.geo.Grid`` used to build output geography.

    Normalises access to a grid built from a grid spec (or an existing
    ``Grid`` instance) and adds helpers to extract flat or distinct
    lat/lon arrays for a given field shape.
    """

    def __init__(self, grid):
        """Initialise the wrapper.

        Parameters
        ----------
        grid_spec : Any
            A grid spec (dict/str) or an existing ``eckit.geo.Grid`` instance.
        """
        self._grid = create_grid_object(grid)
        self._latlons = None

    def __getattr__(self, name):
        """Delegate unknown attribute access to the wrapped ``Grid``."""
        return getattr(self._grid, name)

    def to_latlons(self):
        """Get the flat latitude and longitude arrays for the grid.

        The arrays are computed by the grid once and then cached: the
        geography builder asks for them while deciding the output dimensions
        and again while building the coordinates, and computing them can be
        expensive on a large grid.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            The latitude and longitude arrays.
        """
        if self._latlons is None:
            import numpy as np

            lat, lon = self._grid.to_latlons()
            self._latlons = (np.array(lat), np.array(lon))
        return self._latlons

    def is_spectral(self):
        """bool: Whether the grid is spectral (always False here)."""
        return False

    def to_distinct_latlons(self, field_shape):
        """Get the distinct (1D) latitude and longitude arrays for a 2D field.

        The grid is a regular mesh when its latitudes are constant along
        every row of the field and its longitudes constant along every
        column (which a rotated grid is not, however few points its rows
        and columns have); the distinct arrays are then the latitude of
        each row and the longitude of each column.

        Parameters
        ----------
        field_shape : Tuple[int, int]
            The shape of the field the grid is used for.

        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            The distinct latitude and longitude arrays, or ``(None, None)``
            if the grid is not a regular mesh matching ``field_shape``.
        """
        import numpy as np

        if len(self._grid.shape) == 2:
            lat, lon = self.to_latlons()
            lat = lat.reshape(self._grid.shape)
            lon = lon.reshape(self._grid.shape)
            if (
                np.allclose(lat, lat[:, :1])
                and np.allclose(lon, lon[:1, :])
                and lat.shape[0] == field_shape[0]
                and lon.shape[1] == field_shape[1]
            ):
                return lat[:, 0], lon[0, :]

        return None, None


class XarrayOutputGeographyBuilder:
    """Builds output geography (dims/coords) for a regridded xarray variable.

    Wraps an output grid spec and derives the dimension names, coordinate
    arrays and coordinate-to-dimension mapping to attach to the regridded
    result.
    """

    def __init__(self, grid):
        """Initialise the builder with the given output grid.

        Parameters
        ----------
        grid : Any
            The output grid spec (dict/str) or an ``eckit.geo.Grid`` instance.
        """
        self.grid = GridWrapper(grid)
        self._distinct = None

    @property
    def shape(self):
        """Tuple[int, ...]: The shape of the output grid."""
        return self.grid.shape

    def _distinct_latlons(self):
        """Decide whether the grid carries distinct 1D latitudes and longitudes.

        This is the decision :meth:`coords` makes between latitude/longitude
        coordinates on their own dimensions and 2D ones on ``(y, x)``, cached
        so that :meth:`geo_dims` and :meth:`coords` share it, as well as the
        latitudes and longitudes the grid computes for it.

        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            The distinct latitude and longitude arrays when the grid is a
            regular mesh matching its shape, or ``(None, None)``.
        """
        if self._distinct is None:
            lat = lon = None
            if len(self.shape) == 2:
                try:
                    d_lat, d_lon = self.grid.to_distinct_latlons(self.shape)
                except Exception as e:
                    LOG.debug("Could not extract distinct lat/lon coordinates: %s", e)
                else:
                    if (
                        d_lat is not None
                        and d_lon is not None
                        and len(d_lat) == self.shape[0]
                        and len(d_lon) == self.shape[1]
                    ):
                        lat, lon = d_lat, d_lon
            self._distinct = (lat, lon)
        return self._distinct

    def geo_dims(self):
        """Determine the geographical dimensions of the dataset.

        The decision is the one :meth:`coords` makes, shared through
        :meth:`_distinct_latlons`: a 2D grid has distinct 1D latitudes and
        longitudes (constant along every row and column, which a rotated
        grid is not) only when its coordinates are built on the ``latitude``
        and ``longitude`` dimensions; any other 2D grid carries 2D
        coordinates on ``(y, x)``.

        Returns
        -------
        List[str]
            ``["latitude", "longitude"]`` for a 2D grid with distinct 1D
            latitudes and longitudes, ``["y", "x"]`` for any other 2D grid
            (e.g. a rotated one), and ``["values"]`` for a 1D grid.

        Raises
        ------
        ValueError
            If the grid has no geography this builder supports.
        """
        num = len(self.shape)
        if num >= 2:
            lat, lon = self._distinct_latlons()
            if lat is not None and lon is not None:
                return ["latitude", "longitude"]
            return ["y", "x"]
        if num == 1:
            return ["values"]

        raise ValueError("Geography is not supported.")

    def coords(self):
        """Build the output coordinate arrays for the grid.

        Returns
        -------
        Tuple[Dict[str, int], Dict[str, np.ndarray], Dict[str, Tuple[str, ...]]]
            The output dimension sizes, the coordinate arrays (e.g.
            ``latitude``/``longitude``), and the dimensions each coordinate
            is defined on.
        """
        import math

        field_shape = self.grid.shape

        coords = {}
        dims = {}
        coords_dim = {}

        if self.grid.is_spectral():
            if len(field_shape) == 1:
                dims["values"] = field_shape[0]
        else:
            if len(field_shape) == 1:
                dims["values"] = field_shape[0]
                try:
                    lat, lon = self.grid.to_latlons()
                    if lat is not None and lon is not None:
                        coords["latitude"] = lat
                        coords["longitude"] = lon
                        coords_dim = {k: ("values",) for k in coords}
                except Exception:
                    pass
            elif len(field_shape) == 2:
                lat, lon = self._distinct_latlons()
                if lat is not None and lon is not None:
                    coords["latitude"] = lat
                    coords["longitude"] = lon
                    coords_dim["latitude"] = ("latitude",)
                    coords_dim["longitude"] = ("longitude",)
                    dims["latitude"] = lat.size
                    dims["longitude"] = lon.size

                if not coords or not dims:
                    lat, lon = self.grid.to_latlons()
                    # print("to_latlons:", type(lat), type(lon))
                    if lat is not None and lon is not None:
                        lat = lat.reshape(field_shape)
                        lon = lon.reshape(field_shape)
                        coords["latitude"] = lat
                        coords["longitude"] = lon
                        coords_dim = {k: ("y", "x") for k in coords}
                        dims["y"] = field_shape[0]
                        dims["x"] = field_shape[1]
                        # print("field_shape:", field_shape, lat.shape, lon.shape)
                        assert coords["latitude"].shape == field_shape
                        assert coords["longitude"].shape == field_shape
            else:
                raise ValueError(f"Unsupported field shape {field_shape}")

        for k, v in coords.items():
            assert k in coords_dim, f"{k=}, {coords_dim=}"
            assert all(x in dims for x in coords_dim[k]), f"{k=}, {coords_dim=} {dims=}"
            assert v.size == math.prod([dims[x] for x in coords_dim[k]])

        return dims, coords, coords_dim
