# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import logging

from earthkit.geo.grids._regrid.gridspec import normalise_grid_spec

LOG = logging.getLogger(__name__)


# TODO: This is a temporary wrapper for the Grid interface
class GridWrapper:
    """Thin wrapper around an ``eckit.geo.Grid`` used to build output geography.

    Normalises access to a grid built from a grid spec (or an existing
    ``Grid`` instance) and adds helpers to extract flat or distinct
    lat/lon arrays for a given field shape.
    """

    def __init__(self, grid_spec):
        """Initialise the wrapper.

        Parameters
        ----------
        grid_spec : Any
            A grid spec (dict/str) or an existing ``eckit.geo.Grid`` instance.
        """
        from eckit.geo import Grid

        if isinstance(grid_spec, Grid):
            self._grid = grid_spec
        else:
            self._grid = Grid(grid_spec)
        self._grid_spec = grid_spec

    def __getattr__(self, name):
        """Delegate unknown attribute access to the wrapped ``Grid``."""
        return getattr(self._grid, name)

    def to_latlons(self):
        """Get the flat latitude and longitude arrays for the grid.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            The latitude and longitude arrays.
        """
        import numpy as np

        lat, lon = self._grid.to_latlons()
        return np.array(lat), np.array(lon)

    @property
    def grid_spec(self):
        """Any: The original grid spec passed to the wrapper."""
        # TODO: for grid specs like {'grid': 'O32', 'area': [87.863799, 0.0, -87.863799, 357.5]}
        # The Grid.spec is not correct so we cannot return self.spec
        return self._grid_spec

    def is_spectral(self):
        """bool: Whether the grid is spectral (always False here)."""
        return False

    def to_distinct_latlons(self, field_shape):
        """Get the distinct (1D) latitude and longitude arrays for a 2D field.

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
        if len(self._grid.shape) == 2:
            lat, lon = self.to_latlons()
            lat = lat.reshape(self._grid.shape)
            lon = lon.reshape(self._grid.shape)
            d_lat = self._distinct_lats(lat)
            if d_lat is not None:
                d_lon = self._distinct_lons(lon)
                if d_lon is not None and len(d_lat) == field_shape[0] and len(d_lon) == field_shape[1]:
                    return d_lat, d_lon

        return None, None

    @staticmethod
    def _distinct_lats(lats):
        """Get the distinct (1D) latitude array for a 2D meshed latitude array.

        Parameters
        ----------
        lats : np.ndarray
            2D array of latitudes.

        Returns
        -------
        Optional[np.ndarray]
            The distinct latitudes per row, or None if the rows are not
            regularly spaced.
        """
        import numpy as np

        assert len(lats.shape) == 2
        rows = lats.shape[0]
        r = np.ones(rows)
        if rows > 0:
            for i in range(rows):
                vals = lats[i, :]
                delta = np.diff(vals)
                if np.allclose(delta, delta[0]):
                    r[i] = vals[0]
                else:
                    return None
            return r
        return None

    @staticmethod
    def _distinct_lons(lons):
        """Get the distinct (1D) longitude array for a 2D meshed longitude array.

        Parameters
        ----------
        lons : np.ndarray
            2D array of longitudes.

        Returns
        -------
        Optional[np.ndarray]
            The distinct longitudes per column, or None if the columns are
            not regularly spaced.
        """
        import numpy as np

        assert len(lons.shape) == 2
        cols = lons.shape[1]
        r = np.ones(cols)
        if cols > 0:
            for i in range(cols):
                vals = lons[:, i]
                delta = np.diff(vals)
                if np.allclose(delta, delta[0]):
                    r[i] = vals[0]
                else:
                    return None
            return r
        return None


class XarrayGeographyBuilder:
    """Builds output geography (dims/coords) for a regridded xarray variable.

    Wraps an output grid spec and derives the dimension names, coordinate
    arrays and coordinate-to-dimension mapping to attach to the regridded
    result.
    """

    def __init__(self, grid_spec):
        """Initialise the builder.

        Parameters
        ----------
        grid_spec : Any
            The output grid spec (dict/str) or an ``eckit.geo.Grid`` instance.
        """
        grid_spec = normalise_grid_spec(grid_spec)
        self.grid = GridWrapper(grid_spec)
        self.grid_spec = grid_spec

    @property
    def shape(self):
        """Tuple[int, ...]: The shape of the output grid."""
        return self.grid.shape

    def geo_dims(self):
        """Determine the geographical dimensions of the dataset."""
        num = len(self.shape)
        if num >= 2:
            return ["latitude", "longitude"]
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
                try:
                    lat, lon = self.grid.to_distinct_latlons(field_shape)
                    if (
                        lat is not None
                        and lon is not None
                        and len(lat) == field_shape[0]
                        and len(lon) == field_shape[1]
                    ):
                        coords["latitude"] = lat
                        coords["longitude"] = lon
                        coords_dim["latitude"] = ("latitude",)
                        coords_dim["longitude"] = ("longitude",)
                        dims["latitude"] = lat.size
                        dims["longitude"] = lon.size
                        assert coords["latitude"].size == field_shape[0]
                        assert coords["longitude"].size == field_shape[1]
                        assert dims["latitude"] == field_shape[0]
                        assert dims["longitude"] == field_shape[1]
                except Exception as e:
                    print(e)
                    pass

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
