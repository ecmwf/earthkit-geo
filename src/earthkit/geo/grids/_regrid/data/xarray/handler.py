# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""xarray data handler for the regrid dispatch mechanism.

Provides :class:`XarrayDataHandler`, the entry point used by the generic
regrid machinery (see ``..handler.DataHandler``) when the data to regrid is
an ``xarray.Dataset`` or ``xarray.DataArray``. It uses :mod:`~.loader` to
identify the geographical variables and their input grid, builds the output
geography via :class:`XarrayGeographyBuilder`, and regrids each variable
with ``xarray.apply_ufunc`` (dask-aware) delegating the actual point-to-point
interpolation to ``..numpy.NumpyDataHandler``.
"""

import functools
import logging

from earthkit.geo.grids._regrid.gridspec import normalise_grid_spec
from earthkit.geo.utils import ensure_list

from ..handler import DataHandler

LOG = logging.getLogger(__name__)


# TODO: This is a temporary wrapper to use the grid interface
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


class XarrayDataHandler(DataHandler):
    """Data handler that regrids ``xarray.Dataset``/``xarray.DataArray`` values."""

    @staticmethod
    def match(values):
        """Check whether ``values`` is an xarray object this handler can process.

        Parameters
        ----------
        values : Any
            The data to check.

        Returns
        -------
        bool
            True if ``xarray`` is loaded and ``values`` is an
            ``xr.DataArray`` or ``xr.Dataset``.
        """
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("xarray"):
            return False

        try:
            import xarray as xr

            return isinstance(values, (xr.DataArray, xr.Dataset))
        except Exception:
            return False

    @staticmethod
    def get_out_geo(grid):
        """Get the output geography from the out_grid."""
        out_grid = grid
        if out_grid is None:
            raise ValueError("grid must be provided")

        out_geo = XarrayGeographyBuilder(out_grid)
        return out_geo

    @staticmethod
    def add_geo_coords(ds, out_geo):
        """Add the geographical coordinates to the dataset."""
        dims, coords, coords_dim = out_geo.coords()

        import xarray as xr

        for k, v in coords.items():
            c_dims = {x: dims[x] for x in coords_dim[k]}
            ds.coords[k] = xr.Variable(c_dims, v)

        return ds

    @staticmethod
    def update_attributes(ds, out_geo):
        """Update the earthkit grid_spec attribute of a regridded dataset/array, if present."""
        # TODO: this is a temporary workaround to only set the grid_spec attribute
        # for datasets/arrays created from earthkit-data. The problem is that
        # the earthkit specific attributes cannot be written to NetCDF files. So
        # we avoid adding/updating it if it is not already present in the dataset.
        # This whole approach needs to be rethought once we have a better
        # way to handle the grid spec in xarray.

        # This low level check should be replaced by a more robust way to
        # determine if the dataset/array is created from earthkit-data.
        has_earthkit = False
        if "_earthkit" in ds.attrs:
            has_earthkit = True
        else:
            import xarray as xr

            if isinstance(ds, xr.Dataset):
                for var in ds.data_vars.values():
                    if "_earthkit" in var.attrs:
                        has_earthkit = True
                        break

        if has_earthkit:
            if hasattr(ds, "earthkit"):
                try:
                    ds = ds.earthkit.set({"geography.grid_spec": out_geo.grid_spec})
                except Exception:
                    # TODO: temporary workaround for when storing the new grid spec
                    # fails (e.g. the accessor cannot serialise it). Leaving the
                    # original _earthkit attribute in place would be worse than
                    # having none at all, since it still describes the input grid
                    # and would make the regridded result look like it was on the
                    # source grid. As a safeguard we drop the attribute instead.
                    # This should be revisited once the grid spec is handled
                    # properly in xarray.
                    if isinstance(ds, xr.Dataset):
                        for var in ds.data_vars.values():
                            if "_earthkit" in var.attrs:
                                del var.attrs["_earthkit"]

                    elif isinstance(ds, xr.DataArray):
                        if "_earthkit" in ds.attrs:
                            del ds.attrs["_earthkit"]

        return ds

    def regrid(self, values, in_grid=None, out_grid=None, **kwargs):
        """Regrid an xarray Dataset or DataArray onto a new grid.

        Parameters
        ----------
        values : xr.Dataset or xr.DataArray
            The data to regrid.
        in_grid : Any, optional
            The input grid spec, used when it cannot be determined from
            ``values`` (e.g. from earthkit-data metadata or CF coordinates).
        out_grid : Any, optional
            The output grid spec.
        **kwargs
            Extra keyword arguments, including optional ``in_dims`` and
            ``out_dims`` (input/output geographical dimension names) and any
            arguments forwarded to the underlying point-to-point regridding
            method.

        Returns
        -------
        xr.Dataset or xr.DataArray
            The regridded data, with updated geographical coordinates and
            (where possible) an updated grid_spec attribute. The return type
            matches the type of ``values``.
        """
        ds = values

        kwargs = kwargs.copy()
        in_grid_arg = in_grid
        out_grid_arg = out_grid
        in_dims_arg = kwargs.pop("in_dims", None)
        out_dims_arg = kwargs.pop("out_dims", None)

        import xarray as xr

        from .loader import variables as get_variables

        input_is_dataset = isinstance(ds, xr.Dataset)
        if not input_is_dataset:
            ds = ds.to_dataset()

        variables = get_variables(ds, user_ek_grid=in_grid_arg)

        for v in variables:
            if v.ek_grid is None:
                if in_grid_arg is None:
                    raise ValueError(
                        f"Could not determine grid for variable {v.name} from dataset. Please provide an "
                        f"'in_grid' argument."
                    )
                else:
                    raise ValueError(
                        f"Could not determine grid for variable {v.name} from dataset or from 'in_grid' argument."
                    )
            if not v.geo_dims:
                if in_dims_arg is None:
                    raise ValueError(
                        f"Could not determine grid dimensions for variable {v.name} from dataset. Please provide "
                        "an 'in_dims' argument."
                    )
                else:
                    raise ValueError(
                        f"Could not determine grid dimensions for variable {v.name} from dataset and "
                        "'in_dims' argument is provided but invalid."
                    )

        # the output geography builder which can provide the output grid and the output coordinates
        out_geo = self.get_out_geo(out_grid_arg)
        out_grid = out_geo.grid._grid

        out_dims = out_dims_arg
        if out_dims is None:
            out_dims = out_geo.geo_dims()

        if out_dims is None:
            raise ValueError(f"Could not determine geography related output dimensions: {values.dims}")

        out_dims = ensure_list(out_dims)

        ds_out = xr.Dataset()

        res_out_grid = None
        for v in variables:
            ds_out[v.name], res_out_grid_v = self._regrid_variable(v, out_geo, in_dims_arg, out_dims, **kwargs)
            if res_out_grid is None:
                res_out_grid = res_out_grid_v

        # for a DataArray input will return a single DataArray instead of a Dataset
        if not input_is_dataset:
            ds_out = ds_out[list(ds.keys())[0]]

        # The output geography might have changed, so we need to create a new geography builder
        # with the new grid spec
        out_geo = XarrayGeographyBuilder(res_out_grid)

        ds_out = self.add_geo_coords(ds_out, out_geo)
        ds_out = self.update_attributes(ds_out, out_geo)

        return ds_out

    def _regrid_variable(self, variable, out_geo, in_dims_arg, out_dims, **kwargs):
        """Regrid a single variable using ``xarray.apply_ufunc``.

        Parameters
        ----------
        variable : Variable
            The variable to regrid, as produced by :func:`~.loader.variables`.
        out_geo : XarrayGeographyBuilder
            The output geography.
        in_dims_arg : Optional[List[str]]
            Fallback input grid dimension names, used when they cannot be
            determined from ``variable``.
        out_dims : List[str]
            The output grid dimension names.
        **kwargs
            Extra keyword arguments forwarded to
            ``NumpyDataHandler.regrid``.

        Returns
        -------
        Tuple[xr.DataArray, Any]
            The regridded ``xr.DataArray`` and the (possibly updated) output
            grid spec as returned by the point-to-point regrid method.
        """
        from ..numpy import NumpyDataHandler

        # regrid() can change the specified output gridspec.
        # This is a workaround to get the returned output gridscpec from regrid().
        class _RegridMethod:
            def __init__(self, in_grid, out_grid, **kwargs):
                self.out_grid = out_grid
                self.method = functools.partial(
                    NumpyDataHandler().regrid,
                    in_grid=in_grid,
                    out_grid=out_grid,
                    **kwargs,
                )

            def __call__(self, vals):
                # TODO: ensure it is thread safe
                vals, self.out_grid = self.method(vals)
                return vals

        var = variable.variable
        in_dims = variable.geo_dims
        if in_dims is None:
            in_dims = in_dims_arg

        in_dims = ensure_list(in_dims)

        exclude_dims = set()
        if set(in_dims) == set(out_dims):
            exclude_dims = set(in_dims)

        in_grid = variable.ek_grid

        method = _RegridMethod(in_grid.spec, out_geo.grid_spec, **kwargs)

        import xarray as xr

        res = xr.apply_ufunc(
            method,
            var,
            input_core_dims=[in_dims],
            output_core_dims=[out_dims],
            exclude_dims=exclude_dims,
            vectorize=True,
            dask="parallelized",
            dask_gufunc_kwargs={
                "output_sizes": {dim: out_geo.shape[i] for i, dim in enumerate(out_dims)},
                "allow_rechunk": True,
            },
            output_dtypes=[var.dtype],
            keep_attrs="identical",
        )

        res_out_grid = method.out_grid.copy() if method.out_grid is not None else method.out_grid
        return res, res_out_grid


handler = XarrayDataHandler
