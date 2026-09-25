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

from earthkit.geo.utils import ensure_list

from ..handler import DataHandler
from ..utils import create_grid_object
from .output import XarrayOutputGeographyBuilder

LOG = logging.getLogger(__name__)


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
    def add_geo_coords(ds, out_coords):
        """Add the geographical coordinates to the dataset."""
        dims, coords, coords_dim = out_coords

        import xarray as xr

        for k, v in coords.items():
            c_dims = {x: dims[x] for x in coords_dim[k]}
            ds.coords[k] = xr.Variable(c_dims, v)

        return ds

    @staticmethod
    def update_attributes(ds, out_grid):
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
                    ds = ds.earthkit.set({"geography.grid_spec": out_grid.spec})
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
            The input grid spec. When not provided, it is determined from
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

        Raises
        ------
        ValueError
            If the input grid or input dimensions cannot be determined and are not provided.
        """
        ds = values
        kwargs = kwargs.copy()

        backend = self.backend_from_kwargs(kwargs)
        in_grid = create_grid_object(in_grid)
        out_grid = create_grid_object(out_grid)

        if out_grid is None:
            raise ValueError("out_grid must be specified")

        in_dims_arg = kwargs.pop("in_dims", None)
        out_dims_arg = kwargs.pop("out_dims", None)

        import xarray as xr

        from .loader import variables as get_variables

        input_is_dataset = isinstance(ds, xr.Dataset)
        if not input_is_dataset:
            ds = ds.to_dataset()

        variables = get_variables(ds, user_ek_grid=in_grid)

        print("variables:", variables)

        for v in variables:
            if v.ek_grid is None:
                if in_grid is None:
                    raise ValueError(
                        f"Could not determine grid for variable {v.name} from dataset. Please provide an "
                        f"'in_grid' argument."
                    )
                else:
                    raise ValueError("Could not determine grid from 'in_grid' argument.")
            if not v.geo_dims:
                if in_dims_arg is None:
                    raise ValueError(
                        f"Could not determine grid dimensions for variable {v.name} from dataset. Please provide "
                        "an 'in_dims' argument."
                    )
                else:
                    raise ValueError(
                        f"Dimensions for variable {v.name} from dataset do not match the provided 'in_dims' argument."
                    )

        # the output geography builder which can provide the output grid and the output coordinates
        out_geo = XarrayOutputGeographyBuilder(out_grid)

        out_dims = out_dims_arg
        if out_dims is None:
            out_dims = out_geo.geo_dims()

        if out_dims is None:
            raise ValueError(f"Could not determine geography related output dimensions: {values.dims}")

        out_dims = ensure_list(out_dims)

        ds_out = xr.Dataset()

        res_out_grid = None
        for v in variables:
            ds_out[v.name], res_out_grid_v = self._regrid_variable(
                backend, v, out_grid, in_dims_arg, out_dims, **kwargs
            )
            if res_out_grid is None:
                res_out_grid = res_out_grid_v

        # for a DataArray input will return a single DataArray instead of a Dataset
        if not input_is_dataset:
            ds_out = ds_out[list(ds.keys())[0]]

        print("res_out_grid:", res_out_grid)

        # The output geography might have changed, so we need to create a new geography builder
        # with the new grid spec
        out_grid = create_grid_object(res_out_grid)
        out_geo = XarrayOutputGeographyBuilder(out_grid)
        ds_out = self.add_geo_coords(ds_out, out_geo.coords())
        ds_out = self.update_attributes(ds_out, out_grid)

        return ds_out

    def _regrid_variable(self, backend, variable, out_grid, in_dims_arg, out_dims, **kwargs):
        """Regrid a single variable using ``xarray.apply_ufunc``.

        Parameters
        ----------
        backend : Any
            The backend providing the ``regrid`` method.
        variable : Variable
            The variable to regrid, as produced by :func:`~.loader.variables`.
        out_grid : Any
            The output grid.
        in_dims_arg : Optional[List[str]]
            Fallback input grid dimension names, used when they cannot be
            determined from ``variable``.
        out_dims : List[str]
            The output grid dimension names.
        **kwargs
            Extra keyword arguments forwarded to
            ``XarrayDataHandler.regrid``.

        Returns
        -------
        Tuple[xr.DataArray, Any]
            The regridded ``xr.DataArray`` and the (possibly updated) output
            grid spec as returned by the point-to-point regrid method.
        """

        # regrid() can change the specified output gridspec.
        # This is a workaround to get the returned output gridscpec from regrid().
        class _RegridMethod:
            """Callable wrapping ``backend.regrid`` that records the returned output grid.

            ``xarray.apply_ufunc`` only returns the regridded values, so this
            wrapper is used to also capture the (possibly backend-adjusted)
            output grid spec via its ``out_grid`` attribute, updated on each
            call.
            """

            def __init__(self, in_grid, out_grid, **kwargs):
                """Initialise the wrapper.

                Parameters
                ----------
                in_grid : Any
                    The input grid, bound to ``backend.regrid``.
                out_grid : Any
                    The output grid, bound to ``backend.regrid`` and used as
                    the initial value of :attr:`out_grid`.
                **kwargs : dict
                    Extra keyword arguments bound to ``backend.regrid``.
                """
                self.out_grid = out_grid
                self.method = functools.partial(
                    backend.regrid,
                    in_grid=in_grid,
                    out_grid=out_grid,
                    **kwargs,
                )

            def __call__(self, vals):
                """Regrid ``vals``, updating :attr:`out_grid` with the backend's result.

                Parameters
                ----------
                vals : numpy.ndarray
                    The values to regrid.

                Returns
                -------
                numpy.ndarray
                    The regridded values.
                """
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

        method = _RegridMethod(in_grid, out_grid, **kwargs)

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
                "output_sizes": {dim: out_grid.shape[i] for i, dim in enumerate(out_dims)},
                "allow_rechunk": True,
            },
            output_dtypes=[var.dtype],
            keep_attrs="identical",
        )

        res_out_grid = method.out_grid.copy() if method.out_grid is not None else method.out_grid
        return res, res_out_grid


handler = XarrayDataHandler
