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
an ``xarray.Dataset`` or ``xarray.DataArray``. It uses
:mod:`~earthkit.geo.grids.xarray.loader` to identify the geographical
variables and their input grid, builds the output geography via
:class:`XarrayGeographyBuilder`, and regrids each variable
with ``xarray.apply_ufunc`` (dask-aware) delegating the actual point-to-point
interpolation to ``..numpy.NumpyDataHandler``.
"""

import functools
import logging

from earthkit.geo.utils import ensure_list

from ..handler import DataHandler
from ..utils import create_grid_object, scatter_projected_grid

LOG = logging.getLogger(__name__)


def _cf_indexed_message(names) -> str:
    """The error message for data variables numbered by a reduced Gaussian index the data cannot place.

    Parameters
    ----------
    names : List
        The names of the variables which are indexed by a CF
        ``reduced_gaussian_index`` coordinate whose grid mapping and points per
        latitude the dataset does not hold (decision A7).

    Returns
    -------
    str
        A message naming the way out: the dataset which holds the tables, or the
        explicit ``in_grid`` and ``in_dims``.
    """
    listed = ", ".join(f"'{name}'" for name in names)
    return (
        f"Variable {listed} is indexed by a CF reduced Gaussian point index, but the data does not hold "
        "the grid mapping and the table of points per latitude which place the points it numbers. Pass "
        "the dataset which holds them, or give an 'in_grid' and the point dimension in 'in_dims'."
    )


def _no_grid_message(ds, in_grid, in_dims) -> str:
    """The error message for a regrid input which holds no regriddable variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset the variables were looked for in (a DataArray input already
        turned into one by :meth:`XarrayDataHandler.regrid`).
    in_grid : Any
        The user's ``in_grid``, or None when no grid was given.
    in_dims : Any
        The user's ``in_dims``, or None when no dimensions were given.

    Returns
    -------
    str
        A message naming the way out: the data which describes the grid, or the
        explicit ``in_grid`` and ``in_dims``.
    """
    from earthkit.geo.grids.xarray.loader import reduced_gaussian_indexed

    names = reduced_gaussian_indexed(ds)

    if names and in_grid is None:
        return _cf_indexed_message(names)

    if in_grid is not None:
        if in_dims is None:
            return (
                "The data does not describe its grid, and no 'in_dims' were given with the 'in_grid': the "
                "values can only be regridded once the dimensions they number the grid points with are named."
            )
        return (
            f"None of the data variables is regridded on the given 'in_grid': no variable is indexed by all "
            f"of the given 'in_dims' {ensure_list(in_dims)}. Check the dimension names, or pass data which "
            "describes its grid."
        )

    return (
        "Could not determine the grid of the data: no variable of the input is described by latitude and "
        "longitude coordinates or by a CF grid mapping earthkit-geo recognises. Pass data which describes "
        "its grid (e.g. the dataset holding the coordinates or the CF mapping), or give an 'in_grid' and "
        "its grid dimensions in 'in_dims'."
    )


def _check_user_in_dims(var, ek_grid, in_dims) -> None:
    """Check the user-named grid dimensions of a variable whose data does not describe them.

    Parameters
    ----------
    var : xr.DataArray
        The variable to regrid.
    ek_grid : Any
        The eckit-geo grid the values number, given explicitly by the caller.
    in_dims : List
        The dimension names taken from the user's ``in_dims`` argument.

    Raises
    ------
    ValueError
        If one of ``in_dims`` does not index the variable, or the values laid
        out along them do not fill the grid: index-only data of an explicitly
        given grid is the user's responsibility, but a size which cannot match
        the grid is a mistake worth naming (decision A8).
    """
    import numpy as np

    missing = [str(d) for d in in_dims if d not in var.sizes]
    if missing:
        raise ValueError(
            f"Variable '{var.name}' does not have the dimension(s) {missing} named by 'in_dims'; "
            f"its dimensions are {list(var.dims)}."
        )

    shape = getattr(ek_grid, "shape", None)
    if shape is not None:
        n_grid = int(np.prod(shape))
        n_values = int(np.prod([var.sizes[d] for d in in_dims]))
        if n_values != n_grid:
            raise ValueError(
                f"Variable '{var.name}' holds {n_values} values along the dimensions {list(in_dims)} named "
                f"by 'in_dims', but the given 'in_grid' has {n_grid} points: 'in_grid' and 'in_dims' do "
                "not describe the same grid."
            )


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
        import xarray as xr

        has_earthkit = False
        if "_earthkit" in ds.attrs:
            has_earthkit = True
        else:
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
            matches the type of ``values``. Chunked (dask-backed) input stays
            lazy: the result is only computed by the caller, and its geography
            already describes the grid the backend produces (see
            :meth:`_regrid_variable`). When no ``out_dims`` are given, the
            geographical dimensions follow the grid the backend returns: they
            are derived from the requested grid and renamed to the geography
            of the returned one if the two classify it differently; a
            requested ``out_dims`` is never renamed.

            A target which is a plain reduced Gaussian grid, global or with an
            area, is encoded with the CF ``reduced_gaussian`` grid mapping: its
            points are numbered by a ``reduced_gaussian_index`` dimension
            coordinate holding their index on the global grid, in that increasing
            order, and a Dataset carries the global latitude lines, the points
            per latitude and the grid mapping alongside the regridded variables
            (decisions D2 and D3). A DataArray holds none of those tables, so it
            keeps the latitude and longitude of every point (decision D4). Any
            other target, and a reduced Gaussian grid the encoding does not
            cover, is written as it was: on a ``values`` dimension with
            per-point coordinates. The ``grid_mapping`` attribute the input
            carried is gone from the result either way (decision A6).

        Raises
        ------
        ValueError
            If the input grid or input dimensions cannot be determined and are not provided, if a
            variable of the input is numbered by a CF reduced Gaussian index which the data cannot
            place and no ``in_grid`` is given, or if the CF encoding of a reduced Gaussian target
            finds a name it needs taken by the regridded data (decision A9).
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

        from earthkit.geo.grids.xarray.loader import reduced_gaussian_indexed
        from earthkit.geo.grids.xarray.loader import variables as get_variables
        from earthkit.geo.grids.xarray.output import XarrayOutputGeographyBuilder

        input_is_dataset = isinstance(ds, xr.Dataset)
        if not input_is_dataset:
            ds = ds.to_dataset()

        variables = get_variables(ds, user_ek_grid=in_grid, user_in_dims=in_dims_arg)

        if not variables:
            raise ValueError(_no_grid_message(ds, in_grid, in_dims_arg))

        if in_grid is None:
            # a variable numbered by a CF reduced Gaussian index whose mapping the
            # dataset cannot resolve is not lost among the variables which did:
            # it sits on the points of a grid nobody described, and the result
            # would silently leave it out (decision A7). One with usable per-point
            # latitudes and longitudes was regridded through them, so it is in the
            # variables found; only those which are in none of them are named. A
            # variable the dataset declares auxiliary is never named: the loader
            # skips it as data whatever tables the dataset holds. An explicitly
            # given 'in_grid' takes the data as it stands (decision A8).
            regridded = {v.name for v in variables}
            unresolved = [name for name in reduced_gaussian_indexed(ds) if name not in regridded]
            if unresolved:
                raise ValueError(_cf_indexed_message(unresolved))

        for v in variables:
            if v.ek_grid is None:
                if in_grid is None:
                    raise ValueError(
                        f"Could not determine grid for variable {v.name} from dataset. Please provide an "
                        f"'in_grid' argument."
                    )
                else:
                    raise ValueError("Could not determine grid from 'in_grid' argument.")
            if not v.geo_dims and in_dims_arg is None:
                # without detected grid dimensions the variable can only be regridded
                # through the ones the caller named; a variable with neither is no
                # geographical variable at all, and the loader did not return it
                raise ValueError(
                    f"Could not determine grid dimensions for variable {v.name} from dataset. Please provide "
                    "an 'in_dims' argument."
                )

        # the output geography builder which can provide the output grid and the output coordinates
        out_geo = XarrayOutputGeographyBuilder(out_grid)

        # The output dimensions are taken from the user or, failing that, derived
        # from the REQUESTED grid; they are matched to the geography of the grid
        # the backend really returns at the end of this method.
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

        # a CF grid mapping names a variable of the dataset which describes it, and
        # the regridded result is not that dataset: the attribute of the input would
        # point at nothing, whatever grid the values are now on (decision A6)
        self.drop_grid_mapping(ds_out)

        # for a DataArray input will return a single DataArray instead of a Dataset
        if not input_is_dataset:
            ds_out = ds_out[list(ds.keys())[0]]

        # The output geography might have changed, so we need to create a new geography builder
        # with the new grid spec
        out_grid = create_grid_object(res_out_grid)
        out_geo = XarrayOutputGeographyBuilder(out_grid)

        # The dimensions to regrid onto were derived above from the REQUESTED grid,
        # while the coordinates are built from the grid the backend really returned.
        # The two can classify their geography differently (e.g. the requested
        # regular mesh becomes, after the backend adjusted it, a grid without
        # distinct latitudes/longitudes), so when the caller did not choose the
        # dimensions, follow the final geography and rename them. A rename only
        # touches the dimension labels and works on a lazy result too; a
        # user-provided out_dims always wins and is never renamed.
        if out_dims_arg is None:
            final_dims = out_geo.geo_dims()
            if final_dims != out_dims:
                if len(final_dims) != len(out_dims):
                    raise ValueError(
                        f"Cannot match the output dimensions {out_dims} of the regridded data to the "
                        f"geography dimensions {final_dims} of the grid the backend returned: they "
                        "differ in number. Please provide an 'out_dims' argument."
                    )
                clash = sorted(set(final_dims) & (set(ds_out.sizes) - set(out_dims)))
                if clash:
                    raise ValueError(
                        f"Cannot rename the output dimensions {out_dims} of the regridded data to "
                        f"{final_dims}: the result already has the dimensions {clash}. Please "
                        "provide an 'out_dims' argument."
                    )
                ds_out = ds_out.rename(dict(zip(out_dims, final_dims)))
        else:
            final_dims = out_dims

        # A plain reduced Gaussian target is written as a CF reduced Gaussian grid,
        # which also says how its values are ordered; anything else, and a reduced
        # Gaussian grid the encoding does not cover (decision A10), keeps the
        # per-point coordinates of the geography (decisions D3 and D4).
        if out_geo.cf_reduced_gaussian is not None and len(final_dims) == 1 and final_dims[0] in ds_out.dims:
            ds_out = out_geo.add_cf_reduced_gaussian(ds_out, final_dims[0], as_dataset=input_is_dataset)
        else:
            ds_out = self.add_geo_coords(ds_out, out_geo.coords())

        ds_out = self.update_attributes(ds_out, out_grid)

        return ds_out

    @staticmethod
    def drop_grid_mapping(ds):
        """Remove the CF ``grid_mapping`` attribute of the regridded variables.

        The attribute names the variable of the input dataset which carries the
        attributes of the grid the values came from. The result is on another grid
        and holds no such variable, so the reference would dangle; a Dataset on a
        grid with a CF encoding gets the mapping of its target back from
        :meth:`~earthkit.geo.grids.xarray.output.XarrayOutputGeographyBuilder.add_cf_reduced_gaussian`
        (decision A6).

        Parameters
        ----------
        ds : xr.Dataset or xr.DataArray
            The regridded data, modified in place.

        Returns
        -------
        xr.Dataset or xr.DataArray
            ``ds``.
        """
        import xarray as xr

        if isinstance(ds, xr.Dataset):
            for var in ds.data_vars.values():
                var.attrs.pop("grid_mapping", None)
        else:
            ds.attrs.pop("grid_mapping", None)

        return ds

    def _regrid_variable(self, backend, variable, out_grid, in_dims_arg, out_dims, **kwargs):
        """Regrid a single variable using ``xarray.apply_ufunc``.

        Parameters
        ----------
        backend : Any
            The backend providing the ``regrid`` method.
        variable : Variable
            The variable to regrid, as produced by
            :func:`~earthkit.geo.grids.xarray.loader.variables`.
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

        Notes
        -----
        A backend may interpolate onto another grid than the one that was asked
        for (see :meth:`~earthkit.geo.grids._regrid.backends.mir.MirBackend.regrid`),
        so the grid spec of the result is the one it reports, not the requested
        one. ``apply_ufunc`` only learns it while running the backend, which for
        a dask-backed variable happens when the graph is computed, long after the
        output sizes of the graph have been fixed. For that case the backend is
        therefore run once here on a dummy slice of the geographical dimensions
        before the graph is built, and the output sizes and the result grid come
        from the spec that probe returns. Eager variables are untouched by this:
        their backend is called exactly as often as before.
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

        if variable.geo_dims is None:
            # the grid dimensions are only what the user named: they must index
            # the variable and hold exactly as many values as the grid has
            # points (index-only data, see the loader fallback and decision A8)
            _check_user_in_dims(var, variable.ek_grid, in_dims)

        # The points of a CF reduced Gaussian sub-area are held in the order of
        # the file's index, which the detected eckit grid of its area may scan
        # in another order: gather the values into the order of the grid the
        # backend is given (decision A2). The loader only hands out a point
        # order for a grid it detected itself, so values of an explicitly given
        # grid are regridded as they stand (decision A8); the indexing stays
        # lazy for a dask-backed variable.
        if variable.point_order is not None:
            if len(in_dims) != 1:
                raise ValueError(
                    f"Variable '{variable.name}' is reordered onto its grid along its single grid "
                    f"dimension, but is regridded along {len(in_dims)} dimensions {in_dims}."
                )
            var = var.isel({in_dims[0]: variable.point_order})

        exclude_dims = set()
        if set(in_dims) == set(out_dims):
            exclude_dims = set(in_dims)

        # Backends cannot be expected to parse a PROJ-source projected input grid
        # (e.g. the regular_xy grid the loader detects from a CF grid mapping), so
        # the handler converts it here into the unstructured lat/lon grid of its
        # points, which every backend reads (see scatter_projected_grid). Only the
        # input grid is converted: the output grid must keep its own shape.
        in_grid = scatter_projected_grid(variable.ek_grid)

        method = _RegridMethod(in_grid, out_grid, **kwargs)

        # With a dask-backed variable apply_ufunc() does not run the backend while
        # it builds the graph, so the grid the backend really produces is still
        # unknown at that point (see the Notes above). Ask for it with one dummy
        # call on a slice of the geographical dimensions only, in the in_dims order
        # the backend is handed the data in. The probe costs one full spatial
        # interpolation of one slice; whether MIR reuses the weights it computed
        # for it depends on its configuration and on the worker process, so
        # correctness must not rely on that reuse.
        lazy = var.chunks is not None
        if lazy:
            import numpy as np

            method(np.zeros(tuple(var.sizes[dim] for dim in in_dims), dtype=var.dtype))

        # The shape the output sizes of the graph are taken from: the grid of the
        # probe when there was one, and the requested grid otherwise, whose shape
        # the backend cannot disagree with before apply_ufunc() has run it.
        if lazy and method.out_grid is not None:
            out_grid_shape = create_grid_object(method.out_grid).shape
        else:
            out_grid_shape = out_grid.shape

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
                "output_sizes": {dim: out_grid_shape[i] for i, dim in enumerate(out_dims)},
                "allow_rechunk": True,
            },
            output_dtypes=[var.dtype],
            keep_attrs="identical",
        )

        # The grid spec the backend reported, which the caller builds the output
        # geography from. Only a plain spec is copied: it is the value the method
        # keeps updating, and a Grid has no copy() to take a defensive copy with.
        res_out_grid = method.out_grid
        if isinstance(res_out_grid, dict):
            res_out_grid = res_out_grid.copy()
        return res, res_out_grid


handler = XarrayDataHandler
