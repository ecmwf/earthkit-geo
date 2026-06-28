# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,  # noqa: F401
    Literal,
    TypeAlias,
    overload,
)

if TYPE_CHECKING:
    import io  # type: ignore[import]

    import xarray  # type: ignore[import]
    from earthkit.data import Field, FieldList  # type: ignore[import]
    from earthkit.geo.grid import Grid  # type: ignore[import]

ArrayLike: TypeAlias = Any


def _is_array(values):  # IGNORE
    import numpy as np

    return isinstance(values, np.ndarray)  # IGNORE


@overload
def regrid(
    data: "xarray.DataArray|xarray.Dataset",
    in_grid: "dict|str|Grid|None" = None,
    out_grid: "dict|str|Grid|None" = None,
    in_dims=None,
    out_dims=None,
    interpolation: str = "linear",
    backend: Literal["mir", "precomputed"] = "mir",
    **kwargs,
) -> "xarray.DataArray|xarray.Dataset": ...


@overload
def regrid(
    data: "FieldList|Field",
    in_grid: "dict|str|Grid|None" = None,
    out_grid: "dict|str|Grid|None" = None,
    interpolation: str = "linear",
    backend: Literal["mir", "precomputed"] = "mir",
    **kwargs,
) -> "FieldList|Field": ...


@overload
def regrid(
    data: "bytes|io.BytesIO",
    in_grid: "dict|str|Grid|None" = None,
    out_grid: "dict|str|Grid|None" = None,
    interpolation: str = "linear",
    backend: Literal["mir", "precomputed"] = "mir",
    **kwargs,
) -> "bytes|io.BytesIO": ...


def regrid(
    data: "FieldList|Field|xarray.DataArray|xarray.Dataset|bytes|io.BytesIO",
    in_grid: "dict|str|Grid|None" = None,
    out_grid: "dict|str|Grid|None" = None,
    interpolation: str = "linear",
    backend: Literal["mir", "precomputed"] = "mir",
    **kwargs,
):
    r"""Regrid the high-level ``data`` object with the given backend.

    Parameters
    ----------
    data: FieldList|Field|xarray.DataArray|xarray.Dataset|io.BytesIO
        The input data to be regridded. The supported data types are as follows:

        - an earthkit-data :py:class:`~earthkit.data.core.fieldlist.FieldList`
        - an earthkit-data :py:class:`~earthkit.data.core.field.Field`
        - an :class:`xarray.DataArray` or :class:`xarray.Dataset`
        - a GRIB message as a :py:class:`bytes` or :class:`io.BytesIO` object

    in_grid: dict|str|Grid|None, optional
        The grid spec describing the input grid. It is only needed when the input grid
        cannot be automatically inferred from the input data. This can be the case for
        Xarray, where at present the grid information can only be accessed via the
        "earthkit.grid_spec" attribute, and if it is missing or does not contain the necessary
        information ``in_grid`` needs to be provided. When ``in_grid`` is provided, it takes
        precedence over the metadata of the input data.
        The supported grids depends on the regridding ``backend``.
    out_grid: dict|str|Grid|None
        The grid spec describing the target grid that ``data`` will be interpolated onto.
        The supported grids depends on the regridding ``backend``.
    interpolation: str
        The interpolation method. Please note not all the interpolation methods support all possible
        grid types. The possible values are as follows:

        - "linear": Finite Element based interpolation with linear base functions with supporting triangular mesh
        - "grid-box-average": input/output grid box intersections interpolation preserving input value
          integrals (conservative interpolation).
        - "nearest-neighbour": choose a nearest neighbouring input point to define output point value

    backend: {"mir", "precomputed"}, default="mir"
        The regridding backend to use. The possible values are as follows:

        - "mir": Use the MIR regridding library. This is the default backend. See details in the
           documentation at :ref:`mir-regrid-high`.
        - "precomputed": Use pre-computed interpolation weights for some supported grid pairs and interpolation
           methods. See details in the documentation at :ref:`precomputed-regrid-high`.

    **kwargs
        Additional keyword arguments depending on data type, interpolation or backend.


    Returns
    -------
    "FieldList|Field|xarray.DataArray|xarray.Dataset|bytes|io.BytesIO"
        The regridded data with the same type as ``data`` but with the grid changed to the output grid.


    Examples
    --------
    Using the "mir" backend:

    - :ref:`/how-tos/mir/mir_healpix_fieldlist.ipynb`
    - :ref:`/how-tos/mir/mir_octahedral_fieldlist.ipynb`
    - :ref:`/how-tos/mir/mir_interpolation_types.ipynb`
    - :ref:`/tutorials/mir/mir_regrid_xarray.ipynb`

    Using the "precomputed" backend:

    - :ref:`/how-tos/precomputed/precomp_healpix_fieldlist.ipynb`
    - :ref:`/how-tos/precomputed/precomp_octahedral_fieldlist.ipynb`

    """
    from earthkit.geo.grids._regrid.data import get_data_handler

    grid = kwargs.pop("grid", None)
    if grid is not None:
        if out_grid is not None:
            raise ValueError("Cannot specify both 'grid' and 'out_grid'")
        import warnings

        warnings.warn("'grid' is deprecated in regrid(). Use 'out_grid' instead", DeprecationWarning)
        out_grid = grid

    h = get_data_handler(data)
    if h is None:
        if _is_array(data):
            txt = f"Unsupported data type={type(data)}. Use earthkit.geo.grids.regrid.array.regrid() for arrays"
        else:
            txt = f"Unsupported type={type(data)}"
        raise ValueError(txt)

    kwargs = kwargs.copy()
    return h.regrid(data, in_grid=in_grid, out_grid=out_grid, interpolation=interpolation, backend=backend, **kwargs)
