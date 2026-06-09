# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
#


from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Literal,
)

if TYPE_CHECKING:
    from earthkit.data.grids import Grid  # type: ignore[import]
    from numpy._typing import NDArray  # type: ignore[import]


def regrid(
    data: "NDArray",
    in_grid: "dict|str|Grid|None" = None,
    out_grid: "dict|str|Grid|None" = None,
    *,
    interpolation: str = "linear",
    backend: Literal["mir", "precomputed"] = "mir",
    **kwargs,
):
    r"""Array interface.

    Regrid the high-level ``data`` object with the given backend.

    Parameters
    ----------
    data: NDArray
        The input data to be regridded, represented as a NumPy array, defining a
        single field on the ``in_grid``.
    in_grid: dict|str|Grid
        The grid spec describing the input grid. The supported grids depends on the
        regridding ``backend``.
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
           documentation at :ref:`mir-regrid-array`.
        - "precomputed": Use pre-computed interpolation weights for some supported grid pairs and interpolation
           methods. See details in the documentation at :ref:`precomputed-regrid-array`.

    **kwargs
        Additional keyword arguments depending on data type, interpolation or backend.


    Returns
    -------
    NDArray, dict
        Tuple of the regridded data as a NumPy array and the output grid spec as a dictionary. The output
        grid spec can be different from the input ``out_grid`` argument, depending on the regridding backend
        since normalisation and other adjustments may be applied during the regridding process.

    Examples
    --------
    Using the "mir" backend:

    - :ref:`/how-tos/mir/mir_numpy_array.ipynb`

    Using the "precomputed" backend:

    - :ref:`/how-tos/precomputed/precomp_numpy_array.ipynb`
    - :ref:`/how-tos/precomputed/precomp_.ipynb`

    """
    from earthkit.geo.grids._regrid.data.numpy import handler

    h = handler()
    kwargs = kwargs.copy()
    return h.regrid(data, in_grid=in_grid, out_grid=out_grid, interpolation=interpolation, backend=backend, **kwargs)
