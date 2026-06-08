# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


def _is_array(values):  # IGNORE
    import numpy as np

    return isinstance(values, np.ndarray)  # IGNORE


def regrid(data, in_grid=None, out_grid=None, *, interpolation="linear", backend="mir", **kwargs):
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
