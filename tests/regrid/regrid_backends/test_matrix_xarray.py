# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest

from earthkit.geo import regrid
from earthkit.geo.utils.testing import (
    NO_EKD,  # noqa: E402
    NO_MIR,  # noqa: E402
    compare_dims,
)

if not NO_EKD:
    from earthkit.data import from_source  # noqa


@pytest.mark.matrix_db
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.parametrize(
    "out_grid,dims",
    [
        ({"grid": [10, 10]}, {"step": 2, "latitude": 19, "longitude": 36}),
    ],
)
def test_regrid_matrix_xarray_from_ogg(out_grid, dims):

    ds_in = from_source("sample", "O32_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear", backend="precomputed")

    compare_dims(r, dims, sizes=True)


@pytest.mark.matrix_db
@pytest.mark.skipif(NO_EKD, reason="No earthkit.data available")
@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.parametrize(
    "out_grid,dims",
    [
        ({"grid": [10, 10]}, {"step": 2, "latitude": 19, "longitude": 36}),
    ],
)
def test_regrid_matrix_xarray_from_h_nested(out_grid, dims):

    ds_in = from_source("sample", "H8_nested_t2.grib2").to_fieldlist()
    assert len(ds_in) == 2
    ds = ds_in.to_xarray()

    r = regrid(ds["2t"], out_grid=out_grid, interpolation="linear", backend="precomputed")

    compare_dims(r, dims, sizes=True)


# A Lambert azimuthal equal-area projection on GRS80, described the way CF files
# describe a projection, as in the synthetic dataset of tests/xarray_grid/test_xarray_grid.py
LAEA_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


def _laea_dataset():
    """A small synthetic Dataset on the plane of a LAEA projection.

    Returns
    -------
    xr.Dataset
        The dataset, with its data variable "a" over the 1D ``x``/``y`` axes of
        the projection plane (5 km, ``y`` descending), the 2D
        ``latitude``/``longitude`` coordinates pyproj computes for them, and the
        scalar CF grid mapping variable naming the projection.
    """
    import xarray as xr
    from pyproj import CRS, Transformer

    x_values = 4321000.0 + 5000.0 * np.arange(6)
    y_values = 3210000.0 - 5000.0 * np.arange(4)

    transformer = Transformer.from_crs(CRS.from_proj4(LAEA_PROJ4), "EPSG:4326", always_xy=True)
    xv, yv = np.meshgrid(x_values, y_values)
    lon2d, lat2d = transformer.transform(xv, yv)

    ds = xr.Dataset(
        {
            "a": xr.Variable(
                ["y", "x"],
                np.arange(24, dtype=np.float64).reshape(4, 6),
                attrs={"grid_mapping": "mapping"},
            )
        },
        coords={
            "x": (["x"], x_values, {"standard_name": "projection_x_coordinate", "units": "m"}),
            "y": (["y"], y_values, {"standard_name": "projection_y_coordinate", "units": "m"}),
            "latitude": (["y", "x"], lat2d, {"standard_name": "latitude", "units": "degrees_north"}),
            "longitude": (["y", "x"], lon2d, {"standard_name": "longitude", "units": "degrees_east"}),
        },
    )
    ds["mapping"] = xr.DataArray(
        0.0,
        attrs={"grid_mapping_name": "lambert_azimuthal_equal_area", "proj4_params": LAEA_PROJ4},
    )
    return ds


def test_regrid_matrix_xarray_projected_input_no_weights(tmp_path):
    """A projected xarray input reaches the precomputed backend as its scattered lat/lon grid.

    The xarray regrid handler converts a ``regular_xy`` input grid with a PROJ source
    projection into the unstructured lat/lon grid of its points before calling the
    backend (see :func:`~earthkit.geo.grids._regrid.data.utils.scatter_projected_grid`).
    The precomputed backend then looks that grid up in the (empty, synthetic) matrix
    inventory and raises for the missing weights: the failure is about the inventory,
    not about parsing the ``+proj=...`` grid. The recorded input grid shows the
    conversion happened in the handler: the backend never sees the projected grid.
    """
    import json

    from earthkit.geo.grids._regrid.backends.precomputed import MatrixBackend
    from earthkit.geo.grids._regrid.data.xarray.loader import variables

    class _RecordingBackend(MatrixBackend):
        """Matrix backend recording the input grids its regrid() is called with."""

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.seen_in_grids = []

        def regrid(self, data, **kwargs):
            self.seen_in_grids.append(kwargs["in_grid"])
            return super().regrid(data, **kwargs)

    # an empty, deterministic matrix inventory: no entry can ever be found in it
    with open(tmp_path / "index.json", "w") as f:
        json.dump({"version": 1, "matrix": {}}, f)

    ds = _laea_dataset()

    # the loader detects the projected regular_xy grid that the handler must convert
    assert [v.ek_grid.type for v in variables(ds)] == ["regular_xy"]

    backend = _RecordingBackend(inventory=str(tmp_path))

    with pytest.raises(ValueError, match="No precomputed weights found!"):
        regrid(ds, out_grid={"grid": [10, 10]}, interpolation="nn", backend=backend)

    # the backend was called once, on the unstructured lat/lon grid of the points
    assert len(backend.seen_in_grids) == 1
    in_grid = backend.seen_in_grids[0]
    assert in_grid.type == "unstructured_ll"
    assert in_grid.shape == (24,)
