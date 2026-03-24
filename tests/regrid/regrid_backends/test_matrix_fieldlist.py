# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import numpy as np
import pytest
from eckit.geo import Grid

from earthkit.geo.regrid import regrid
from earthkit.geo.utils.testing import (
    NO_EKD,  # noqa: E402
    SYSTEM_MATRIX_BACKEND_NAME,  # noqa: E402
    get_test_data,  # noqa: E402
    get_test_data_path,  # noqa: E402
)

if not NO_EKD:
    from earthkit.data import from_source  # noqa


def _create_fieldlist(filename, field_type):
    ds = from_source("url", get_test_data_path(filename)).to_fieldlist()
    if field_type == "array":
        return ds.to_fieldlist()
    elif field_type == "grib":
        return ds
    else:
        raise ValueError(f"Unknown field type: {field_type}")


@pytest.mark.download
@pytest.mark.tmp_cache
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
@pytest.mark.parametrize(
    "_kwarg,interpolation",
    [
        ({}, "linear"),
        ({"interpolation": "linear"}, "linear"),
        ({"interpolation": "nearest-neighbour"}, "nearest-neighbour"),
        ({"interpolation": "nn"}, "nearest-neighbour"),
        ({"interpolation": "nearest-neighbor"}, "nearest-neighbour"),
    ],
)
@pytest.mark.parametrize("field_type", ["grib", "array"])
@pytest.mark.parametrize("out_grid", [{"grid": [10, 10]}, Grid({"grid": [10, 10]}), "'grid': [10, 10]"])
def test_regrid_matrix_fieldlist_reg_ll(_kwarg, interpolation, field_type, out_grid):
    ds = _create_fieldlist("5x5.grib", field_type)

    f_ref = get_test_data(f"out_5x5_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = ds.metadata(["param", "level", "date", "time"])

    r = regrid(ds, grid=out_grid, backend=SYSTEM_MATRIX_BACKEND_NAME, **_kwarg)

    assert len(r) == 1
    assert r[0].geography.shape() == (19, 36)
    assert np.allclose(r[0].values, v_ref)
    assert r.metadata(["param", "level", "date", "time"]) == metadata_ref

    grid_ref = Grid({"grid": [10, 10]}).spec

    for f in r:
        assert f.geography.grid().spec == grid_ref


@pytest.mark.download
@pytest.mark.tmp_cache
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
@pytest.mark.parametrize(
    "_kwarg,interpolation",
    [
        ({}, "linear"),
        ({"interpolation": "linear"}, "linear"),
        ({"interpolation": "nearest-neighbour"}, "nearest-neighbour"),
        ({"interpolation": "nn"}, "nearest-neighbour"),
        ({"interpolation": "nearest-neighbor"}, "nearest-neighbour"),
    ],
)
@pytest.mark.parametrize("field_type", ["grib", "array"])
@pytest.mark.parametrize("out_grid", [{"grid": [10, 10]}, Grid({"grid": [10, 10]}), "'grid': [10, 10]"])
def test_regrid_matrix_fieldlist_gg(_kwarg, interpolation, field_type, out_grid):
    ds = _create_fieldlist("O32.grib", field_type)

    f_ref = get_test_data(f"out_O32_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = ds.metadata(["param", "level", "date", "time"])

    r = regrid(ds, grid=out_grid, backend=SYSTEM_MATRIX_BACKEND_NAME, **_kwarg)

    assert len(r) == 1
    assert r[0].geography.shape() == (19, 36)
    assert np.allclose(r[0].values, v_ref)
    assert r.metadata(["param", "level", "date", "time"]) == metadata_ref

    grid_ref = Grid({"grid": [10, 10]}).spec

    for f in r:
        assert f.geography.grid().spec == grid_ref


@pytest.mark.download
@pytest.mark.tmp_cache
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
@pytest.mark.parametrize(
    "_kwarg,interpolation",
    [
        ({}, "linear"),
        ({"interpolation": "linear"}, "linear"),
        ({"interpolation": "nearest-neighbour"}, "nearest-neighbour"),
        ({"interpolation": "nn"}, "nearest-neighbour"),
        ({"interpolation": "nearest-neighbor"}, "nearest-neighbour"),
    ],
)
@pytest.mark.parametrize("field_type", ["grib", "array"])
@pytest.mark.parametrize("out_grid", [{"grid": [10, 10]}, Grid({"grid": [10, 10]}), "'grid': [10, 10]"])
def test_regrid_matrix_single_field(_kwarg, interpolation, field_type, out_grid):
    ds = _create_fieldlist("5x5.grib", field_type)

    f_ref = get_test_data(f"out_5x5_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    field = ds[0]
    metadata_ref = field.metadata(["param", "level", "date", "time"])

    r = regrid(field, grid=out_grid, backend=SYSTEM_MATRIX_BACKEND_NAME, **_kwarg)

    assert r.geography.shape() == (19, 36)
    assert np.allclose(r.values, v_ref)
    assert r.metadata(["param", "level", "date", "time"]) == metadata_ref

    grid_ref = Grid({"grid": [10, 10]}).spec

    assert r.geography.grid().spec == grid_ref
