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
from earthkit.geo.grids import Grid
from earthkit.geo.utils.testing import (
    NO_EKD,  # noqa: E402
    NO_MIR,  # noqa: E402
    compare_global_ll_results,
    get_test_data,  # noqa: E402
    get_test_data_path,  # noqa: E402
)

if not NO_EKD:
    from earthkit.data import from_source  # noqa
    from earthkit.data import create_fieldlist  # noqa


def _create_fieldlist(filename, subfolder="global_0_360", field_type="grib"):
    ds = from_source("url", get_test_data_path(filename, subfolder=subfolder)).to_fieldlist()
    if field_type == "array":
        return ds.to_fieldlist()
    elif field_type == "grib":
        return ds
    else:
        raise ValueError(f"Unknown field type: {field_type}")


@pytest.mark.skipif(NO_MIR, reason="No mir available")
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
def test_regrid_fieldlist_reg_ll(_kwarg, interpolation, field_type):
    ds = _create_fieldlist("5x5.grib", field_type=field_type)

    f_ref = get_test_data(f"out_5x5_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = ds.metadata(["param", "level", "date", "time", "gridType"])

    r = regrid(ds, out_grid={"grid": [10, 10]}, **_kwarg)

    assert len(r) == 1
    assert r[0].shape == (19, 36)
    compare_global_ll_results(r[0].to_numpy(), v_ref, interpolation, rtol=1e-4)
    assert r.metadata(["param", "level", "date", "time", "gridType"]) == metadata_ref

    grid_ref = {"iDirectionIncrementInDegrees": 10.0, "jDirectionIncrementInDegrees": 10.0}
    for f in r:
        for k, v in grid_ref.items():
            assert np.isclose(f.metadata(k), v), k


@pytest.mark.skipif(NO_MIR, reason="No mir available")
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
def test_regrid_fieldlist_gg(_kwarg, interpolation, field_type):
    ds = _create_fieldlist("O32.grib", field_type=field_type)

    f_ref = get_test_data(f"out_O32_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = ds.metadata(["param", "level", "date", "time"])

    r = regrid(ds, out_grid={"grid": [10, 10]}, **_kwarg)

    assert len(r) == 1
    assert r[0].shape == (19, 36)
    compare_global_ll_results(r[0].to_numpy(), v_ref, interpolation, rtol=1e-4)
    assert r.metadata(["param", "level", "date", "time"]) == metadata_ref

    grid_ref = {"iDirectionIncrementInDegrees": 10.0, "jDirectionIncrementInDegrees": 10.0}
    for f in r:
        for k, v in grid_ref.items():
            assert np.isclose(f.metadata(k), v), k

    assert r.metadata("gridType") == ["regular_ll"]


@pytest.mark.skipif(NO_MIR, reason="No mir available")
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
def test_regrid_single_field(_kwarg, interpolation, field_type):
    ds = _create_fieldlist("O32.grib", field_type=field_type)
    field = ds[0]

    f_ref = get_test_data(f"out_O32_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = field.metadata(["param", "level", "date", "time"])

    r = regrid(field, out_grid={"grid": [10, 10]}, **_kwarg)

    assert r.shape == (19, 36)
    compare_global_ll_results(r.to_numpy(), v_ref, interpolation, rtol=1e-4)
    assert r.metadata(["param", "level", "date", "time"]) == metadata_ref

    grid_ref = {"iDirectionIncrementInDegrees": 10.0, "jDirectionIncrementInDegrees": 10.0}
    for k, v in grid_ref.items():
        assert np.isclose(r.metadata(k), v), k

    assert r.metadata("gridType") == "regular_ll"


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
@pytest.mark.download
@pytest.mark.tmp_cache
def test_regrid_single_field_non_grib():
    interpolation = "linear"
    out_grid = {"grid": [10, 10]}
    field_type = "array"

    fl = _create_fieldlist("5x5.grib", field_type=field_type)
    field = fl[0].set({"parameter.variable": "msl"}, labels={"my_label": "my_value"})

    # the field is now decoupled from the original grib message
    assert field.message() is None

    f_ref = get_test_data(f"out_5x5_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]

    r = regrid(field, out_grid=out_grid, interpolation=interpolation)

    assert r.geography.shape() == (19, 36)
    assert np.allclose(r.values, v_ref)
    assert r.get("parameter.variable") == "msl"
    assert r.get("labels.my_label") == "my_value"

    grid_ref = Grid({"grid": [10, 10]}).spec

    assert r.geography.grid_spec() == grid_ref


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
def test_regrid_fieldlist_deprec_grid_kwarg():
    ds = _create_fieldlist("5x5.grib", field_type="grib")

    interpolation = "linear"

    f_ref = get_test_data(f"out_5x5_10x10_{interpolation}.npz")
    v_ref = np.load(f_ref)["arr_0"]
    metadata_ref = ds.metadata(["param", "level", "date", "time", "gridType"])

    r = regrid(ds, grid={"grid": [10, 10]}, interpolation=interpolation)

    assert len(r) == 1
    assert r[0].shape == (19, 36)
    compare_global_ll_results(r[0].to_numpy(), v_ref, interpolation, rtol=1e-4)
    assert r.metadata(["param", "level", "date", "time", "gridType"]) == metadata_ref

    grid_ref = {"iDirectionIncrementInDegrees": 10.0, "jDirectionIncrementInDegrees": 10.0}
    for f in r:
        for k, v in grid_ref.items():
            assert np.isclose(f.metadata(k), v), k


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
# @pytest.mark.parametrize("field_type", ["grib", "array"])
@pytest.mark.parametrize("field_type", ["grib"])
def test_regrid_grib_1_fieldlist_ll_to_points(field_type):
    ds = _create_fieldlist("5x5_multi.grib1", subfolder="grib", field_type=field_type)

    lats = [40.0, 50.0]
    lons = [10.0, 20.0]
    out_grid = {"latitudes": lats, "longitudes": lons}

    r = regrid(ds, out_grid=out_grid, interpolation="nn")

    metadata_ref = ds.get(["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"])
    points_num = 2

    ref_vals = np.array([[288.44410706, 291.83082581], [289.20581055, 277.6784668]])

    assert len(r) == 2

    for i, f in enumerate(r):
        assert f.shape == (points_num,)
        lats_res, lons_res = f.geography.latlons()
        assert np.allclose(lats_res, lats)
        assert np.allclose(lons_res, lons)
        assert np.allclose(f.values, ref_vals[i])

    assert r.get(["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"]) == metadata_ref


@pytest.mark.skipif(NO_MIR, reason="No mir available")
@pytest.mark.skipif(NO_EKD, reason="No access to earthkit-data")
# @pytest.mark.parametrize("field_type", ["grib", "array"])
@pytest.mark.parametrize("field_type", ["grib"])
def test_regrid_grib_2_fieldlist_ll_to_points_round_trip(field_type):
    ds = _create_fieldlist("5x5_multi.grib2", subfolder="grib", field_type=field_type)

    # Regrid from lat-lon to points
    lats = [40.0, 50.0]
    lons = [10.0, 20.0]
    out_grid_points = {"latitudes": lats, "longitudes": lons}

    r = regrid(ds, out_grid=out_grid_points, interpolation="nn")

    metadata_ref = ds.get(["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"])
    points_num = 2

    ref_vals = np.array([[288.44410706, 291.83082581], [289.20581055, 277.6784668]])

    assert len(r) == 2

    for i, f in enumerate(r):
        assert f.shape == (points_num,)
        lats_res, lons_res = f.geography.latlons()
        assert np.allclose(lats_res, lats)
        assert np.allclose(lons_res, lons)
        assert np.allclose(f.values, ref_vals[i])

    assert r.get(["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"]) == metadata_ref

    # Encode each regridded field to a new GRIB message and create a new fieldlist from them
    # TODO: this is currently failing
    # fields = []
    # for f in r:
    #     f_grib = from_source("memory", f.sync().message()).to_fieldlist()[0]
    #     fields.append(f_grib)

    # ds_r = create_fieldlist(fields)

    # assert len(ds_r) == 2
    # for i, f in enumerate(ds_r):
    #     assert f.shape == (points_num,)
    #     # TODO: this is currently failing
    #     # lats_res, lons_res = f.geography.latlons()
    #     # assert np.allclose(lats_res, lats)
    #     # assert np.allclose(lons_res, lons)
    #     assert np.allclose(f.values, ref_vals[i])

    # assert ds_r.get(["parameter.variable", "vertical.level", "time.valid_datetime", "time.step"]) == metadata_ref

    # TODO: this is currently failing
    # Regrid from points back to lat-lon

    # r_ll = regrid(ds_r, in_grid=out_grid_points, out_grid={"grid": [5,5]}, interpolation="nn")

    # assert len(r_ll) == 2
    # for i, f in enumerate(r_ll):
    #     assert f.shape == (5, 5)
    #     lats_res, lons_res = f.geography.latlons()
    #     # assert np.allclose(lats_res, np.linspace(40.0, 50.0, 5)[:, None])
    #     # assert np.allclose(lons_res, np.linspace(10.0, 20.0, 5)[None, :])
    #     # assert np.allclose(f.values, ref_vals[i].reshape(5, 5))
