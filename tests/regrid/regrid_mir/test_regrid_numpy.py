# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Numpy array regridding tests (``grids.array.regrid``) against MIR reference files.

The shared pieces -- interpolation name lists, the kwarg alias parametrisation, the
``mir`` skip mark and the loading of the legacy npz data -- come from
``regrid_common``. The reference comparisons deliberately keep the strict
``compare_global_ll_results`` helper of ``earthkit.geo.utils.testing`` or bare numpy
assertions instead of the registry compare ``regrid_common.get_interpolation(...).compare``,
which is used by the grid matrix tests (``test_regrid_numpy_grids.py``). These legacy
checks keep the exact tolerance, pole-row and tie handling each test already had (the
npz references of ``global_0_360`` are shared with other test folders). The per-test
comments say which.
"""

import numpy as np
import pytest
from regrid_common import (
    BASE_INTERPOLATION_NAMES,
    INTERPOLATION_NAMES,
    KWARG_ALIASES,
    legacy_values,
    requires_mir,
)

from earthkit.geo.grids.array import regrid as regrid_array
from earthkit.geo.utils.testing import compare_global_ll_results, proportion_close

# Every spelling of an interpolation MIR accepts in these tests: the canonical names
# from the registry plus four further MIR spellings (aliases / internal names).
MIR_INTERPOLATION_SPELLINGS = INTERPOLATION_NAMES + [
    "grid-box-statistics",
    "nn",
    "bilinear",
    "k-nearest-neighbours",
]

LATLON_REFS_FULL = [
    (0.05, (3601, 7200)),
    (0.1, (1801, 3600)),
    (0.125, (1441, 2880)),
    (0.15, (1201, 2400)),
    (0.2, (901, 1800)),
    (0.25, (721, 1440)),
    (0.3, (601, 1200)),
    (0.35, (515, 1029)),
    (0.4, (451, 900)),
    (0.5, (361, 720)),
    (0.6, (301, 600)),
    (0.7, (257, 515)),
    (0.8, (225, 450)),
    (0.9, (201, 400)),
    (1, (181, 360)),
    (1.2, (151, 300)),
    (1.25, (145, 288)),
    (1.4, (129, 258)),
    (1.5, (121, 240)),
    (1.6, (113, 225)),
    (1.8, (101, 200)),
    (2, (91, 180)),
    (2.5, (73, 144)),
    (5, (37, 72)),
    (10, (19, 36)),
]

LATLON_REFS = [
    (0.5, (361, 720)),
]

# Larger resolutions (O1280, N320) and the ORCA grids below are too slow for now.
RED_GG_O_REFS = [
    ("O48", (10944,)),
]

RED_GG_N_REFS = [
    ("N48", (13280,)),
]

REG_GG_REFS = [
    ("F32", (64, 128)),
]

H_RING_REFS = [
    ("H4", (192,)),
]

H_NESTED_REFS = [
    ("H4", (192,)),
]

ORCA_REFS = []  # e.g. ("eORCA025_T", (1442, 1207)) -- too slow for now


def make_refs():
    """Create references for the tests."""
    refs = []
    for res_dx, res_shape in LATLON_REFS:
        refs.append(({"grid": [res_dx, res_dx]}, res_shape))

    for out_grid, res_shape in RED_GG_O_REFS:
        refs.append(({"grid": out_grid}, res_shape))

    for out_grid, res_shape in RED_GG_N_REFS:
        refs.append(({"grid": out_grid}, res_shape))

    for out_grid, res_shape in REG_GG_REFS:
        refs.append(({"grid": out_grid}, res_shape))

    for out_grid, res_shape in H_RING_REFS:
        refs.append(({"grid": out_grid, "order": "ring"}, res_shape))

    for out_grid, res_shape in H_NESTED_REFS:
        refs.append(({"grid": out_grid, "order": "nested"}, res_shape))

    for out_grid, res_shape in ORCA_REFS:
        refs.append(({"grid": out_grid}, res_shape))

    return refs


REFS = make_refs()

# Input grids of the any-to-any / any-to-points tests (grid spec, shape of the input
# values). The eORCA025_T grid (1442, 1207) is left out: this test is too slow for now.
IN_GRID_SHAPES = [
    ({"grid": [5, 5]}, (37, 72)),
    ({"grid": "O32"}, (5248,)),
    ({"grid": "N32"}, (6114,)),
    ({"grid": "F48"}, (96, 192)),
    ({"grid": "H8", "order": "ring"}, (768,)),
    ({"grid": "H8", "order": "nested"}, (768,)),
]


@requires_mir
# KWARG_ALIASES has no grid-box-average entry: there is no out_O32_10x10_grid-box-average
# reference file in the legacy global_0_360 folder, so that alias was already commented
# out in the old parametrisation.
@pytest.mark.parametrize("_kwarg,interpolation", KWARG_ALIASES)
def test_regrid_numpy_kwarg(_kwarg, interpolation):
    v_in = legacy_values("in_O32.npz")
    v_ref = legacy_values(f"out_O32_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(v_in, {"grid": "O32"}, {"grid": [10, 10]}, **_kwarg)

    assert v_res.shape == (19, 36)
    # Strict comparison kept: the registry compare would let up to 1% of the nearest
    # values differ, this reference check is exact except for the pole rows.
    compare_global_ll_results(v_res, v_ref, interpolation, rtol=1e-4)


def _ll_to_ll(interpolation):
    v_in = legacy_values("in_5x5.npz")
    v_ref = legacy_values(f"out_5x5_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(v_in, {"grid": [5, 5]}, {"grid": [10, 10]}, interpolation=interpolation)

    assert v_res.shape == (19, 36)
    # Strict comparison kept: np.allclose over every value, stricter than the registry.
    assert np.allclose(v_res.flatten(), v_ref)

    # repeated use
    v_res, _ = regrid_array(v_in, {"grid": [5, 5]}, {"grid": [10, 10]}, interpolation=interpolation)

    assert v_res.shape == (19, 36)
    assert np.allclose(v_res.flatten(), v_ref)


@requires_mir
@pytest.mark.parametrize("interpolation", BASE_INTERPOLATION_NAMES)
def test_regrid_numpy_ll_to_ll_1(interpolation):
    _ll_to_ll(interpolation)


@requires_mir
@pytest.mark.parametrize("interpolation", ["linear"])
@pytest.mark.parametrize("res_dx,res_shape", LATLON_REFS_FULL)
def test_regrid_numpy_ll_to_ll_2(res_dx, res_shape, interpolation):
    values = np.random.random((37, 72))
    res_v, _ = regrid_array(
        values,
        in_grid={"grid": [5, 5]},
        out_grid={"grid": [res_dx, res_dx]},
        interpolation=interpolation,
    )
    assert res_v.shape == res_shape


@requires_mir
@pytest.mark.parametrize("interpolation", BASE_INTERPOLATION_NAMES)
def test_regrid_numpy_ogg_to_ll_1(interpolation):
    v_in = legacy_values("in_O32.npz")
    v_ref = legacy_values(f"out_O32_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(v_in, {"grid": "O32"}, {"grid": [10, 10]}, interpolation=interpolation)

    assert v_res.shape == (19, 36)
    # Strict comparison kept, see test_regrid_numpy_kwarg.
    compare_global_ll_results(v_res, v_ref, interpolation, rtol=1e-4)


@requires_mir
@pytest.mark.parametrize("interpolation", INTERPOLATION_NAMES)
@pytest.mark.parametrize("in_grid", ["O32", "o32"])
def test_regrid_numpy_ogg_to_ll_2(interpolation, in_grid):
    values_in = legacy_values("in_O32.npz")

    values, gridspec = regrid_array(
        values_in,
        {"grid": in_grid},
        {"grid": [30, 30]},
        interpolation=interpolation,
    )

    assert values.shape == (7, 12)
    assert gridspec == dict(grid=[30, 30])


@requires_mir
@pytest.mark.parametrize("interpolation", MIR_INTERPOLATION_SPELLINGS)
def test_regrid_numpy_ogg_to_ll_interpolations(interpolation):
    # O32
    values = np.random.random(5248)
    res_v, _ = regrid_array(
        values,
        in_grid={"grid": "O32"},
        out_grid={"grid": [10, 10]},
        interpolation=interpolation,
    )
    assert res_v.shape == (19, 36), f"{interpolation=} failed"


@requires_mir
@pytest.mark.parametrize("interpolation", BASE_INTERPOLATION_NAMES)
def test_regrid_numpy_ngg_to_ll(interpolation):
    v_in = legacy_values("in_N32.npz")
    v_ref = legacy_values(f"out_N32_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(v_in, {"grid": "N32"}, {"grid": [10, 10]}, interpolation=interpolation)

    assert v_res.shape == (19, 36)

    # Compare only the valid values (FIXME inaccurate reference results, ingnore first and last row)
    mask = np.ones_like(v_res, dtype=bool)
    mask[0] = False
    mask[-1] = False
    # Strict comparison kept (default tolerance of compare_global_ll_results): the
    # registry compare would let up to 1% of the nearest values differ.
    compare_global_ll_results(v_res[mask], v_ref[mask.flatten()], interpolation)


@requires_mir
@pytest.mark.parametrize("interpolation", BASE_INTERPOLATION_NAMES)
@pytest.mark.parametrize("in_grid", ["H4", "h4"])
def test_regrid_numpy_healpix_ring_to_ll(interpolation, in_grid):
    v_in = legacy_values("in_H4_ring.npz")
    v_ref = legacy_values(f"out_H4_ring_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(
        v_in,
        {"grid": in_grid, "order": "ring"},
        {"grid": [10, 10]},
        interpolation=interpolation,
    )

    assert v_res.shape == (19, 36)

    # Strict comparisons kept (see the module docstring): the old rules apply to the
    # whole flattened result, while the registry compare exempts the pole rows.
    # nearest-neighbour can produce machine-dependent results, so allow for a few
    # values to be different
    if interpolation == "nearest-neighbour":
        allowable_prop_close = 0.99
        prop_close = proportion_close(v_res.flatten(), v_ref, allowable_prop_close)
        assert prop_close >= allowable_prop_close, f"Only {prop_close:.1%} values are close"
    else:
        np.testing.assert_allclose(v_res.flatten(), v_ref, verbose=False)


@requires_mir
@pytest.mark.skip(reason="This test is currently failing, the reference has to be updated")
@pytest.mark.parametrize("interpolation", BASE_INTERPOLATION_NAMES)
def test_regrid_numpy_healpix_nested_to_ll(interpolation):
    v_in = legacy_values("in_H4_nested.npz")
    v_ref = legacy_values(f"out_H4_nested_10x10_{interpolation}.npz")

    v_res, _ = regrid_array(
        v_in,
        {"grid": "H4", "order": "nested"},
        {"grid": [10, 10]},
        interpolation=interpolation,
    )

    assert v_res.shape == (19, 36)
    v_ref = v_ref.reshape(v_res.shape)

    # Whole-array assert_allclose, stricter than the registry compare (see the module
    # docstring), kept as it was.
    np.testing.assert_allclose(v_res.flatten(), v_ref.flatten(), verbose=False)


@requires_mir
@pytest.mark.skip(reason="This test is too slow for now")
@pytest.mark.parametrize("interpolation", ["linear"])
def test_regrid_numpy_orca_to_ogg(interpolation):
    v_in = legacy_values("in_eORCA025_T.npz", subfolder="orca")
    v_ref = legacy_values(f"out_eORCA025_T_O96_{interpolation}.npz", subfolder="local")

    out_grid = {"grid": "O96"}
    v_res, grid_res = regrid_array(
        v_in,
        {
            "grid": "eORCA025_T",
        },
        out_grid,
        interpolation=interpolation,
    )

    assert v_res.shape == (40320,)
    assert grid_res == out_grid

    # Compare only the valid values (FIXME inaccurate reference results)
    # (assert not np.any(np.isnan(v_res) ^ np.isnan(v_ref)),
    #  "Mismatch of missing values between result and reference")
    mask = ~np.isnan(v_res) & ~np.isnan(v_ref)
    np.testing.assert_allclose(v_res[mask], v_ref[mask], atol=5e-1, rtol=5e-3, verbose=False)


@requires_mir
@pytest.mark.parametrize("interpolation", ["linear"])
@pytest.mark.parametrize("in_grid,in_shape", IN_GRID_SHAPES)
@pytest.mark.parametrize("out_grid,res_shape", REFS)
def test_regrid_numpy_any_to_any(interpolation, in_grid, in_shape, out_grid, res_shape):
    values = np.random.random(in_shape)
    res_v, _ = regrid_array(values, in_grid=in_grid, out_grid=out_grid, interpolation=interpolation)
    assert res_v.shape == res_shape, f"Expected shape {res_shape}, got {res_v.shape}"


@requires_mir
@pytest.mark.parametrize("interpolation", ["nn"])
@pytest.mark.parametrize("in_grid,in_shape", IN_GRID_SHAPES)
def test_regrid_numpy_any_to_points(interpolation, in_grid, in_shape):
    values = np.random.random(in_shape)

    lats = [40.0, 50.0]
    lons = [10.0, 20.0]
    out_grid = {"latitudes": lats, "longitudes": lons}
    res_v, _ = regrid_array(values, in_grid=in_grid, out_grid=out_grid, interpolation=interpolation)
    assert res_v.shape == (len(lats),), f"Expected shape {(len(lats), len(lons))}, got {res_v.shape}"


@requires_mir
def test_regrid_numpy_rotation_not_supported():
    v_in = np.ones((37, 72))
    out_grid = {"grid": [10, 10], "rotation": [-43, 10]}
    with pytest.raises(ValueError, match="only supported for GRIB"):
        regrid_array(v_in, {"grid": [5, 5]}, out_grid)
