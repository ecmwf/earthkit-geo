# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

import pytest

from earthkit.geo import Grid

ARAKAWA_C_UM_N96_T_SPEC = {"grid": [1.875, 1.25], "order": "i+j+", "reference": [0.9375, 0.625]}


class TestArakawaCUMConstruction:
    """Equivalent constructions of the arakawa_c_um N=96 T-arrangement grid."""

    def test_kwargs(self):
        assert Grid(type="arakawa_c_um", N=96).spec == ARAKAWA_C_UM_N96_T_SPEC

    def test_yaml_string_explicit_arrangement_T(self):
        assert Grid("{type: arakawa_c_um, N: 96, arrangement: T}").spec == ARAKAWA_C_UM_N96_T_SPEC

    def test_yaml_string_arakawa_c_with_grid_factor(self):
        # grid-factor (hyphenated key) is only expressible via the YAML string form
        g = Grid("{type: arakawa_c, N: 96, grid-factor: [2, 1.3333333333], order: i+j+}")
        assert g.spec == ARAKAWA_C_UM_N96_T_SPEC

    def test_yaml_string_explicit_grid_and_reference(self):
        g = Grid("{grid: [1.875, 1.25], reference: [0.9375, 0.625], order: i+j+}")
        assert g.spec == ARAKAWA_C_UM_N96_T_SPEC
        assert g.shape == (144, 192)

    def test_yaml_string_explicit_grid_and_area(self):
        # area is resolved to an equivalent reference in the spec
        g = Grid("{grid: [1.875, 1.25], area: [89.375, 0.9375, -89.375, 359.0625], order: i+j+}")
        assert g.spec == ARAKAWA_C_UM_N96_T_SPEC


class TestArakawaCUMEquality:
    """Grid equality across construction methods."""

    def test_all_equivalent_constructions_are_equal(self):
        a = Grid(type="arakawa_c_um", N=96)
        b = Grid("{type: arakawa_c, N: 96, grid-factor: [2, 1.3333333333], order: i+j+}")
        c = Grid("{grid: [1.875, 1.25], reference: [0.9375, 0.625], order: i+j+}")
        d = Grid("{grid: [1.875, 1.25], area: [89.375, 0.9375, -89.375, 359.0625], order: i+j+}")
        assert a == b == c == d

    def test_kwargs_grid_factor_precision_matters(self):
        # grid_factor=[2, 1.33333333] (8 decimal places) differs from the exact fraction
        a = Grid(type="arakawa_c_um", N=96)
        b = Grid(type="arakawa_c", N=96, grid_factor=[2, 1.33333333], order="i+j+")
        assert a != b

    def test_arrangement_T_equals_default(self):
        a = Grid(type="arakawa_c_um", N=96)
        e = Grid("{type: arakawa_c_um, N: 96, arrangement: T}")
        assert a == e

    def test_arrangement_U_differs_from_default(self):
        a = Grid(type="arakawa_c_um", N=96)
        e = Grid("{type: arakawa_c_um, N: 96, arrangement: U}")
        assert a != e


class TestArakawaCUMArrangements:
    """Arrangement matching is case-insensitive."""

    @pytest.mark.parametrize("arrangement", ["T", "t", "U", "u", "V", "v"])
    def test_arrangement_case_insensitive(self, arrangement):
        Grid(f"{{type: arakawa_c_um, N: 96, arrangement: {arrangement}}}")


class TestArakawaCUMLatLons:
    """to_latlons() point layout matches C++ to_points() for each Arakawa-C arrangement."""

    EPS = 1e-3

    def test_arrangement_T(self):
        g = Grid(type="arakawa_c_um", N=96, arrangement="T")
        lat, lon = g.to_latlons()
        shape = (144, 192)
        size = shape[0] * shape[1]
        assert len(lat) == size
        assert lon[0] == pytest.approx(0.9375)
        assert lat[0] == pytest.approx(-89.375)
        assert lon[1] == pytest.approx(2.8125)
        assert lat[1] == pytest.approx(-89.375)
        assert lon[-2] == pytest.approx(357.188, abs=self.EPS)
        assert lat[-2] == pytest.approx(89.375, abs=self.EPS)
        assert lon[-1] == pytest.approx(359.062, abs=self.EPS)
        assert lat[-1] == pytest.approx(89.375, abs=self.EPS)

    def test_arrangement_U(self):
        g = Grid(type="arakawa_c_um", N=96, arrangement="U")
        lat, lon = g.to_latlons()
        shape = (144, 192)
        size = shape[0] * shape[1]
        assert len(lat) == size
        assert lon[0] == pytest.approx(0.0)
        assert lat[0] == pytest.approx(-89.375)
        assert lon[1] == pytest.approx(1.875)
        assert lat[1] == pytest.approx(-89.375)
        assert lon[-2] == pytest.approx(356.25, abs=self.EPS)
        assert lat[-2] == pytest.approx(89.375, abs=self.EPS)
        assert lon[-1] == pytest.approx(358.125, abs=self.EPS)
        assert lat[-1] == pytest.approx(89.375, abs=self.EPS)

    def test_arrangement_V(self):
        g = Grid(type="arakawa_c_um", N=96, arrangement="V")
        lat, lon = g.to_latlons()
        shape = (144, 192)
        size = (shape[0] + 1) * shape[1]
        assert len(lat) == size
        assert lon[0] == pytest.approx(0.9375)
        assert lat[0] == pytest.approx(-90.0)
        assert lon[1] == pytest.approx(2.8125)
        assert lat[1] == pytest.approx(-90.0)
        assert lon[-2] == pytest.approx(357.188, abs=self.EPS)
        assert lat[-2] == pytest.approx(90.0, abs=self.EPS)
        assert lon[-1] == pytest.approx(359.062, abs=self.EPS)
        assert lat[-1] == pytest.approx(90.0, abs=self.EPS)
