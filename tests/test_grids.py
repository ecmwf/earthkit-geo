#!/usr/bin/env python3

# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

import pytest

from earthkit.geo.grids import ICON_CH1_EPS
from earthkit.geo.grids import ICON_CH2_EPS
from earthkit.geo.grids import Equirectangular
from earthkit.geo.grids import Grid
from earthkit.geo.grids import HEALPix
from earthkit.geo.grids import Octahedral
from earthkit.geo.grids import ReducedGaussian


class TestEquirectangular:
    """Tests for Equirectangular grid."""

    def test_creation(self):
        grid = Equirectangular(dlon=1.0, dlat=1.0)
        assert grid.dlon == 1.0
        assert grid.dlat == 1.0
        assert grid.NAME == "Equirectangular"

    def test_creation_with_integers(self):
        grid = Equirectangular(dlon=5, dlat=5)
        assert grid.dlon == 5.0
        assert grid.dlat == 5.0

    def test_grid_spec(self):
        grid = Equirectangular(dlon=2.5, dlat=2.5)
        assert grid.grid_spec["grid"] == [2.5, 2.5]

    def test_code_property(self):
        grid = Equirectangular(dlon=1.0, dlat=2.0)
        assert grid.code == "1.0° × 2.0°"


class TestHEALPix:
    """Tests for HEALPix grid."""

    def test_creation_default_order(self):
        grid = HEALPix(nside=8)
        assert grid.nside == 8
        assert grid.order == "ring"
        assert grid.NAME == "HEALPix"

    def test_creation_nested_order(self):
        grid = HEALPix(nside=16, order="nested")
        assert grid.nside == 16
        assert grid.order == "nested"

    def test_grid_spec(self):
        grid = HEALPix(nside=32, order="ring")
        assert grid.grid_spec["grid"] == "H32"
        assert grid.grid_spec["order"] == "ring"

    def test_code_property(self):
        grid = HEALPix(nside=64)
        assert grid.code == "H64"

    @pytest.mark.parametrize("nside", [8, 16, 32, 64, 128])
    def test_various_nside_values(self, nside):
        grid = HEALPix(nside=nside)
        assert grid.nside == nside
        assert grid.code == f"H{nside}"


class TestReducedGaussian:
    """Tests for ReducedGaussian grid."""

    def test_creation(self):
        grid = ReducedGaussian(nlat=32)
        assert grid.nlat == 32
        assert grid.NAME == "Reduced Gaussian"

    def test_grid_spec(self):
        grid = ReducedGaussian(nlat=48)
        assert grid.grid_spec["grid"] == "N48"

    def test_code_property(self):
        grid = ReducedGaussian(nlat=320)
        assert grid.code == "N320"

    @pytest.mark.parametrize("nlat", [32, 48, 64, 128, 320, 640])
    def test_various_nlat_values(self, nlat):
        grid = ReducedGaussian(nlat=nlat)
        assert grid.nlat == nlat
        assert grid.code == f"N{nlat}"


class TestOctahedral:
    """Tests for Octahedral grid."""

    def test_creation(self):
        grid = Octahedral(nlat=32)
        assert grid.nlat == 32
        assert grid.NAME == "Octahedral"

    def test_grid_spec(self):
        grid = Octahedral(nlat=48)
        assert grid.grid_spec["grid"] == "O48"

    def test_code_property(self):
        grid = Octahedral(nlat=1280)
        assert grid.code == "O1280"

    @pytest.mark.parametrize("nlat", [32, 48, 64, 128, 320, 640, 1280])
    def test_various_nlat_values(self, nlat):
        grid = Octahedral(nlat=nlat)
        assert grid.nlat == nlat
        assert grid.code == f"O{nlat}"


class TestICONGrids:
    """Tests for ICON grid types."""

    def test_icon_ch1_eps(self):
        grid = ICON_CH1_EPS()
        assert grid.NAME == "ICON_CH1_EPS"
        assert grid.code == "ICON_CH1_EPS"

    def test_icon_ch2_eps(self):
        grid = ICON_CH2_EPS()
        assert grid.NAME == "ICON_CH2_EPS"
        assert grid.code == "ICON_CH2_EPS"


class TestFactoryFromDict:
    """Tests for Grid.from_dict factory method."""

    def test_from_dict_healpix(self):
        grid = Grid.from_dict({"grid": "H32"})
        assert isinstance(grid, HEALPix)
        assert grid.nside == 32

    def test_from_dict_healpix_with_order(self):
        grid = Grid.from_dict({"grid": "H64", "order": "nested"})
        assert isinstance(grid, HEALPix)
        assert grid.nside == 64
        assert grid.order == "nested"

    def test_from_dict_reduced_gaussian(self):
        grid = Grid.from_dict({"grid": "N320"})
        assert isinstance(grid, ReducedGaussian)
        assert grid.nlat == 320

    def test_from_dict_octahedral(self):
        grid = Grid.from_dict({"grid": "O1280"})
        assert isinstance(grid, Octahedral)
        assert grid.nlat == 1280

    def test_from_dict_equirectangular_list(self):
        grid = Grid.from_dict({"grid": [1, 1]})
        assert isinstance(grid, Equirectangular)
        assert grid.dlon == 1.0
        assert grid.dlat == 1.0

    def test_from_dict_equirectangular_list_different_values(self):
        grid = Grid.from_dict({"grid": [2.5, 5.0]})
        assert isinstance(grid, Equirectangular)
        assert grid.dlon == 2.5
        assert grid.dlat == 5.0

    def test_from_dict_icon_ch1(self):
        grid = Grid.from_dict({"grid": "ICON_CH1_EPS"})
        assert isinstance(grid, ICON_CH1_EPS)

    def test_from_dict_icon_ch2(self):
        grid = Grid.from_dict({"grid": "ICON_CH2_EPS"})
        assert isinstance(grid, ICON_CH2_EPS)

    def test_from_dict_generic_grid(self):
        with pytest.raises(RuntimeError):
            Grid.from_dict({"grid": "unknown_grid"})


class TestFactoryFromString:
    """Tests for Grid.from_string factory method."""

    def test_from_string_healpix(self):
        grid = Grid.from_string("H32")
        assert isinstance(grid, HEALPix)
        assert grid.nside == 32

    def test_from_string_healpix_with_kwargs(self):
        grid = Grid.from_string("H64", order="nested")
        assert isinstance(grid, HEALPix)
        assert grid.nside == 64
        assert grid.order == "nested"

    def test_from_string_reduced_gaussian(self):
        grid = Grid.from_string("N320")
        assert isinstance(grid, ReducedGaussian)
        assert grid.nlat == 320

    def test_from_string_octahedral(self):
        grid = Grid.from_string("O640")
        assert isinstance(grid, Octahedral)
        assert grid.nlat == 640

    def test_from_string_icon_ch1(self):
        grid = Grid.from_string("ICON_CH1_EPS")
        assert isinstance(grid, ICON_CH1_EPS)

    def test_from_string_icon_ch2(self):
        grid = Grid.from_string("ICON_CH2_EPS")
        assert isinstance(grid, ICON_CH2_EPS)

    def test_from_string_generic_grid(self):
        with pytest.raises(RuntimeError):
            Grid.from_string("unknown_grid")


class TestSubclassFactoryMethods:
    """Tests for calling factory methods on subclasses directly."""

    def test_healpix_from_string(self):
        # Calling from_string on a subclass should use that subclass
        grid = HEALPix.from_string("H32")
        assert isinstance(grid, HEALPix)
        assert grid.nside == 32

    def test_healpix_from_dict(self):
        grid = HEALPix.from_dict({"grid": "H32"})
        assert isinstance(grid, HEALPix)
        assert grid.nside == 32

    def test_wrong_code_for_subclass_from_string(self):
        # Using wrong grid code on a subclass should raise TypeError
        with pytest.raises(TypeError, match="resolves to ReducedGaussian, but called on HEALPix"):
            HEALPix.from_string("N320")

    def test_wrong_code_for_subclass_from_dict(self):
        # Using wrong grid spec on a subclass should raise TypeError
        with pytest.raises(TypeError, match="resolves to Octahedral, but called on HEALPix"):
            HEALPix.from_dict({"grid": "O640"})


class TestGridProperties:
    """Tests for common Grid properties and methods."""

    def test_grid_name(self):
        grid = HEALPix(nside=8)
        assert grid.name == "HEALPix"

    def test_grid_spec_property(self):
        grid = ReducedGaussian(nlat=32)
        spec = grid.grid_spec
        assert "grid" in spec
        assert spec["grid"] == "N32"

    def test_size_method_exists(self):
        # Test that size method exists (actual value depends on eckit)
        grid = HEALPix(nside=8)
        assert hasattr(grid, "size")
        assert callable(grid.size)

    def test_shape_property_exists(self):
        # Test that shape property exists (actual value depends on eckit)
        grid = ReducedGaussian(nlat=32)
        assert hasattr(grid, "shape")

    def test_to_latlon_method_exists(self):
        # Test that to_latlon method exists
        grid = Equirectangular(dlon=5, dlat=5)
        assert hasattr(grid, "to_latlon")
        assert callable(grid.to_latlon)

    def test_bounding_box_method_exists(self):
        # Test that bounding_box method exists
        grid = Octahedral(nlat=32)
        assert hasattr(grid, "bounding_box")
        assert callable(grid.bounding_box)

    def test_plot_method_exists(self):
        # Test that plot method exists
        grid = HEALPix(nside=8)
        assert hasattr(grid, "plot")
        assert callable(grid.plot)


class TestTypeConsistency:
    """Tests to ensure type consistency across different creation methods."""

    def test_healpix_consistency(self):
        grid1 = HEALPix(nside=32)
        grid2 = Grid.from_string("H32")
        grid3 = Grid.from_dict({"grid": "H32"})

        assert type(grid1) is type(grid2) is type(grid3)
        assert grid1.nside == grid2.nside == grid3.nside == 32

    def test_reduced_gaussian_consistency(self):
        grid1 = ReducedGaussian(nlat=48)
        grid2 = Grid.from_string("N48")
        grid3 = Grid.from_dict({"grid": "N48"})

        assert type(grid1) is type(grid2) is type(grid3)
        assert grid1.nlat == grid2.nlat == grid3.nlat == 48

    def test_octahedral_consistency(self):
        grid1 = Octahedral(nlat=64)
        grid2 = Grid.from_string("O64")
        grid3 = Grid.from_dict({"grid": "O64"})

        assert type(grid1) is type(grid2) is type(grid3)
        assert grid1.nlat == grid2.nlat == grid3.nlat == 64

    def test_equirectangular_consistency(self):
        grid1 = Equirectangular(dlon=5, dlat=5)
        grid2 = Grid.from_dict({"grid": [5, 5]})

        assert type(grid1) is type(grid2)
        assert grid1.dlon == grid2.dlon == 5.0
        assert grid1.dlat == grid2.dlat == 5.0
