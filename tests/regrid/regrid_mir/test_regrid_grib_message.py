# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Regridding a GRIB message (a buffer, not an earthkit-data field) with MIR.

The interpolation name list and the ``mir`` / ``earthkit-data`` skip marks come from
``regrid_common``. The input is the tracked ``tests/data/o32.grib2``, read through
``earthkit_test_data_path`` (never downloaded, so not a ``regrid_common.data_path`` case).
"""

import warnings
from io import BytesIO

import numpy as np
import pytest
from regrid_common import INTERPOLATION_NAMES, requires_ekd, requires_mir

from earthkit.geo import regrid
from earthkit.geo.utils.testing import earthkit_test_data_path


@requires_mir
@requires_ekd
@pytest.mark.parametrize("input_format", ["BytesIO", "bytes"])
@pytest.mark.parametrize("interpolation", INTERPOLATION_NAMES)
def test_regrid_grib_message_to_ll(input_format, interpolation):
    from earthkit.data.field.grib.create import create_grib_field_from_buffer

    def check_header(contents):
        assert contents.startswith(b"GRIB") and contents.endswith(b"7777")

    with open(earthkit_test_data_path("o32.grib2"), "rb") as fh:
        in_grib = BytesIO(fh.read())
        check_header(in_grib.getvalue())
        if input_format == "bytes":
            in_grib = in_grib.getvalue()

    out = regrid(in_grib, out_grid={"grid": [30, 30]}, interpolation=interpolation)

    field = create_grib_field_from_buffer(out)
    assert field.metadata("gridType") == "regular_ll"
    assert field.shape == (7, 12)


@requires_mir
@requires_ekd
def test_regrid_grib_message_to_rotated_ll():
    from earthkit.data.field.grib.create import create_grib_field_from_buffer

    with open(earthkit_test_data_path("o32.grib2"), "rb") as fh:
        in_grib = fh.read()

    out_grid = {"grid": [1, 1], "area": [10, -10, -10, 10], "rotation": [-43, 10]}
    with pytest.warns(DeprecationWarning):
        out = regrid(in_grib, out_grid=out_grid, interpolation="nearest-neighbour")

    field = create_grib_field_from_buffer(out)
    assert field.metadata("gridType") == "rotated_ll"
    assert field.metadata("latitudeOfSouthernPoleInDegrees") == -43
    assert field.metadata("longitudeOfSouthernPoleInDegrees") == 10
    assert field.shape == (21, 21)


def _regrid_o32(out_grid, interpolation="nearest-neighbour", data=None):
    """Regrid the tracked o32.grib2 and return the decoded output field."""
    from earthkit.data.field.grib.create import create_grib_field_from_buffer

    if data is None:
        with open(earthkit_test_data_path("o32.grib2"), "rb") as fh:
            data = fh.read()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        out = regrid(data, out_grid=out_grid, interpolation=interpolation)
    return create_grib_field_from_buffer(out)


@requires_mir
@requires_ekd
@pytest.mark.parametrize("interpolation", ["nearest-neighbour", "linear"])
def test_regrid_grib_message_identity_rotation(interpolation):
    # A south pole at [-90, 0] is the identity rotation
    base = {"grid": [1, 1], "area": [10, -10, -10, 10]}
    ref = _regrid_o32(base, interpolation)
    rot = _regrid_o32({**base, "rotation": [-90, 0]}, interpolation)

    assert ref.metadata("gridType") == "regular_ll"
    assert rot.metadata("gridType") == "rotated_ll"
    np.testing.assert_array_equal(rot.to_numpy(), ref.to_numpy())


@requires_mir
@requires_ekd
@pytest.mark.parametrize("interpolation", ["nearest-neighbour", "linear"])
def test_regrid_grib_message_rotation_global(interpolation):
    field = _regrid_o32({"grid": [10, 10], "rotation": [-43, 10]}, interpolation)

    assert field.metadata("gridType") == "rotated_ll"
    assert field.metadata("latitudeOfSouthernPoleInDegrees") == -43
    assert field.metadata("longitudeOfSouthernPoleInDegrees") == 10
    assert field.shape == (19, 36)
    assert np.isfinite(field.to_numpy()).all()


@requires_mir
@requires_ekd
def test_regrid_grib_message_rotation_string():
    out_grid = {"grid": [1, 1], "area": [10, -10, -10, 10], "rotation": "-43/10"}
    field = _regrid_o32(out_grid)
    ref = _regrid_o32({**out_grid, "rotation": [-43, 10]})

    assert field.metadata("gridType") == "rotated_ll"
    assert field.metadata("latitudeOfSouthernPoleInDegrees") == -43
    assert field.metadata("longitudeOfSouthernPoleInDegrees") == 10
    np.testing.assert_array_equal(field.to_numpy(), ref.to_numpy())


@requires_mir
@requires_ekd
def test_regrid_grib_message_rotation_bytesio():
    with open(earthkit_test_data_path("o32.grib2"), "rb") as fh:
        data = BytesIO(fh.read())
    out_grid = {"grid": [1, 1], "area": [10, -10, -10, 10], "rotation": [-43, 10]}
    field = _regrid_o32(out_grid, data=data)

    assert field.metadata("gridType") == "rotated_ll"
    assert field.shape == (21, 21)


@requires_mir
@requires_ekd
@pytest.mark.parametrize("as_fieldlist", [False, True])
def test_regrid_grib_field_rotation(as_fieldlist):
    from earthkit.data import from_source

    fl = from_source("file", earthkit_test_data_path("o32.grib2")).to_fieldlist()
    data = fl if as_fieldlist else fl[0]
    out_grid = {"grid": [1, 1], "area": [10, -10, -10, 10], "rotation": [-43, 10]}
    with pytest.warns(DeprecationWarning):
        out = regrid(data, out_grid=out_grid, interpolation="nearest-neighbour")

    field = out[0] if as_fieldlist else out
    assert field.metadata("gridType") == "rotated_ll"
    assert field.metadata("latitudeOfSouthernPoleInDegrees") == -43
    assert field.metadata("longitudeOfSouthernPoleInDegrees") == 10
    assert field.shape == (21, 21)

    # same values as the GRIB-bytes path
    ref = _regrid_o32(out_grid)
    np.testing.assert_array_equal(field.to_numpy(), ref.to_numpy())
