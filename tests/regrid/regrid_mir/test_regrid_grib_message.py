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

from io import BytesIO

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
