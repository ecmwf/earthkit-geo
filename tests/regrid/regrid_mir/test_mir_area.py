# (C) Copyright 2026 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

"""Regression test for an ``area`` inside the output grid spec of the MIR backend.

The backend once moved ``area`` out of the grid spec into the job options and
warned about it (``MirBackend.adjust_options``, a ``DeprecationWarning``). That
workaround is gone: the spec reaches MIR as it is, so regridding onto the
regional 1 x 1 degree target of the matrix case ``ll-1x1-area``
(``regrid_cases.py``) must return the 48 x 71 grid with no warning and leave
the dict it was given unchanged -- on both MIR paths, the array and the GRIB
message one. The input is the small global 10 x 10 degree lat-lon grid and the
tracked ``tests/data/o32.grib2``, never downloaded (see the skip marks, which
come from ``regrid_common`` as elsewhere in this folder).
"""

import warnings

import numpy as np
import pytest
from regrid_common import NO_EKD, requires_mir

from earthkit.geo.grids._regrid.backends.mir import MirBackend
from earthkit.geo.utils.testing import earthkit_test_data_path

try:
    import eccodes
except ImportError:
    eccodes = None

# the regional target of the matrix case ll-1x1-area: lats 72..25 and lons
# -25..45 in 1 degree steps, ie (48, 71) points
OUT_SPEC = {"grid": [1, 1], "area": [72, -25, 25, 45]}
OUT_SHAPE = (48, 71)


@requires_mir
def test_regrid_array_area_in_output_spec():
    """An ``area`` in the array regridding's output spec passes to MIR unchanged."""
    # the global 19 x 36 point 10 x 10 degree input grid
    values = np.random.rand(19 * 36)

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        result, _ = MirBackend().regrid(values, {"grid": [10, 10]}, OUT_SPEC)

    assert result.shape == OUT_SHAPE
    # the backend must not touch the dict it was given (no popping the area out)
    assert OUT_SPEC == {"grid": [1, 1], "area": [72, -25, 25, 45]}


@requires_mir
@pytest.mark.skipif(
    eccodes is None and NO_EKD,
    reason="Decoding the output GRIB needs eccodes or earthkit-data",
)
def test_regrid_grib_area_in_output_spec():
    """An ``area`` in the GRIB regridding's output spec passes to MIR unchanged."""
    with open(earthkit_test_data_path("o32.grib2"), "rb") as fh:
        in_grib = fh.read()

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        out = MirBackend().regrid_grib(in_grib, OUT_SPEC)

    assert out, "regrid_grib returned no GRIB bytes"
    assert out.startswith(b"GRIB") and out.endswith(b"7777")

    if eccodes is not None:
        handle = eccodes.codes_new_from_message(out)
        try:
            assert eccodes.codes_get(handle, "gridType") == "regular_ll"
            assert (eccodes.codes_get(handle, "Nj"), eccodes.codes_get(handle, "Ni")) == OUT_SHAPE
            assert eccodes.codes_get(handle, "numberOfDataPoints") == OUT_SHAPE[0] * OUT_SHAPE[1]
        finally:
            eccodes.codes_release(handle)
    else:
        from earthkit.data.field.grib.create import create_grib_field_from_buffer

        field = create_grib_field_from_buffer(out)
        assert field.shape == OUT_SHAPE

    assert OUT_SPEC == {"grid": [1, 1], "area": [72, -25, 25, 45]}
