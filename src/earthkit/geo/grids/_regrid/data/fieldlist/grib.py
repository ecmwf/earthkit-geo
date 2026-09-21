# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import json
import logging

from .grid import get_grid

LOG = logging.getLogger(__name__)


def _encode_grib(message, v_res, out_grid_spec):
    # convert gridspec to string as this is what ecCodes expects in the "gridSpec" key when encoding GRIB messages
    out_spec = json.dumps(out_grid_spec)

    # when the field has an associated GRIB message, we
    # use a GribEncoder to encode the resulting data and preserve
    # the original GRIB metadata as much as possible (bar the grid)

    # TODO: avoid using a GribEncoder once the grid(spec) handling in the
    # earthkit-data field is improved
    from earthkit.data.encoders.grib import GribEncoder

    encoder = GribEncoder()
    d = encoder.encode(template=message, values=v_res, gridSpec=out_spec)
    res_message = d.to_bytes()

    return res_message


def regrid_grib_field(field, message, index, in_grid=None, out_grid=None, backend=None, **kwargs):
    from earthkit.data.field.grib.create import create_grib_field_from_message

    assert message is not None, "GRIB message must be provided"

    if in_grid is None:
        in_grid = get_grid(field, index)

    if in_grid is None:
        raise ValueError(
            f"Input grid could not be determined from data for field[{index}]."
            " Please provide an explicit input grid via the 'in_grid' parameter."
        )

    assert out_grid is not None, "Output grid must be provided"

    in_is_unstructured = in_grid.type in ("unstructured", "unstructured_ll")
    out_is_unstructured = out_grid.type in ("unstructured", "unstructured_ll")

    if hasattr(backend, "regrid_grib") and not in_is_unstructured and not out_is_unstructured:
        res_message = backend.regrid_grib(message, out_grid, **kwargs)
        res = create_grib_field_from_message(res_message, template_field=field)
    else:
        vv = field.to_numpy(flatten=True)
        v_res, out_grid_spec_res = backend.regrid(
            vv,
            in_grid,
            out_grid,
            **kwargs,
        )

        if out_is_unstructured:
            res = field.set({"values": v_res, "geography.grid_spec": out_grid_spec_res})
        else:
            res_message = _encode_grib(message, v_res, out_grid_spec_res)
            res = create_grib_field_from_message(res_message, template_field=field)

    return res
