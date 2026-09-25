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
    """Encode regridded values into a new GRIB message.

    Parameters
    ----------
    message : eccodes GRIB message
        Template GRIB message whose metadata (bar the grid) is preserved.
    v_res : numpy.ndarray
        Regridded values to encode.
    out_grid_spec : dict
        Grid spec of the output grid, encoded into the message's "gridSpec" key.

    Returns
    -------
    bytes
        The encoded GRIB message.
    """
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
    """Regrid a GRIB-backed ``earthkit.data`` field.

    Uses the backend's ``regrid_grib`` method when available and applicable
    (structured input grid), which lets the backend regrid the GRIB message
    directly and preserve its metadata. Otherwise falls back to array-based
    regridding followed by re-encoding into a GRIB message.

    Parameters
    ----------
    field : earthkit.data.Field
        The field to regrid.
    message : eccodes GRIB message
        GRIB message associated with ``field``. Must be provided.
    index : int
        Position of ``field`` in its parent :class:`~earthkit.data.FieldList`, used
        only for error/log messages.
    in_grid : eckit.geo.Grid or earthkit.geo.grids.Grid, optional
        Input grid. When not provided it is determined from ``field``'s metadata.
    out_grid : eckit.geo.Grid or earthkit.geo.grids.Grid
        Output grid. Must be provided.
    backend : object
        Regrid backend exposing ``regrid(values, in_grid, out_grid, **kwargs)`` and,
        optionally, ``regrid_grib(message, out_grid, **kwargs)``.
    **kwargs : dict
        Additional keyword arguments passed to the backend.

    Returns
    -------
    earthkit.data.Field
        A new field with the regridded values and updated grid spec.

    Raises
    ------
    ValueError
        If ``in_grid`` cannot be determined or ``out_grid`` is not provided.
    """
    from earthkit.data.field.grib.create import create_grib_field_from_message

    assert message is not None, "GRIB message must be provided"

    if in_grid is None:
        in_grid = get_grid(field, index)

    if in_grid is None:
        raise ValueError(
            f"Input grid could not be determined from field[{index}]."
            " Please provide an explicit input grid via the 'in_grid' parameter."
        )

    assert out_grid is not None, "Output grid must be provided"

    in_is_unstructured = in_grid.type in ("unstructured", "unstructured_ll")
    out_is_unstructured = out_grid.type in ("unstructured", "unstructured_ll")

    # for certain cases use grib-specific regridding via the backend
    if hasattr(backend, "regrid_grib"):
        if not in_is_unstructured:
            res_message = backend.regrid_grib(message, out_grid, **kwargs)
            if not out_is_unstructured:
                res = create_grib_field_from_message(res_message, template_field=field)
            else:
                # When the output grid is unstructured, the mir backend changes the grib encoding,
                # from edition 1 to edition 2.
                res = create_grib_field_from_message(res_message, template_field=field)
                res = res.set({"geography.grid_spec": out_grid})
            return res

    # otherwise use array-based regridding via the backend
    vv = field.to_numpy(flatten=True)
    v_res, out_grid_spec_res = backend.regrid(
        vv,
        in_grid,
        out_grid,
        **kwargs,
    )

    if out_is_unstructured:
        # When the output grid is unstructured, at present we cannot encode it back into
        # a standard GRIB message, so we modify the field to store the values and grid spec
        # directly. With this the original field metadata is kept.
        res = field.set({"values": v_res, "geography.grid_spec": out_grid_spec_res})
    else:
        res_message = _encode_grib(message, v_res, out_grid_spec_res)
        res = create_grib_field_from_message(res_message, template_field=field)

    return res
