# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import json
import logging

from .grid import get_grid, get_latlon_grid

LOG = logging.getLogger(__name__)

# The GRIB gridTypes whose specification describes a map projection.
PROJECTED_GRID_TYPES = ("lambert", "polar_stereographic", "mercator", "lambert_azimuthal_equal_area")


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

    # The grid of the result is the gridSpec the encoder is given: it resets the grid keys the
    # template brought along, so the message reads back as that grid.
    # TODO: writing the message from the template here, without the encoder, may be the simpler
    # way once it needs no more keys than the ones we hold.
    from earthkit.data.encoders.grib import GribEncoder

    encoder = GribEncoder()
    d = encoder.encode(template=message, values=v_res, gridSpec=out_spec)
    res_message = d.to_bytes()

    return res_message


def _is_spectral(field) -> bool:
    """Return True if ``field``'s GRIB grid is a spectral truncation (ecCodes gridType "sh").

    A spectral field stores spherical-harmonic coefficients, not grid values: it has
    no grid points, so its grid spec is None and its latitudes/longitudes cannot be
    computed. Any other grid type, and any field that does not answer the metadata
    key, returns False.
    """
    try:
        return field.metadata("gridType") == "sh"
    except Exception:
        return False


def _is_projected(field) -> bool:
    """Return True if ``field``'s GRIB grid is a map projection.

    A projected grid is one whose ecCodes gridType is in :data:`PROJECTED_GRID_TYPES`
    (Lambert conformal, polar stereographic, Mercator or Lambert azimuthal equal-area).
    Any other grid type, and any field that does not answer the metadata key, returns
    False.
    """
    try:
        return field.metadata("gridType") in PROJECTED_GRID_TYPES
    except Exception:
        return False


def regrid_grib_field(field, message, index, in_grid=None, out_grid=None, backend=None, **kwargs):
    """Regrid a GRIB-backed ``earthkit.data`` field.

    Uses the backend's ``regrid_grib`` method when available and applicable
    (structured input grid), which lets the backend regrid the GRIB message
    directly and preserve its metadata. Otherwise falls back to array-based
    regridding followed by re-encoding into a GRIB message.

    A spectral field (GRIB gridType "sh") has no grid points, so its input grid
    cannot be determined from latitudes/longitudes and the array path is not
    available at all. When no ``in_grid`` is given for such a field, the message
    is regridded through ``regrid_grib`` directly (which never uses the input
    grid), and a backend without that method is refused with a clear error.

    A projected field (GRIB gridType lambert, polar_stereographic, mercator or
    lambert_azimuthal_equal_area) goes the other way: when no ``in_grid`` is given,
    its input grid is built from its latitudes/longitudes, not from its grid spec,
    which keeps it on the array path while MIR's GRIB path cannot regrid these
    messages (the regrid-grib findings, §E).

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
        Input grid. When not provided it is determined from ``field``'s metadata,
        except for a spectral field, which goes straight to the GRIB path, and a
        projected field, whose grid comes from its latitudes/longitudes and which
        therefore goes straight to the array path.
    out_grid : eckit.geo.Grid or earthkit.geo.grids.Grid
        Output grid. Must be provided.
    backend : object
        Regrid backend exposing ``regrid(values, in_grid, out_grid, **kwargs)`` and,
        optionally, ``regrid_grib(message, out_grid, **kwargs)``. A spectral field
        requires the latter.
    **kwargs : dict
        Additional keyword arguments passed to the backend.

    Returns
    -------
    earthkit.data.Field
        A new field with the regridded values and updated grid spec.

    Raises
    ------
    ValueError
        If ``in_grid`` cannot be determined, ``out_grid`` is not provided, or the
        input is spectral and the backend has no ``regrid_grib`` method.
    """
    from earthkit.data.field.grib.create import create_grib_field_from_message

    assert message is not None, "GRIB message must be provided"

    # A spectral field has no grid points: get_grid() below would ask earthkit-data for the
    # latitudes/longitudes of a grid that has none and get (None, None) back, which is no input
    # grid at all. MIR regrids spectral GRIB messages natively, and regrid_grib() does not use
    # in_grid, so such a field goes straight to it.
    spectral = in_grid is None and _is_spectral(field)

    if spectral and not hasattr(backend, "regrid_grib"):
        raise ValueError(
            f"Spectral input field[{index}] (GRIB gridType 'sh') has no grid points: regridding it"
            " needs a backend that regrids the GRIB message directly, i.e. one with a"
            " 'regrid_grib' method, such as the MIR backend."
        )

    if in_grid is None and not spectral:
        if _is_projected(field):
            # A projected field takes the array path even when earthkit-data answers a grid
            # spec for it, so its input grid comes from its latitudes/longitudes: MIR 1.30.2
            # cannot regrid projected messages through its GRIB path (mercator "GRIB
            # packing=complex requires spectral data", polar "type missing in Params",
            # lambert a wrong output shape -- the regrid-grib findings, §E).
            # With the spec the in_grid below would be structured and the field would be
            # handed to regrid_grib and fail; the latlons keep it on the array path that
            # worked before earthkit-data started returning projected grid specs.
            # TODO: remove this routing once MIR's GRIB path supports projected grids, when
            # the plain get_grid() below takes them through the spec as it does every other
            # structured grid.
            in_grid = get_latlon_grid(field, index)
        else:
            in_grid = get_grid(field, index)

        if in_grid is None:
            raise ValueError(
                f"Input grid could not be determined from field[{index}]."
                " Please provide an explicit input grid via the 'in_grid' parameter."
            )

    assert out_grid is not None, "Output grid must be provided"

    # for a spectral field in_grid stays None: it takes the GRIB path below, which never reads it
    in_is_unstructured = in_grid is not None and in_grid.type in ("unstructured", "unstructured_ll")
    out_is_unstructured = out_grid.type in ("unstructured", "unstructured_ll")

    # for certain cases use grib-specific regridding via the backend
    if hasattr(backend, "regrid_grib"):
        if not in_is_unstructured:
            res_message = backend.regrid_grib(message, out_grid, **kwargs)
            res = create_grib_field_from_message(res_message, template_field=field)
            if out_is_unstructured:
                # The mir backend re-encodes unstructured output as GRIB2 (even for GRIB1 input).
                # GRIB cannot hold the points, so the grid spec must be set explicitly.
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
