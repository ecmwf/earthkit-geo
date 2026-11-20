# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


import logging

from ..utils import create_grid_object

LOG = logging.getLogger(__name__)


def get_latlon_grid(field, index):
    """Determine the grid of an ``earthkit.data`` field from its own latitudes/longitudes.

    Builds the unstructured grid that holds the field's points exactly as ecCodes
    placed them on the sphere, whatever projection the message is written on. This
    is what :func:`get_grid` falls back to when the field has no grid spec, and
    what callers use to force the latlons even for a field that has one (see
    ``regrid_grib_field``'s routing of projected grids).

    Parameters
    ----------
    field : earthkit.data.Field
        The field to inspect.
    index : int
        Position of ``field`` in its parent :class:`~earthkit.data.FieldList`, used
        only for log messages.

    Returns
    -------
    eckit.geo.Grid or None
        The field's latlons grid, or None if it has no latitudes/longitudes.

    Raises
    ------
    Exception
        Whatever ``field.geography.latlons`` raises when the points cannot be
        computed (re-logged before being re-raised).
    """
    try:
        lat, lon = field.geography.latlons(flatten=True)
        if lat is None or lon is None:
            return None
        grid_spec = {"latitudes": lat.tolist(), "longitudes": lon.tolist()}
    except Exception as e:
        LOG.exception(f"Cannot get latlons for field[{index}]: {e}")
        raise

    return create_grid_object(grid_spec)


def get_grid(field, index):
    """Determine the grid of an ``earthkit.data`` field.

    Tries the field's grid spec first, falling back to its latitude/longitude
    values when no grid spec is available (see :func:`get_latlon_grid`).

    Parameters
    ----------
    field : earthkit.data.Field
        The field to inspect.
    index : int
        Position of ``field`` in its parent :class:`~earthkit.data.FieldList`, used
        only for log messages.

    Returns
    -------
    eckit.geo.Grid or None
        The field's grid, or None if it could not be determined.
    """
    try:
        grid_spec = field.geography.grid_spec()
    except Exception as e:
        LOG.warning(f"Cannot get input grid_spec from metadata for field[{index}]: {e}")
        grid_spec = None

    if grid_spec is None:
        return get_latlon_grid(field, index)

    if not grid_spec:
        return None

    return create_grid_object(grid_spec)
