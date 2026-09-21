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


def get_grid(field, index):
    try:
        grid_spec = field.geography.grid_spec()
    except Exception as e:
        LOG.warning(f"Cannot get input grid_spec from metadata for field[{index}]: {e}")
        grid_spec = None

    if grid_spec is None:
        try:
            lat, lon = field.geography.latlon()
            if lat is not None and lon is not None:
                grid_spec = {"latitudes": lat.tolist(), "longitudes": lon.tolist()}
        except Exception as e:
            LOG.exception(f"Cannot get latlons for field[{index}]: {e}")
            raise

    if not grid_spec:
        return None

    return create_grid_object(grid_spec)
