# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import logging

from ..handler import DataHandler
from ..utils import create_grid_object

LOG = logging.getLogger(__name__)


class FieldListDataHandler(DataHandler):
    @staticmethod
    def match(data):
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.FieldList) and hasattr(data, "get")
        except Exception:
            return False

    def regrid(self, data, in_grid=None, out_grid=None, **kwargs):
        backend = self.backend_from_kwargs(kwargs)
        in_grid = create_grid_object(in_grid)
        out_grid = create_grid_object(out_grid)

        if out_grid is None:
            raise ValueError("out_grid must be specified")

        # in_grid, out_grid = self.parse_grid_kwargs(backend, in_grid, out_grid)

        fields = []
        for i, field in enumerate(data):
            r = self._regrid_field(field, i, backend, in_grid, out_grid, **kwargs)
            fields.append(r)

        from earthkit.data import FieldList

        return FieldList.from_fields(fields)

    def _regrid_field(self, field, index, backend, in_grid, out_grid, **kwargs):
        message = None
        if field._get_grib():
            message = field.message()

        if message is not None:
            from .grib import regrid_grib_field

            return regrid_grib_field(field, message, index, in_grid, out_grid, backend, **kwargs)
        else:
            from .generic import regrid_generic_field

            return regrid_generic_field(field, index, in_grid, out_grid, backend, **kwargs)


class FieldDataHandler(DataHandler):
    @staticmethod
    def match(data):
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.Field) and hasattr(data, "get")
        except Exception:
            return False

    def regrid(self, data, **kwargs):
        from earthkit.data import FieldList

        ds = FieldList.from_fields([data])
        return FieldListDataHandler().regrid(ds, **kwargs)[0]


handler = [FieldListDataHandler, FieldDataHandler]
