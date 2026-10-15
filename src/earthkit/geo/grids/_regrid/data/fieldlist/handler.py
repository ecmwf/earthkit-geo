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
    """Regrid data handler for ``earthkit.data.FieldList`` objects."""

    @staticmethod
    def match(data):
        """Return True if ``data`` is an ``earthkit.data.FieldList``."""
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.FieldList) and hasattr(data, "get")
        except Exception:
            return False

    def regrid(self, data, in_grid=None, out_grid=None, **kwargs):
        """Regrid every field in a ``FieldList``.

        Parameters
        ----------
        data : earthkit.data.FieldList
            The fields to regrid.
        in_grid : Any, optional
            Input grid or grid spec. When not provided it is determined per-field
            from its metadata.
        out_grid : Any
            Output grid or grid spec. Must be provided.
        **kwargs : dict
            Must include ``backend`` (name or object of the regrid backend to use).
            Any other keyword arguments are forwarded to the backend.

        Returns
        -------
        earthkit.data.FieldList
            A new field list with the regridded fields.

        Raises
        ------
        ValueError
            If ``out_grid`` is not provided.
        """
        backend = self.backend_from_kwargs(kwargs)
        in_grid = create_grid_object(in_grid)
        out_grid = create_grid_object(out_grid)

        if out_grid is None:
            raise ValueError("out_grid must be specified")

        fields = []
        for i, field in enumerate(data):
            r = self._regrid_field(field, i, backend, in_grid, out_grid, **kwargs)
            fields.append(r)

        from earthkit.data import FieldList

        return FieldList.from_fields(fields)

    def _regrid_field(self, field, index, backend, in_grid, out_grid, **kwargs):
        """Regrid a single field, dispatching to the GRIB or generic implementation.

        Parameters
        ----------
        field : earthkit.data.Field
            The field to regrid.
        index : int
            Position of ``field`` in its parent :class:`~earthkit.data.FieldList`.
        backend : object
            Regrid backend to use.
        in_grid : Any, optional
            Input grid or grid spec.
        out_grid : Any
            Output grid or grid spec.
        **kwargs : dict
            Additional keyword arguments forwarded to the backend.

        Returns
        -------
        earthkit.data.Field
            The regridded field.
        """
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
    """Regrid data handler for a single ``earthkit.data.Field``."""

    @staticmethod
    def match(data):
        """Return True if ``data`` is an ``earthkit.data.Field``."""
        from earthkit.geo.utils import is_module_loaded

        if not is_module_loaded("earthkit.data"):
            return False

        try:
            import earthkit.data

            return isinstance(data, earthkit.data.Field) and hasattr(data, "get")
        except Exception:
            return False

    def regrid(self, data, **kwargs):
        """Regrid a single field by wrapping it in a one-element ``FieldList``.

        Parameters
        ----------
        data : earthkit.data.Field
            The field to regrid.
        **kwargs : dict
            Forwarded to :meth:`FieldListDataHandler.regrid`.

        Returns
        -------
        earthkit.data.Field
            The regridded field.
        """
        from earthkit.data import FieldList

        ds = FieldList.from_fields([data])
        return FieldListDataHandler().regrid(ds, **kwargs)[0]


handler = [FieldListDataHandler, FieldDataHandler]
