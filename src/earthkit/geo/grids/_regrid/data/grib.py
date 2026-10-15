# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Raw GRIB message backend for the regrid data-handler dispatch mechanism.

Exposes :data:`handler`, the :class:`GribMessageDataHandler` class used to
regrid a raw GRIB message given as ``bytes`` or a ``BytesIO``.
"""

import logging

from earthkit.geo.grids._regrid.data.handler import DataHandler
from earthkit.geo.utils import is_module_loaded

LOG = logging.getLogger(__name__)


class GribMessageDataHandler(DataHandler):
    """Regrid data handler for a raw GRIB message (``bytes`` or ``BytesIO``)."""

    @staticmethod
    def match(values):
        """Return True if ``values`` is a ``bytes`` or ``BytesIO`` object.

        Parameters
        ----------
        values : Any
            The data to check.

        Returns
        -------
        bool
            True if ``values`` looks like a GRIB message.
        """
        if isinstance(values, bytes):
            return True
        else:
            from io import BytesIO

            # TODO: add further checks to see if the object is a GRIB message
            if isinstance(values, BytesIO):
                return True
        return False

    def regrid(self, values, in_grid=None, out_grid=None, **kwargs):
        """Regrid a raw GRIB message.

        If the selected backend implements ``regrid_grib``, the message is passed to it directly and
        the input grid is determined by the backend from the message itself (``in_grid`` is unused in
        this case). Otherwise, when ``earthkit.data`` is available, the message is wrapped into an
        earthkit-data ``Field`` and regridded via :class:`~.fieldlist.FieldDataHandler` (which does use
        ``in_grid`` when the grid cannot be inferred from the message), then re-encoded back into a GRIB
        message. This fallback also brings the unstructured-grid handling of the ``Field``/``FieldList``
        path (see :mod:`~.fieldlist.grib`) to raw GRIB message input.

        Parameters
        ----------
        values : bytes or io.BytesIO
            The GRIB message to regrid.
        in_grid : Any, optional
            The input grid or grid spec. Only used as a fallback, when the backend has no
            ``regrid_grib`` method and the input grid cannot be inferred from the message.
        out_grid : Any
            The output grid or grid spec. Must be provided.
        **kwargs : dict
            Must include ``backend`` (name or object of the regrid backend
            to use). Any other keyword arguments are forwarded to the
            backend's ``regrid_grib`` method (or to the ``Field`` regridding
            fallback).

        Returns
        -------
        bytes
            The regridded GRIB message.

        Raises
        ------
        ValueError
            If the selected backend does not support GRIB message input
            (i.e. has no ``regrid_grib`` method) and ``earthkit.data`` is not available to fall back to.
        """
        backend = self.backend_from_kwargs(kwargs)
        if not isinstance(values, bytes):
            from io import BytesIO

            if isinstance(values, BytesIO):
                values = values.getvalue()

        if hasattr(backend, "regrid_grib"):
            return backend.regrid_grib(values, out_grid, **kwargs)
        else:
            if is_module_loaded("earthkit.data"):
                from earthkit.data.field.grib.create import create_grib_field_from_message

                field = create_grib_field_from_message(values)
                from .fieldlist import FieldDataHandler

                res_field = FieldDataHandler().regrid(
                    field, in_grid=in_grid, out_grid=out_grid, backend=backend, **kwargs
                )
                return res_field.sync().message()

            raise ValueError(
                f"regrid() does not support GRIB message input for {backend=} when no earthkit.data is available!"
            )


handler = GribMessageDataHandler
