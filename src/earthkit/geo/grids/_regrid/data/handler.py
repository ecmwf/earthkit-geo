# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#

"""Abstract base class for the regrid data-handler dispatch mechanism.

A :class:`DataHandler` subclass regrids one supported data type (e.g. a
``FieldList``, an xarray object, a raw GRIB message) and is registered via a
``handler`` attribute in its module (see
:mod:`earthkit.geo.grids._regrid.data`).
"""

import logging
from abc import ABCMeta, abstractmethod

LOG = logging.getLogger(__name__)


OPTIONAL_BACKENDS_KWARGS = ["inventory"]


class DataHandler(metaclass=ABCMeta):
    """Abstract base class for a data-type-specific regrid handler.

    Subclasses must implement :meth:`regrid` and are expected to expose a
    ``match(values)`` static method used by
    :func:`~earthkit.geo.grids._regrid.data.get_data_handler` to select the
    handler for a given data object.
    """

    @abstractmethod
    def regrid(self, values, in_grid=None, out_grid=None, **kwargs):
        """Regrid ``values`` from ``in_grid`` onto ``out_grid``.

        Parameters
        ----------
        values : Any
            The data to regrid, of the type this handler supports.
        in_grid : Any, optional
            The input grid or grid spec.
        out_grid : Any, optional
            The output grid or grid spec.
        **kwargs : dict
            Must include ``backend`` (name or object of the regrid backend
            to use, see :meth:`backend_from_kwargs`). Any other keyword
            arguments are handler/backend-specific.

        Returns
        -------
        Any
            The regridded data, of the same type as ``values``.
        """
        pass

    def backend_from_kwargs(self, kwargs):
        """Pop and resolve the ``backend`` (and its options) from ``kwargs``.

        Parameters
        ----------
        kwargs : dict
            The keyword arguments to pop from, modified in place. The
            ``"backend"`` key is required; any keys listed in
            :data:`OPTIONAL_BACKENDS_KWARGS` (e.g. ``"inventory"``) are also
            popped and forwarded to the backend constructor.

        Returns
        -------
        Backend
            The resolved backend instance.

        Raises
        ------
        ValueError
            If ``kwargs`` has no ``"backend"`` key.
        """
        backend = kwargs.pop("backend", None)
        if backend is None:
            raise ValueError("Missing 'backend' keyword argument")

        from ..backends import Backend

        if isinstance(backend, Backend):
            return backend

        b_kwargs = {}
        for k in OPTIONAL_BACKENDS_KWARGS:
            if k in kwargs:
                b_kwargs[k] = kwargs.pop(k)
        return self.get_backend(backend, **b_kwargs)

    def get_backend(self, backend, **backend_kwargs):
        """Resolve a backend name (or object) into a backend instance.

        Parameters
        ----------
        backend : str or Backend
            The registered backend name, or an already-constructed backend.
        **backend_kwargs : dict
            Keyword arguments passed to the backend constructor.

        Returns
        -------
        Backend
            The resolved backend instance.

        Raises
        ------
        ValueError
            If no backend can be resolved.
        """
        from earthkit.geo.grids._regrid.backends import get_backend

        backend = get_backend(backend, **backend_kwargs)

        if not backend:
            raise ValueError(f"No backend={backend} found")

        return backend
