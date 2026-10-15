# (C) Copyright 2023 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#


"""Data-handler discovery and dispatch for the high-level ``regrid()`` entry point.

Imports the ``fieldlist``, ``grib`` and ``xarray`` submodules (each exposing
a ``handler`` attribute: a :class:`~.handler.DataHandler` subclass or a list
of them) into :data:`DATA_HANDLERS`, and provides :func:`get_data_handler`
to pick the handler whose :meth:`~.handler.DataHandler.match` accepts a
given data object.
"""

from importlib import import_module

_modules = [
    "fieldlist",
    "grib",
    "xarray",
]


DATA_HANDLERS = []
for name in _modules:
    module = import_module(f".{name}", package=__name__)
    lst = getattr(module, "handler", [])
    if not isinstance(lst, list):
        lst = [lst]
    assert isinstance(lst, list), f"Expected a list of handlers in {module.__name__}"
    DATA_HANDLERS.extend(lst)


def get_data_handler(values):
    """Return the registered data handler that can regrid ``values``.

    Parameters
    ----------
    values : Any
        The data object to find a handler for (e.g. an
        ``earthkit.data.FieldList``, an ``xarray.Dataset``, a GRIB message).

    Returns
    -------
    DataHandler or None
        A new instance of the first matching handler in
        :data:`DATA_HANDLERS`, or None if no handler matches.
    """
    for h in DATA_HANDLERS:
        if h.match(values):
            # TODO: rethink if we need to create a new handler instance each time
            return h()
